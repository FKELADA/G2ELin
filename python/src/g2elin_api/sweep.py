"""Parameter sweeps for root-locus plots: re-solve the power flow and the
closed-loop eigenvalues of a network while one of its parameters steps
across a range.

Each step is a full rebuild -- the parameter can be anything that shapes the
operating point (a load, a line impedance, a setpoint) as well as a pure
control gain -- so a sweep costs one power flow plus one linearisation per
value. Results stream as newline-delimited JSON, one line per value, so the
client can draw the loci while the sweep runs.

Eigenvalues are reordered from one step to the next by a minimum-cost
assignment (Hungarian algorithm on a scale-aware distance), so index ``k``
refers to the same mode at every step and the client can join the points of
a mode into one locus.
"""

from __future__ import annotations

import json
import math
from typing import AsyncIterator

import numpy as np
from fastapi import HTTPException, Request
from scipy.optimize import linear_sum_assignment
from starlette.concurrency import run_in_threadpool

from g2elin_core.modal import ModalAnalysisResult, analyze
from g2elin_core.network.breakers import energized_network
from g2elin_core.network.schema import Network
from g2elin_core import tuning
from g2elin_core.operating_point import overridable_param_keys, unit_params
from g2elin_core.pipeline import linearize_network
from g2elin_core.powerflow import run_power_flow

from .schemas import SweepRequest, SweepTarget

# How many states of each mode a sweep streams (the page shows this many bars).
PARTICIPATION_TOP = 15

SWEEP_MAX_POINTS = 201

# Numeric fields a sweep may vary, per element kind (ids and bus references
# are topology, not parameters).
SWEEPABLE_FIELDS = {
    "network": {"f_hz", "sn_mva"},
    "bus": {"vn_kv"},
    "line": {"r_pu", "x_pu", "b_pu", "length_km"},
    "transformer": {"r_pu", "x_pu", "sn_mva"},
    "load": {"p_mw", "q_mvar"},
    "unit": {"v_set_pu", "p_set_mw", "q_set_mvar", "p_cons_mw", "q_cons_mvar", "xd_pu"},
}
_LIST_KEYS = {"line": "lines", "transformer": "transformers", "load": "loads"}


def sweep_values(req: SweepRequest) -> list[float]:
    """``start`` to ``stop`` inclusive in increments of ``step`` (its sign is
    taken from the direction of the range)."""
    if not (math.isfinite(req.start) and math.isfinite(req.stop) and math.isfinite(req.step)) or req.step == 0:
        raise HTTPException(status_code=422, detail="start, stop and step must be finite, and step non-zero")
    span = req.stop - req.start
    step = math.copysign(abs(req.step), span) if span else abs(req.step)
    n = int(math.floor(abs(span) / abs(step) + 1e-9)) + 1
    if n > SWEEP_MAX_POINTS:
        raise HTTPException(status_code=422, detail=f"{n} values requested; at most {SWEEP_MAX_POINTS} -- use a larger step")
    values = [req.start + i * step for i in range(n)]
    if abs(values[-1] - req.stop) > 1e-9 * max(1.0, abs(req.stop)):
        values.append(req.stop)  # the range end is always included
    return values


def _locate(data: dict, element: str, key: int | None) -> dict:
    if element == "network":
        return data
    if element == "bus":
        found = next((b for b in data["buses"] if b["id"] == key), None)
    elif element == "unit":
        found = next((d for d in data["der_units"] if d["id"] == key), None)
    elif element in _LIST_KEYS:
        items = data[_LIST_KEYS[element]]
        found = items[key] if key is not None and 0 <= key < len(items) else None
    else:
        raise HTTPException(status_code=422, detail=f"unknown element {element!r}; one of {sorted(SWEEPABLE_FIELDS)}")
    if found is None:
        raise HTTPException(status_code=422, detail=f"no {element} {key!r} in this network")
    return found


def _check_one(network: Network, t: SweepTarget) -> None:
    el = _locate(network.model_dump(), t.element, t.key)
    if t.field.startswith("tune."):
        if t.element != "unit":
            raise HTTPException(status_code=422, detail="tune.* fields only exist on units")
        _tuning_target(el, t.field)
    elif t.field.startswith("params."):
        if t.element != "unit":
            raise HTTPException(status_code=422, detail="params.* fields only exist on units")
        name = t.field.split(".", 1)[1]
        if name not in overridable_param_keys(
            el["unit_type"], exciter=el.get("exciter"), pss=el.get("pss"),
            governor=el.get("governor"),
        ):
            raise HTTPException(status_code=422, detail=f"a {el['unit_type']} unit has no parameter {name!r}")
    elif t.field not in SWEEPABLE_FIELDS.get(t.element, set()):
        raise HTTPException(
            status_code=422,
            detail=f"{t.element}.{t.field} can't be swept; one of {sorted(SWEEPABLE_FIELDS[t.element])} (or params.* on a unit)",
        )


def check_target(network: Network, req: SweepRequest) -> None:
    """422 on an unknown element/field -- for the main target and every
    extra one -- or on the same parameter listed twice, before any step runs."""
    targets = [req.target] + [x.target for x in req.extra]
    for t in targets:
        _check_one(network, t)
    keys = [(t.element, t.key, t.field) for t in targets]
    if len(set(keys)) != len(keys):
        raise HTTPException(status_code=422, detail="the same parameter is listed more than once")
    for x in req.extra:
        if not (math.isfinite(x.start) and math.isfinite(x.stop)):
            raise HTTPException(status_code=422, detail="start and stop must be finite")


def extra_values(req: SweepRequest, values: list[float]) -> list[list[float]]:
    """Each extra parameter's value at every step: it moves from its start to
    its stop in proportion to the main parameter's progress, so all of them
    reach their end together."""
    span = values[-1] - values[0]
    frac = [(v - values[0]) / span if span else 0.0 for v in values]
    return [[x.start + (x.stop - x.start) * f for f in frac] for x in req.extra]


def _tuning_target(el: dict, field: str) -> tuple[tuning.Loop, str, float]:
    """``tune.<loop>.<quantity>`` -> (loop, target keyword, SI scale); 422 if
    the unit has no such loop or quantity."""
    parts = field.split(".")
    if len(parts) != 3:
        raise HTTPException(status_code=422, detail="tuning fields look like tune.<loop>.<quantity>, e.g. tune.cl.tr_ms")
    _, loop_id, quantity = parts
    try:
        lp = tuning.loop(el["unit_type"], loop_id)
    except KeyError as e:
        raise HTTPException(status_code=422, detail=str(e.args[0])) from None
    if quantity not in tuning.quantities_for(lp):
        raise HTTPException(
            status_code=422, detail=f"{lp.name} can be swept by {tuning.quantities_for(lp)}, not {quantity!r}"
        )
    keyword, scale = tuning.TUNING_QUANTITIES[quantity]
    return lp, keyword, scale


def _set(data: dict, network: Network, target: SweepTarget, value: float) -> None:
    """Sets one swept parameter in ``data`` (a ``Network.model_dump()``);
    ``network`` is the state the value applies to (tuning targets derive
    their gains from its parameters)."""
    el = _locate(data, target.element, target.key)
    if target.field.startswith("tune."):
        lp, keyword, scale = _tuning_target(el, target.field)
        der = next(d for d in network.der_units if d.id == el["id"])
        el.setdefault("params", {}).update(tuning.gains_for(lp, unit_params(network, der), **{keyword: value * scale}))
    elif target.field.startswith("params."):
        el.setdefault("params", {})[target.field.split(".", 1)[1]] = value
    elif target.element == "line" and target.field == "length_km":
        # Line impedances are stored as totals; a longer line at the same
        # per-km constants has proportionally larger R, X and B.
        if value <= 0:
            raise ValueError("length must be positive")
        scale = value / el["length_km"]
        for k in ("r_pu", "x_pu", "b_pu"):
            el[k] *= scale
        el["length_km"] = value
    else:
        el[target.field] = value


def network_with(network: Network, req: SweepRequest, value: float, extras: list[float] = ()) -> Network:
    """A copy of ``network`` with the swept parameter(s) set -- the main one
    to ``value``, each extra one to its entry in ``extras`` -- re-validated,
    so e.g. a non-positive reactance fails that step.

    A ``tune.<loop>.<quantity>`` target (a loop's response time in ms,
    damping or inertia) becomes that loop's controller gains, the loop's
    other quantity held at its current value. Tuning targets are applied
    after every other one, so they see the other swept values (e.g. a swept
    filter inductance) when deriving the gains."""
    pairs = [(req.target, value)] + [(x.target, v) for x, v in zip(req.extra, extras)]
    pairs.sort(key=lambda tv: tv[0].field.startswith("tune."))
    data = network.model_dump()
    current = network
    for target, v in pairs:
        if target.field.startswith("tune."):
            current = Network(**data)
        _set(data, current, target, v)
    return Network(**data)


def modal_at(network: Network) -> "ModalAnalysisResult":
    """The linearised model's modal analysis at one point of a sweep."""
    network = energized_network(network)
    result = run_power_flow(network)
    if not result.converged:
        raise ValueError("power flow did not converge")
    system = linearize_network(network, result)
    return analyze(system.A, system.state_names)


def eigenvalues_at(network: Network) -> np.ndarray:
    return modal_at(network).eigenvalues


def match_order(previous: np.ndarray, current: np.ndarray) -> np.ndarray:
    """The permutation putting ``current`` in ``previous``'s order, i.e. so
    that ``current[order][k]`` continues ``previous[k]``'s locus. Distance is
    relative to magnitude, since these eigenvalues span many orders of
    magnitude (a fixed absolute distance would pair every slow mode with its
    nearest neighbour and ignore the fast ones entirely).
    """
    if len(previous) != len(current):
        return np.arange(len(current))
    a, b = previous[:, None], current[None, :]
    cost = np.abs(a - b) / (1.0 + np.minimum(np.abs(a), np.abs(b)))
    _, cols = linear_sum_assignment(cost)
    return cols


def match_to(previous: np.ndarray, current: np.ndarray) -> np.ndarray:
    """``current`` reordered so ``current[k]`` continues ``previous[k]``'s locus."""
    return current[match_order(previous, current)]


def top_participation(participation: np.ndarray, top: int, floor: float = 0.005) -> list[list[list[float]]]:
    """Per mode (column), its ``top`` most participating states as
    ``[state index, factor]`` pairs, dropping the negligible ones -- the whole
    matrix is (states x modes) and far too big to stream at every step, while
    what a reader looks at is the handful of states that drive a mode.
    """
    out: list[list[list[float]]] = []
    for j in range(participation.shape[1]):
        col = participation[:, j]
        idx = np.argsort(col)[::-1][:top]
        out.append([[int(i), round(float(col[i]), 4)] for i in idx if col[i] >= floor])
    return out


async def sweep_stream(network: Network, req: SweepRequest, request: Request) -> AsyncIterator[str]:
    """NDJSON lines: ``{"i", "value", "ok": true, "eig": [[re, im], ...]}`` per
    value (``"ok": false, "error"`` for a value that doesn't solve), then
    ``{"done": true, "n": ...}``.

    With ``req.participation`` (the default), a solved line also carries
    ``"part"``: each mode's most participating states at that value, in the
    same order as ``"eig"``, so the page can show how a mode's composition
    moves along its locus. The state names they index come with the first
    solved line (``"state_names"``).
    """
    values = sweep_values(req)
    extras = extra_values(req, values)
    yield json.dumps({"start": True, "values": values, "extra_values": extras}) + "\n"
    previous: np.ndarray | None = None
    sent_names = False
    for i, value in enumerate(values):
        if await request.is_disconnected():
            return
        step_extras = [xs[i] for xs in extras]
        try:
            modal = await run_in_threadpool(lambda v=value, xs=step_extras: modal_at(network_with(network, req, v, xs)))
        except HTTPException as e:
            yield json.dumps({"i": i, "value": value, "ok": False, "error": str(e.detail)}) + "\n"
            continue
        except Exception as e:  # noqa: BLE001 -- one bad value shouldn't end the sweep
            yield json.dumps({"i": i, "value": value, "ok": False, "error": f"{type(e).__name__}: {e}"}) + "\n"
            continue
        eig = modal.eigenvalues
        order = np.lexsort((eig.imag, eig.real)) if previous is None else match_order(previous, eig)
        eig = eig[order]
        previous = eig
        line = {"i": i, "value": value, "ok": True, "eig": [[float(z.real), float(z.imag)] for z in eig]}
        if req.participation:
            if not sent_names:
                line["state_names"] = list(modal.state_names)
                sent_names = True
            line["part"] = top_participation(modal.participation[:, order], PARTICIPATION_TOP)
        yield json.dumps(line) + "\n"
    yield json.dumps({"done": True, "n": len(values)}) + "\n"
