"""Network-based analysis cores, shared by the preset-id-keyed endpoints in
``main.py`` and the arbitrary-``Network``-JSON endpoints in
``network_routes.py`` -- both are "run this analysis on this ``Network``,"
differing only in how the ``Network`` is obtained (a fixed preset vs. an
arbitrary POSTed body). Every function here takes an already-built
``Network`` and returns the exact response model its endpoint returns, so
neither caller re-implements request validation or error handling.
"""

from __future__ import annotations

import json
import math
from collections import OrderedDict
from dataclasses import asdict
from typing import AsyncIterator, Callable

import numpy as np
import scipy.signal
from fastapi import HTTPException, Request

from g2elin_core import reduction
from g2elin_core.interconnect import AssembledSystem
from g2elin_core.modal import (
    unregulated_frequency_mode,
    ModalAnalysisResult, analyze, check_adequacy, classify_modes, eigenvalue_sensitivity,
    free_response, mode_shape, parameter_sensitivity, reference_angle_modes, step_response,
)
from g2elin_core.network.breakers import TYPE_LABEL, NoReferenceUnit, energized_network, service_state
from g2elin_core.network.schema import Network
from g2elin_core.network.topology import compute_topology_layout
from g2elin_core.network.validation import validate_network
from g2elin_core.pipeline import linearize_network
from g2elin_core.powerflow import run_power_flow
from g2elin_core.powerflow.pandapower_adapter import PowerFlowResult
from g2elin_core.timedomain import (
    DEFAULT_SOLVER, SOLVERS, NonlinearNetworkModel, build_nonlinear_network, find_state_index,
    simulate, simulate_fixed_step, simulate_steps,
)
from g2elin_core.timedomain.events import EventError, NetworkEvent, apply_event
from g2elin_core.timedomain.measurements import MeasurementSet
from g2elin_core.timeseries import Snapshot, apply_snapshot, run_time_series, scale_loads

from .schemas import (
    AdequacyResponse,
    BatchPowerFlowRequest,
    BatchPowerFlowResponse,
    BatchSnapshot,
    BusNodeRow,
    BusRow,
    EmtRequest,
    EmtResponse,
    FreeResponseRequest,
    FreeResponseResponse,
    ElementModelInfo,
    LinearOverlay,
    MeasurementInfo,
    ModalResponse,
    ModeCategoryInfo,
    ModelLevelsResponse,
    ModelSummaryResponse,
    ModePairRow,
    ModeRow,
    ModeShapeRequest,
    ModeShapeResponse,
    EntryParametersRow,
    ParameterEffectRow,
    ParameterSensitivityRequest,
    ParameterSensitivityResponse,
    NetworkIssueRow,
    PowerFlowOptions,
    PowerFlowResponse,
    IslandRow,
    SensitivityEntryRow,
    SensitivityRequest,
    SensitivityResponse,
    SolverInfo,
    SolversResponse,
    StateRiskRow,
    ServiceInfo,
    StatesResponse,
    StepResponseRequest,
    StepResponseResponse,
    TimeSeriesResponse,
    TimeSeriesSnapshot,
    TopologyEdgeRow,
    TopologyResponse,
    UnitModelRow,
    ValidateResponse,
)

# Bounds on user-controlled simulation cost, enforced server-side rather than
# just documented -- these endpoints run a real nonlinear DAE integration
# synchronously inside one HTTP request.
#
# `t_final` deliberately has no ceiling: how long a run is worth watching is
# the user's call, not this file's, and a governor or an inter-area mode
# takes tens of seconds to say anything. What is still bounded is the part
# that is genuinely a *cost*: the number of samples returned, and, for a live
# run, the number of solver steps streamed. A long run is a long request, so
# the live endpoint (which streams as it goes) is the better way to watch
# one.
EMT_MIN_N_POINTS = 10
# One-shot runs: samples per trajectory. Each extra one costs a Newton solve
# when inputs/outputs/measurements are plotted, so this is a cost bound.
EMT_MAX_N_POINTS = 10000
# Live runs: solver steps streamed to the browser. Not a physical limit -- how
# many an adaptive stiff solver takes depends on the network, and a stiff one
# (CIGRE) takes far more than a small one -- just a bound on what one stream
# may send. Hitting it ends the run cleanly with what was traced
# (see emt_live_stream), it never fails it.
EMT_LIVE_MAX_STEPS = 50000
EMT_DEFAULT_N_POINTS = 200  # -> default dt = t_final/199 when req.dt isn't given
MODAL_MAX_T_FINAL = 20.0  # free/step response are cheap (linear algebra, no ODE solve) -- a looser bound
# The solver settings a request doesn't override. ``first_step`` is small on
# purpose: the full EMT model's fastest modes are ~1e7 rad/s, and an adaptive
# solver that opens with too big a step spends its first dozen steps being
# rejected. With the analytic Jacobian in place it matters much less than it
# did, but it costs nothing to start carefully.
EMT_FIRST_STEP = 1e-8
EMT_SIMULATE_KWARGS = dict(rtol=1e-4, atol=1e-6, first_step=EMT_FIRST_STEP)
# What a fixed-step run will not do: one that would take more steps than this
# is refused up front rather than left to time the request out.
EMT_MAX_FIXED_STEPS = 200_000
EMT_MAX_T_PRE = 1.0
BATCH_MAX_STEPS = 100
BATCH_MAX_SCALE = 5.0

# pandapower.runpp's ``algorithm`` values this API exposes -- id -> label
# (the label is what the web UI shows in its solver picker).
POWERFLOW_ALGORITHMS = {
    "nr": "Newton-Raphson",
    "iwamoto_nr": "Newton-Raphson with Iwamoto multiplier",
    "fdbx": "Fast-decoupled (BX)",
    "fdxb": "Fast-decoupled (XB)",
    "gs": "Gauss-Seidel",
    "bfsw": "Backward/forward sweep (radial networks)",
}
POWERFLOW_INITS = ("auto", "flat", "dc")


# --- Model cache -------------------------------------------------------------
# Every modal view (sensitivity, mode shape, free/step response...) and every
# EMT run re-derives the same model from the same Network -- a power flow plus
# a full linearization/nonlinear build each time. The web UI fires several of
# these per page, so memoize the last few models per network content. Keyed by
# the Network's own JSON: an edited network is a different key, never a stale
# hit. Both model types are only ever read after construction.
_MODEL_CACHE_SIZE = 8
_model_cache: "OrderedDict[tuple[str, str], object]" = OrderedDict()


def _cached(kind: str, network: Network, build: Callable[[], object]) -> object:
    key = (kind, network.model_dump_json())
    if key in _model_cache:
        _model_cache.move_to_end(key)
        return _model_cache[key]
    value = build()
    _model_cache[key] = value
    if len(_model_cache) > _MODEL_CACHE_SIZE:
        _model_cache.popitem(last=False)
    return value


def _finite(rows: list[dict]) -> list[dict]:
    """NaN (pandapower's value on de-energized buses/elements) -> None, which JSON can carry."""
    return [
        {k: (None if isinstance(v, float) and not math.isfinite(v) else v) for k, v in row.items()} for row in rows
    ]


def bus_rows(result: PowerFlowResult) -> list[BusRow]:
    table = result.bus_table()
    return [BusRow(**row) for row in _finite(table.to_dict(orient="records"))]


def powerflow_tables(result: PowerFlowResult) -> dict:
    """Shared by the powerflow and timeseries responses below -- both wrap
    the exact same six other-pandapower-tables extraction, one per
    snapshot for the latter.
    """
    return dict(
        lines=_finite(result.line_table().to_dict(orient="records")),
        transformers=_finite(result.trafo_table().to_dict(orient="records")),
        loads=_finite(result.load_table().to_dict(orient="records")),
        generators=_finite(result.gen_table().to_dict(orient="records")),
        static_generators=_finite(result.sgen_table().to_dict(orient="records")),
        external_grid=_finite(result.ext_grid_table().to_dict(orient="records")),
    )


def topology_response(network: Network) -> TopologyResponse:
    """Bus/line/transformer layout only -- no power flow needed, cheap."""
    topo = compute_topology_layout(network)
    return TopologyResponse(
        nodes=[BusNodeRow(**asdict(n)) for n in topo.nodes],
        edges=[TopologyEdgeRow(**asdict(e)) for e in topo.edges],
    )


def validate_network_response(network: Network) -> ValidateResponse:
    """Every structural issue found (network.validation.validate_network()),
    not just the first -- a pre-flight check the web UI can call before
    even attempting power flow/modal/EMT, so a hand-built network's
    problems can all be shown at once instead of discovered one crash at a
    time. Doesn't run power flow itself, so this is cheap.
    """
    issues = validate_network(network)
    st = service_state(network)
    return ValidateResponse(
        ok=not any(i.severity == "error" for i in issues),
        issues=[NetworkIssueRow(severity=i.severity, message=i.message, affects=list(i.affects)) for i in issues],
        service=ServiceInfo(
            references=list(st.references), energized_buses=sorted(st.energized_buses),
            islands=[IslandRow(buses=sorted(i.buses), reference=i.reference) for i in st.islands],
            lines=list(st.lines), transformers=list(st.transformers), loads=list(st.loads),
            der_units=st.der_units,
        ),
    )


def energized_or_422(network: Network) -> Network:
    """The network the dynamic models are built from: without whatever open
    breakers switch out (network.breakers.energized_network)."""
    try:
        return energized_network(network)
    except NoReferenceUnit as e:
        raise HTTPException(status_code=422, detail=str(e)) from e


def runpp_kwargs(options: PowerFlowOptions | None) -> dict:
    """Validates ``options`` and maps them to ``pandapower.runpp`` keyword
    arguments (422 on anything out of range, before any solve runs).
    """
    if options is None:
        return {}
    if options.algorithm not in POWERFLOW_ALGORITHMS:
        raise HTTPException(status_code=422, detail=f"algorithm must be one of {sorted(POWERFLOW_ALGORITHMS)}")
    if options.init not in POWERFLOW_INITS:
        raise HTTPException(status_code=422, detail=f"init must be one of {list(POWERFLOW_INITS)}")
    if not (0 < options.tolerance_mva <= 1.0):
        raise HTTPException(status_code=422, detail="tolerance_mva must be in (0, 1]")
    max_iteration = options.max_iteration
    if max_iteration != "auto":
        try:
            max_iteration = int(max_iteration)
        except (TypeError, ValueError):
            raise HTTPException(status_code=422, detail="max_iteration must be 'auto' or an integer") from None
        if not (1 <= max_iteration <= 10000):
            raise HTTPException(status_code=422, detail="max_iteration must be in [1, 10000]")
    return dict(
        algorithm=options.algorithm, max_iteration=max_iteration,
        tolerance_mva=options.tolerance_mva, init=options.init,
    )


def solve_powerflow_or_422(network: Network, kwargs: dict) -> PowerFlowResult:
    """``run_power_flow`` already maps non-convergence to ``converged=False``;
    anything else pandapower raises here is a solver/network mismatch (e.g.
    backward/forward sweep on a meshed network), reported as a 422 with the
    solver's own reason instead of a bare 500.
    """
    try:
        return run_power_flow(network, **kwargs)
    except Exception as e:  # noqa: BLE001 -- pandapower raises many unrelated types here
        algo = POWERFLOW_ALGORITHMS.get(kwargs.get("algorithm", "nr"), "power flow")
        raise HTTPException(status_code=422, detail=f"{algo} solver failed: {type(e).__name__}: {e}") from e


def _solver_diagnostics(result: PowerFlowResult) -> dict:
    ppc = getattr(result.net, "_ppc", None)
    if not isinstance(ppc, dict):
        return dict(iterations=None, solve_time_s=None)
    iterations, et = ppc.get("iterations"), ppc.get("et")
    return dict(
        iterations=int(iterations) if iterations is not None else None,
        solve_time_s=float(et) if et is not None else None,
    )


def powerflow_response(network: Network, options: PowerFlowOptions | None = None) -> PowerFlowResponse:
    kwargs = runpp_kwargs(options)
    result = solve_powerflow_or_422(network, kwargs)
    algorithm = kwargs.get("algorithm", "nr")
    if not result.converged:
        return PowerFlowResponse(
            converged=False, buses=[], total_losses_mw=0.0,
            lines=[], transformers=[], loads=[], generators=[], static_generators=[], external_grid=[],
            algorithm=algorithm,
        )
    return PowerFlowResponse(
        converged=True, buses=bus_rows(result), total_losses_mw=result.total_losses_mw(),
        **powerflow_tables(result), algorithm=algorithm, **_solver_diagnostics(result),
    )


def batch_powerflow_response(network: Network, req: BatchPowerFlowRequest) -> BatchPowerFlowResponse:
    """A ramp of ``req.steps + 1`` power flows from the base operating point
    (every scale factor 1.0) to the requested targets, in equal increments --
    each snapshot an independent steady-state solve (no dynamics between
    them), like ``timeseries_response`` but with user-chosen load and DER
    setpoint scaling instead of a fixed load-only sweep.
    """
    kwargs = runpp_kwargs(req.options)
    if not (1 <= req.steps <= BATCH_MAX_STEPS):
        raise HTTPException(status_code=422, detail=f"steps must be in [1, {BATCH_MAX_STEPS}]")
    factors = [req.load_p_scale, req.load_q_scale, *req.der_scale.values()]
    if any(not (0.0 <= f <= BATCH_MAX_SCALE) for f in factors):
        raise HTTPException(status_code=422, detail=f"scale factors must be in [0, {BATCH_MAX_SCALE}]")
    unknown = set(req.der_scale) - {"sm", "gfm", "gfl", "infinite_bus"}
    if unknown:
        raise HTTPException(status_code=422, detail=f"unknown DER unit type(s) in der_scale: {sorted(unknown)}")

    snapshots: list[BatchSnapshot] = []
    for k in range(req.steps + 1):
        frac = k / req.steps
        lp = 1.0 + (req.load_p_scale - 1.0) * frac
        lq = 1.0 + (req.load_q_scale - 1.0) * frac
        ds = {typ: 1.0 + (f - 1.0) * frac for typ, f in req.der_scale.items()}
        scaled = [d for d in network.der_units if d.unit_type.value in ds]
        snap = Snapshot(
            label=f"step {k}",
            load_p_mw={i: ld.p_mw * lp for i, ld in enumerate(network.loads)},
            load_q_mvar={i: ld.q_mvar * lq for i, ld in enumerate(network.loads)},
            der_p_mw={d.id: d.p_set_mw * ds[d.unit_type.value] for d in scaled},
            der_q_mvar={d.id: d.q_set_mvar * ds[d.unit_type.value] for d in scaled},
        )
        result = solve_powerflow_or_422(apply_snapshot(network, snap), kwargs)
        common = dict(label=snap.label, load_p_scale=lp, load_q_scale=lq, der_scale=ds)
        if not result.converged:
            snapshots.append(BatchSnapshot(
                converged=False, buses=[], total_losses_mw=0.0,
                lines=[], transformers=[], loads=[], generators=[], static_generators=[], external_grid=[], **common,
            ))
        else:
            snapshots.append(BatchSnapshot(
                converged=True, buses=bus_rows(result), total_losses_mw=result.total_losses_mw(),
                **powerflow_tables(result), **common,
            ))
    return BatchPowerFlowResponse(snapshots=snapshots)


def build_modal_from_network(network: Network) -> tuple[AssembledSystem, ModalAnalysisResult]:
    """Shared by every modal-analysis response below -- all of them are
    different views over the same (system, modal) pair, not separate models.
    Memoized per network content (see ``_cached``).
    """
    return _cached("modal", network, lambda: _build_modal_uncached(network))


def _build_modal_uncached(network: Network) -> tuple[AssembledSystem, ModalAnalysisResult]:
    network = energized_or_422(network)
    result = run_power_flow(network)
    if not result.converged:
        raise HTTPException(status_code=422, detail="power flow did not converge; can't run modal analysis")
    try:
        system = linearize_network(network, result)
    except NotImplementedError as e:
        raise HTTPException(status_code=501, detail=str(e)) from e
    except ValueError as e:
        # A structural problem with the network itself (e.g. a DER with no
        # transformer connecting it to the rest of the network -- see
        # operating_point.compute_operating_point's own checks), not a
        # solver failure -- a 422 with the actual reason, not a bare 500.
        raise HTTPException(status_code=422, detail=str(e)) from e
    modal = analyze(system.A, system.state_names)
    return system, modal


def named_index_or_422(names: list[str], name: str, what: str) -> int:
    try:
        return names.index(name)
    except ValueError:
        raise HTTPException(status_code=422, detail=f"unknown {what} {name!r}; known: {names}") from None


def mode_or_422(modal: ModalAnalysisResult, mode: int) -> None:
    if not (0 <= mode < len(modal.eigenvalues)):
        raise HTTPException(status_code=422, detail=f"mode must be in [0, {len(modal.eigenvalues) - 1}]")


def modal_response(network: Network) -> ModalResponse:
    system, modal = build_modal_from_network(network)
    table = modal.summary_table().sort_values("real", ascending=False)
    # What kind of mode each one is -- the eigenvalue map groups by it.
    kinds = {c.mode: c for c in classify_modes(modal)}
    # The reference-angle modes are the model's own coordinates, marginal by
    # construction (see modal.reference_angle_modes), so the verdict and the
    # worst real part are read off the physical ones.
    reference = reference_angle_modes(modal)
    # A fleet with nothing regulating frequency has one more marginal
    # direction of the same kind: the common frequency is a free integrator,
    # so it drifts rather than returning. Left in, it reads as a very slow
    # instability, which is the absence of a governor being reported as a
    # dynamic fault.
    if not network.frequency_is_regulated():
        free = unregulated_frequency_mode(modal, reference)
        if free is not None:
            reference = [*reference, free]
    physical = np.array([z for j, z in enumerate(modal.eigenvalues) if j not in set(reference)])
    if not len(physical):
        physical = modal.eigenvalues
    return ModalResponse(
        n_states=system.A.shape[0],
        stable=bool((physical.real < 1e-6).all()),
        max_real_part=float(physical.real.max()),
        reference_modes=reference,
        modes=[
            ModeRow(
                **row,
                category=kinds[row["mode"]].category,
                category_share=kinds[row["mode"]].share,
                category_shares=kinds[row["mode"]].shares,
            )
            for row in table.to_dict(orient="records")
        ],
        categories=[
            ModeCategoryInfo(id=k, label=v["label"], note=v["note"])
            for k, v in reduction.CATEGORIES.items()
        ],
        state_names=system.state_names,
        input_names=system.input_names,
        output_names=system.output_names,
        participation=modal.participation.tolist(),
    )


def modal_sensitivity_response(network: Network, req: SensitivityRequest) -> SensitivityResponse:
    system, modal = build_modal_from_network(network)
    mode_or_422(modal, req.mode)
    sens = eigenvalue_sensitivity(modal, system.A, req.mode)
    return SensitivityResponse(
        mode=sens.mode,
        state_names=modal.state_names,
        matrix=sens.matrix.tolist(),
        top=[SensitivityEntryRow(row_state=e.row_state, col_state=e.col_state, value=e.value) for e in sens.top],
    )


def modal_parameter_sensitivity_response(
    network: Network, req: ParameterSensitivityRequest
) -> ParameterSensitivityResponse:
    """Which physical parameter moves the chosen mode, and which parameters
    the entries of A it is most sensitive to are built from.

    Costs two *component* linearisations per parameter scanned, well under a
    millisecond each: a full 118-bus scan of 1674 parameters is a few seconds,
    less than this request's own eigendecomposition. The request can still be
    narrowed to particular units or parameters when only those are of
    interest.
    """
    net = energized_or_422(network)
    result = run_power_flow(net)
    if not result.converged:
        raise HTTPException(status_code=422, detail="power flow did not converge")
    _, modal = build_modal_from_network(network)
    mode_or_422(modal, req.mode)
    try:
        report = parameter_sensitivity(
            net, result, modal, req.mode,
            units=req.units or None, parameters=req.parameters or None,
            n_entries=req.n_entries,
        )
    except ValueError as e:
        raise HTTPException(status_code=422, detail=str(e)) from e
    data = report.as_dict()
    return ParameterSensitivityResponse(
        **{k: v for k, v in data.items() if k not in ("effects", "entries")},
        effects=[ParameterEffectRow(**e) for e in data["effects"]],
        entries=[EntryParametersRow(**e) for e in data["entries"]],
    )


def modal_mode_shape_response(network: Network, req: ModeShapeRequest) -> ModeShapeResponse:
    _, modal = build_modal_from_network(network)
    mode_or_422(modal, req.mode)
    ms = mode_shape(modal, req.mode)
    return ModeShapeResponse(mode=ms.mode, states=ms.states, angles_deg=ms.angles_deg)


def modal_free_response_response(network: Network, req: FreeResponseRequest) -> FreeResponseResponse:
    if not (0 < req.t_final <= MODAL_MAX_T_FINAL):
        raise HTTPException(status_code=422, detail=f"t_final must be in (0, {MODAL_MAX_T_FINAL}]")
    _, modal = build_modal_from_network(network)
    idx = named_index_or_422(modal.state_names, req.perturb_state, "state")
    t = np.linspace(0.0, req.t_final, 200)
    x_t = free_response(modal, idx, req.offset, t)
    if req.plot_states:
        matches = [named_index_or_422(modal.state_names, n, "state") for n in req.plot_states]
    else:
        matches = [i for i, n in enumerate(modal.state_names) if req.state_filter in n]
    if not matches:
        raise HTTPException(status_code=422, detail=f"state_filter {req.state_filter!r} matches no state names")
    series = {modal.state_names[i]: x_t[i, :].tolist() for i in matches}
    return FreeResponseResponse(perturb_state=modal.state_names[idx], t=t.tolist(), series=series)


def modal_step_response_response(network: Network, req: StepResponseRequest) -> StepResponseResponse:
    if not (0 < req.t_final <= MODAL_MAX_T_FINAL):
        raise HTTPException(status_code=422, detail=f"t_final must be in (0, {MODAL_MAX_T_FINAL}]")
    output_names = list(req.output_names) or ([req.output_name] if req.output_name else [])
    if not output_names:
        raise HTTPException(status_code=422, detail="give output_name or a non-empty output_names")
    system, _ = build_modal_from_network(network)
    named_index_or_422(system.input_names, req.input_name, "input")
    for name in output_names:
        named_index_or_422(system.output_names, name, "output")
    t = np.linspace(0.0, req.t_final, 200)
    series = {name: step_response(system, req.input_name, name, req.amplitude, t).tolist() for name in output_names}
    return StepResponseResponse(
        t=t.tolist(), y=series[output_names[0]], series=series,
        note=_inert_input_note(network, system, req.input_name),
    )


def _inert_input_note(network: Network, system, input_name: str) -> str:
    """Why this input moves nothing, when it moves nothing.

    A flat response is a real answer, but on its own it looks like a broken
    plot. The usual cause is an input whose regulator is not fitted: a
    machine's P_ref is its *governor's* setpoint, so with no governor there
    is nothing for it to act on and its column of B is exactly zero.
    """
    k = system.input_names.index(input_name)
    if np.any(system.B[:, k]) or np.any(system.D[:, k]):
        return ""
    unit = input_name.split("_{")[-1].rstrip("}")
    der = next((d for d in network.der_units if f"{d.unit_type.value.upper()}_{d.id}" == unit), None)
    if der is not None and der.unit_type.value == "sm" and der.governor_model == "none":
        if input_name.startswith("P_ref"):
            return (f"{unit} has no governor fitted, so its P_ref has nothing to act on and the "
                    f"response is exactly zero. Fit a governor to make the power reference live.")
    return (f"{input_name} does not reach the model: its column of B is zero, so every output "
            f"stays flat. This is the model saying the path is not there, not a failed solve.")


def timeseries_response(network: Network) -> TimeSeriesResponse:
    snapshots = scale_loads(network, {"70% load": 0.7, "100% load": 1.0, "130% load": 1.3})
    result = run_time_series(network, snapshots)

    def _snapshot(snap) -> TimeSeriesSnapshot:
        pf_result = result.results[snap.label]
        if not result.converged[snap.label]:
            return TimeSeriesSnapshot(
                label=snap.label, converged=False, buses=[], total_losses_mw=0.0,
                lines=[], transformers=[], loads=[], generators=[], static_generators=[], external_grid=[],
            )
        return TimeSeriesSnapshot(
            label=snap.label, converged=True, buses=bus_rows(pf_result),
            total_losses_mw=pf_result.total_losses_mw(), **powerflow_tables(pf_result),
        )

    return TimeSeriesResponse(snapshots=[_snapshot(snap) for snap in snapshots])


def build_nonlinear_model_from_network(network: Network) -> NonlinearNetworkModel:
    """Memoized per network content (see ``_cached``)."""
    return _cached("nonlinear", network, lambda: _build_nonlinear_uncached(network))


def _build_nonlinear_uncached(network: Network) -> NonlinearNetworkModel:
    network = energized_or_422(network)
    result = run_power_flow(network)
    if not result.converged:
        raise HTTPException(status_code=422, detail="power flow did not converge; can't build a time-domain model")
    try:
        return build_nonlinear_network(network, result)
    except NotImplementedError as e:
        raise HTTPException(status_code=501, detail=str(e)) from e
    except ValueError as e:
        # Same structural-network-problem case build_modal_from_network()
        # handles -- a 422 with the actual reason, not a bare 500.
        raise HTTPException(status_code=422, detail=str(e)) from e


def find_state_or_422(model: NonlinearNetworkModel, name_contains: str) -> int:
    try:
        return find_state_index(model, name_contains)
    except ValueError as e:
        raise HTTPException(status_code=422, detail=str(e)) from e


def states_response(network: Network) -> StatesResponse:
    """State/input/output names for the nonlinear (EMT) model -- lets
    the frontend populate its pickers without hard-coding per-network names
    (which differ: e.g. only synchronous machines have ``dw_r_*`` states,
    and which DER id is the slack varies).
    """
    model = build_nonlinear_model_from_network(network)
    return StatesResponse(
        state_names=model.state_names, input_names=model.input_names, output_names=model.output_names,
        measurements=[MeasurementInfo(**vars(m)) for m in measurement_set(network).catalog()],
    )


def measurement_set(network: Network) -> MeasurementSet:
    """The measurement outputs of ``network``'s EMT model (memoized with it)."""
    return _cached("measurements", network, lambda: MeasurementSet(build_nonlinear_model_from_network(network)))


def _check_measurements(network: Network, names: list[str]) -> MeasurementSet | None:
    if not names:
        return None
    ms = measurement_set(network)
    unknown = set(names) - set(ms.names())
    if unknown:
        raise HTTPException(status_code=422, detail=f"unknown measurement name(s): {sorted(unknown)}")
    return ms


class _EmtRun:
    """A resolved EMT disturbance: what to integrate from t=0, and its
    linear equivalent. ``model`` is the undisturbed (pre-event) model;
    ``sim`` is the one integrated -- the same object, except for a network
    event that changes the model (see timedomain.events)."""

    def __init__(
        self, model: NonlinearNetworkModel, sim: NonlinearNetworkModel, x0: np.ndarray,
        u_exo_fn: Callable[[float], np.ndarray] | None, name: str,
        dx0: np.ndarray | None, du: np.ndarray | None, linear_note: str | None = None,
    ) -> None:
        self.model, self.sim, self.x0, self.u_exo_fn, self.name = model, sim, x0, u_exo_fn, name
        # Linear equivalent: an initial-state offset and/or an input step
        # (both None when the event has none -- linear_note says why).
        self.dx0, self.du, self.linear_note = dx0, du, linear_note

    @property
    def changes_model(self) -> bool:
        return self.sim is not self.model


def _resolve_emt_run(model: NonlinearNetworkModel, req: EmtRequest) -> _EmtRun:
    """Shared by :func:`emt_response` and :func:`prepare_emt_live` so the two
    never resolve a disturbance differently."""
    x0 = model.initial_state()
    n_x, n_u = len(x0), len(model.input_names)

    if req.perturb_kind == "state":
        idx = find_state_or_422(model, req.perturb_name)
        x0 = x0.copy()
        x0[idx] += req.perturb_offset
        dx0 = np.zeros(n_x)
        dx0[idx] = req.perturb_offset
        return _EmtRun(model, model, x0, None, model.state_names[idx], dx0, np.zeros(n_u))
    if req.perturb_kind == "input":
        # A permanent step in one exogenous reference, held from t=0 for the
        # whole run (not just an initial-condition offset) -- the system
        # starts at its own equilibrium and that equilibrium itself moves,
        # the more standard "P_ref step" kind of disturbance test.
        idx = named_index_or_422(model.input_names, req.perturb_name, "input")
        perturbed_u_exo = model.default_u_exo()
        perturbed_u_exo[idx] += req.perturb_offset
        du = np.zeros(n_u)
        du[idx] = req.perturb_offset
        return _EmtRun(
            model, model, x0, lambda t: perturbed_u_exo, model.input_names[idx], np.zeros(n_x), du,
        )
    if req.perturb_kind == "event":
        if req.event is None:
            raise HTTPException(status_code=422, detail="perturb_kind 'event' needs an event")
        try:
            applied = apply_event(model, NetworkEvent(**req.event.model_dump()))
        except EventError as e:
            raise HTTPException(status_code=422, detail=str(e)) from e
        except RuntimeError as e:
            raise HTTPException(status_code=422, detail=f"couldn't start from the post-event network ({e})") from e
        if applied.dx0 is not None:
            return _EmtRun(model, applied.model, applied.x0, None, applied.description, applied.dx0, np.zeros(n_u))
        return _EmtRun(
            model, applied.model, applied.x0, None, applied.description, None, None,
            linear_note="no linear overlay for this event: it changes the network itself (topology or load), "
            "which a model linearized around the pre-event point can't represent",
        )
    raise HTTPException(status_code=422, detail="perturb_kind must be 'state', 'input' or 'event'")


def _plot_state_names(model: NonlinearNetworkModel, req: EmtRequest) -> list[str]:
    if req.plot_states:
        unknown = set(req.plot_states) - set(model.state_names)
        if unknown:
            raise HTTPException(status_code=422, detail=f"unknown state name(s): {sorted(unknown)}")
        return list(req.plot_states)
    return [n for n in model.state_names if "dw_r" in n] or model.state_names


def _check_names(names: list[str], known: list[str], what: str) -> list[str]:
    unknown = set(names) - set(known)
    if unknown:
        raise HTTPException(status_code=422, detail=f"unknown {what} name(s): {sorted(unknown)}")
    return list(names)


def _removed_fill(name: str) -> float | None:
    """What a measurement of an element an event removed reads: no power
    flows through it; anything else is undefined (a gap in the plot)."""
    return 0.0 if name.startswith(("P_", "Q_")) else None


def _post_measurements(run: _EmtRun, ms: MeasurementSet | None, names: list[str]) -> tuple[MeasurementSet | None, list[str]]:
    """The measurement set to evaluate after t=0, and which of ``names`` it has."""
    if ms is None or not names:
        return None, []
    post = MeasurementSet(run.sim) if run.changes_model else ms
    known = set(post.names())
    return post, [n for n in names if n in known]


def _validate_emt_t_final_and_n_points(req: EmtRequest) -> int:
    """Shared ``t_final``/``dt`` bounds check -- returns the resolved
    sample count (one-shot endpoint) / step-count safety bound (live
    endpoint), both keyed off the same ``EMT_MIN/MAX_N_POINTS`` constants.
    """
    if not (req.t_final > 0 and math.isfinite(req.t_final)):
        raise HTTPException(status_code=422, detail="t_final must be a positive, finite number")
    if not (0 <= req.t_pre <= EMT_MAX_T_PRE):
        raise HTTPException(status_code=422, detail=f"t_pre must be in [0, {EMT_MAX_T_PRE}]")

    if req.dt is None:
        return EMT_DEFAULT_N_POINTS
    if req.dt <= 0:
        raise HTTPException(status_code=422, detail="dt must be > 0")
    n_points = max(2, round(req.t_final / req.dt) + 1)
    if not (EMT_MIN_N_POINTS <= n_points <= EMT_MAX_N_POINTS):
        raise HTTPException(
            status_code=422,
            detail=f"dt={req.dt} gives {n_points} samples over t_final={req.t_final}; "
            f"need between {EMT_MIN_N_POINTS} and {EMT_MAX_N_POINTS}",
        )
    return n_points


def _equilibrium_signals(model: NonlinearNetworkModel) -> tuple[np.ndarray, dict[str, float], dict[str, float]]:
    """``(x_eq, inputs_eq, outputs_eq)`` at the undisturbed operating point --
    what every signal reads before the disturbance at t=0.
    """
    x_eq = model.initial_state()
    z_eq, u_eq = model.solve_algebraic(x_eq, model.default_u_exo())
    inputs_eq, outputs_eq = model._inputs_and_outputs(x_eq, z_eq, u_eq)
    return x_eq, inputs_eq, outputs_eq


PreSamples = tuple[list[float], dict[str, list[float]], dict[str, list[float]], dict[str, list[float]]]


def _pre_disturbance_samples(
    model: NonlinearNetworkModel, t_pre: float, plot_states: list[str], plot_inputs: list[str], plot_outputs: list[str],
    n: int = 2,
) -> PreSamples:
    """Equilibrium samples from ``t=-t_pre`` to ``t=0`` (the instant just
    before the disturbance), for prepending to a trajectory: a plot then shows
    the flat pre-disturbance x0 and a vertical jump/step at t=0. Two samples
    are enough for constant signals; ``n`` more are used when measurements
    need them (3-phase voltages keep oscillating before T0).
    """
    x_eq, inputs_eq, outputs_eq = _equilibrium_signals(model)
    t = np.linspace(-t_pre, 0.0, max(2, n)).tolist()
    states = {name: [float(x_eq[model.state_names.index(name)])] * len(t) for name in plot_states}
    inputs = {name: [inputs_eq[name]] * len(t) for name in plot_inputs}
    outputs = {name: [outputs_eq[name]] * len(t) for name in plot_outputs}
    return t, states, inputs, outputs


def _pre_points(t_pre: float, dt: float, measurements: list[str]) -> int:
    """How many pre-disturbance samples: 2, or the trajectory's own spacing
    when measurements are plotted."""
    return max(2, int(round(t_pre / dt)) + 1) if measurements and dt > 0 else 2


def _prepend(pre: PreSamples, t: list[float], series: dict, inputs: dict, outputs: dict) -> PreSamples:
    pre_t, pre_s, pre_i, pre_o = pre
    return (
        pre_t + list(t),
        {n: pre_s.get(n, []) + list(v) for n, v in series.items()},
        {n: pre_i.get(n, []) + list(v) for n, v in inputs.items()},
        {n: pre_o.get(n, []) + list(v) for n, v in outputs.items()},
    )


def linear_overlay(
    network: Network, run: _EmtRun, t: np.ndarray,
    plot_states: list[str], plot_inputs: list[str], plot_outputs: list[str],
) -> LinearOverlay:
    """The linearized model's response to the same disturbance the EMT run
    applied, around the same equilibrium, in absolute (not deviation) units so
    it overlays the nonlinear trajectories directly: ``x = x_eq + dx``,
    ``y = y_eq + C dx + D du``. A state perturbation is an initial condition
    ``dx(0)``; an input perturbation is a step ``du`` held from t=0 -- both
    one ``lsim`` of the full ``(A, B, C, D)``.

    Relies on the linear and nonlinear models sharing their state/input/
    output naming and ordering (both are assembled from the same blocks by
    the same interconnection code) -- checked, not assumed; outputs missing
    from the linear model are left out rather than guessed.
    """
    if run.dx0 is None:
        raise HTTPException(status_code=422, detail=run.linear_note or "no linear equivalent")
    model = run.model
    system, _ = build_modal_from_network(network)
    x_eq, inputs_eq, outputs_eq = _equilibrium_signals(model)
    n_x, n_u = system.A.shape[0], system.B.shape[1]
    if len(x_eq) != n_x or list(system.state_names) != list(model.state_names):
        raise HTTPException(
            status_code=501, detail="linearized and nonlinear models don't share a state vector for this network"
        )

    t = np.asarray(t, dtype=float)
    if list(system.input_names) != list(model.input_names):
        raise HTTPException(status_code=501, detail="linearized and nonlinear models don't share their inputs")
    dx0, du = run.dx0, run.du

    sys = scipy.signal.StateSpace(system.A, system.B, system.C, system.D)
    U = np.tile(du, (len(t), 1))
    _, dy, dx = scipy.signal.lsim(sys, U=U, T=t - t[0], X0=dx0, interp=False)
    dx = np.asarray(dx).reshape(len(t), n_x)
    dy = np.asarray(dy).reshape(len(t), -1)

    out_names = list(system.output_names)
    in_names = list(system.input_names)
    series = {n: (x_eq[model.state_names.index(n)] + dx[:, model.state_names.index(n)]).tolist() for n in plot_states}
    inputs = {n: [inputs_eq[n] + (du[in_names.index(n)] if n in in_names else 0.0)] * len(t) for n in plot_inputs}
    outputs = {n: (outputs_eq[n] + dy[:, out_names.index(n)]).tolist() for n in plot_outputs if n in out_names}
    return LinearOverlay(t=t.tolist(), series=series, inputs=inputs, outputs=outputs)


def _trajectories(
    run: _EmtRun, t: np.ndarray, x: np.ndarray, names: list[str],
) -> dict[str, list[float | None]]:
    """State trajectories of ``run.sim`` by name; a state of an element the
    event removed reads None (a gap) after t=0."""
    idx = {n: i for i, n in enumerate(run.sim.state_names)}
    return {n: x[idx[n], :].tolist() if n in idx else [None] * len(t) for n in names}


def _solver_kwargs(req: EmtRequest, with_max_step: bool = True) -> dict:
    """The integrator settings one request asks for.

    ``max_step`` is left out entirely when the request doesn't set it: SciPy
    reads it as "no bound", and passing a bound the user didn't ask for makes
    the solver take steps it didn't need. ``with_max_step=False`` leaves it
    out either way, for the live stream, which passes its own (see
    :func:`prepare_emt_live`) and would otherwise get it twice.
    """
    if req.solver not in SOLVERS:
        raise HTTPException(
            status_code=422,
            detail=f"unknown solver {req.solver!r}; available: {sorted(SOLVERS)}",
        )
    kwargs = dict(method=req.solver, rtol=req.rtol, atol=req.atol, first_step=EMT_FIRST_STEP)
    if with_max_step and req.max_step is not None:
        kwargs["max_step"] = req.max_step
    return kwargs


def _fixed_step_or_422(req: EmtRequest, n_points: int) -> float:
    """The fixed step a request asks for, defaulting to one step per plotted
    point, and refused up front if it would take an unreasonable number."""
    step = req.fixed_step if req.fixed_step is not None else req.t_final / max(n_points - 1, 1)
    n_steps = req.t_final / step
    if n_steps > EMT_MAX_FIXED_STEPS:
        raise HTTPException(
            status_code=422,
            detail=(
                f"a {step:g} s step over {req.t_final:g} s is {n_steps:.0f} steps, past this "
                f"endpoint's {EMT_MAX_FIXED_STEPS:,} limit. Use a larger step, a shorter run, or "
                "variable stepping, which spends steps only where the trajectory needs them."
            ),
        )
    return step


def solvers_response() -> SolversResponse:
    """The integrator catalogue -- what the time-domain page's solver picker
    is built from, so it never hard-codes a method name or a default."""
    return SolversResponse(
        solvers=[
            SolverInfo(id=k, label=v["label"], note=v["note"], implicit=v["implicit"])
            for k, v in SOLVERS.items()
        ],
        default=DEFAULT_SOLVER,
        default_rtol=EMT_SIMULATE_KWARGS["rtol"],
        default_atol=EMT_SIMULATE_KWARGS["atol"],
    )


def emt_response(network: Network, req: EmtRequest) -> EmtResponse:
    n_points = _validate_emt_t_final_and_n_points(req)
    model = build_nonlinear_model_from_network(network)
    ms = _check_measurements(network, req.plot_measurements)
    run = _resolve_emt_run(model, req)
    plot_states = _plot_state_names(model, req)
    plot_inputs = _check_names(req.plot_inputs, model.input_names, "input")
    plot_outputs = _check_names(req.plot_outputs, model.output_names, "output")

    t_eval = np.linspace(0.0, req.t_final, n_points)
    try:
        if req.stepping == "fixed":
            step = _fixed_step_or_422(req, n_points)
            sim = simulate_fixed_step(
                run.sim, (0.0, req.t_final), step, x0=run.x0, u_exo_fn=run.u_exo_fn
            )
        else:
            sim = simulate(
                run.sim, (0.0, req.t_final), x0=run.x0, u_exo_fn=run.u_exo_fn, t_eval=t_eval,
                **_solver_kwargs(req),
            )
    except RuntimeError as e:
        # The coupled Newton solve (see timedomain/emt.py) can fail to
        # converge for a large enough perturbation -- a trajectory that
        # fails mid-integration has no partial result worth returning, so
        # this surfaces as a clean 422 instead of an unhandled 500.
        raise HTTPException(
            status_code=422,
            detail=f"nonlinear solver did not converge for this perturbation ({e}); try a smaller offset",
        ) from e
    series = _trajectories(run, sim.t, sim.x, plot_states)

    # Inputs/outputs are opt-in (empty list = skip) -- recovering them costs
    # about as much again as the integration itself (see
    # recover_inputs_and_outputs's docstring). They come from the same
    # u_exo_fn the trajectory was integrated with, and measurements from
    # the same per-sample algebraic solve.
    inputs: dict[str, list] = {}
    outputs: dict[str, list] = {}
    measurements: dict[str, list] = {}
    ms_post, meas_present = _post_measurements(run, ms, req.plot_measurements)
    if plot_inputs or plot_outputs or ms is not None:
        all_inputs, all_outputs, all_meas = run.sim.recover_signals(
            sim.t, sim.x, u_exo_fn=run.u_exo_fn, measure=ms_post.evaluator(meas_present) if ms_post else None,
        )
        n = len(sim.t)
        if ms is not None:
            smoothed = ms_post.smooth(sim.t, all_meas) if meas_present else {}
            measurements = {
                m: smoothed[m].tolist() if m in smoothed else [_removed_fill(m)] * n for m in req.plot_measurements
            }
        inputs = {m: all_inputs[m].tolist() if m in all_inputs else [None] * n for m in plot_inputs}
        outputs = {m: all_outputs[m].tolist() if m in all_outputs else [None] * n for m in plot_outputs}

    actual_dt = float(sim.t[1] - sim.t[0]) if len(sim.t) > 1 else req.t_final
    linear, linear_note = None, None
    if req.linear_overlay:
        if run.dx0 is None:
            linear_note = run.linear_note
        else:
            linear = linear_overlay(network, run, sim.t, list(series), list(inputs), list(outputs))
    t_out = sim.t.tolist()
    if req.t_pre > 0:
        n_pre = _pre_points(req.t_pre, actual_dt, req.plot_measurements)
        pre = _pre_disturbance_samples(model, req.t_pre, list(series), list(inputs), list(outputs), n=n_pre)
        if ms is not None:
            pre_meas = ms.before_disturbance(req.plot_measurements, np.array(pre[0]))
            measurements = {n: pre_meas[n].tolist() + v for n, v in measurements.items()}
        t_out, series, inputs, outputs = _prepend(pre, t_out, series, inputs, outputs)
        if linear is not None:
            lt, ls, li, lo = _prepend(pre, linear.t, linear.series, linear.inputs, linear.outputs)
            linear = LinearOverlay(t=lt, series=ls, inputs=li, outputs=lo)
    return EmtResponse(
        perturbed=run.name, perturb_kind=req.perturb_kind, dt=actual_dt,
        state_names=list(series.keys()), t=t_out, series=series, inputs=inputs, outputs=outputs,
        measurements=measurements, linear=linear, linear_note=linear_note,
    )


class _EmtLivePlan:
    """Everything :func:`emt_live_stream` needs, already validated --
    kept apart from the stream itself so every possible 422 is raised
    synchronously, in the route handler, *before* ``StreamingResponse``
    starts sending anything. Raising from inside an async generator after
    it has already started streaming can't turn into a normal JSON error
    response any more; the client just sees a broken stream.
    """

    def __init__(
        self, run: _EmtRun, plot_states: list[str], plot_inputs: list[str], plot_outputs: list[str],
        max_step: float, t_final: float, network: Network, req: EmtRequest,
        measurements: MeasurementSet | None = None,
    ) -> None:
        self.run, self.model = run, run.model
        self.plot_states, self.plot_inputs, self.plot_outputs = plot_states, plot_inputs, plot_outputs
        self.max_step, self.t_final = max_step, t_final
        self.network, self.req = network, req
        self.measurements = measurements


def prepare_emt_live(network: Network, req: EmtRequest) -> _EmtLivePlan:
    """Resolves and validates everything a live EMT run needs -- same
    ``t_final``/``dt`` bounds, same perturbation resolution, same
    plot_states/plot_inputs/plot_outputs name-checking as
    :func:`emt_response`, reusing the exact same helpers so the two paths
    can never validate differently. ``req.dt``, when given, is read as a
    *step-size ceiling* here (``max_step``) rather than an exact output
    sample spacing -- see ``timedomain.emt.simulate_steps``'s own
    docstring for why that mostly controls redraw granularity, not
    runtime, for a stiff system.
    """
    _validate_emt_t_final_and_n_points(req)  # same bounds check; the resolved n_points isn't used here
    model = build_nonlinear_model_from_network(network)
    ms = _check_measurements(network, req.plot_measurements)
    run = _resolve_emt_run(model, req)
    # A live run caps its step so the trace arrives smoothly rather than in
    # a few big jumps -- that is a display concern, not the user's solver
    # setting, so when they have set one too, the tighter of the two wins.
    max_step = req.dt if req.dt else req.t_final / EMT_DEFAULT_N_POINTS
    if req.max_step is not None:
        max_step = min(max_step, req.max_step)
    return _EmtLivePlan(
        run, _plot_state_names(model, req),
        _check_names(req.plot_inputs, model.input_names, "input"),
        _check_names(req.plot_outputs, model.output_names, "output"),
        max_step, req.t_final, network=network, req=req, measurements=ms,
    )


async def emt_live_stream(plan: _EmtLivePlan, perturb_kind: str, request: Request) -> AsyncIterator[str]:
    """Newline-delimited JSON, one line per accepted solver step:
    ``{"t": ..., "states": {...}, "inputs": {...}, "outputs": {...}}``,
    a final ``{"done": true, ...}`` line on success, or a single
    ``{"error": "..."}`` line (in place of ``done``) if the coupled
    Newton solve fails partway through -- there's no partial result worth
    keeping at that point, matching :func:`emt_response`'s own 422 on the
    same failure, but a stream that's already started can't become an
    HTTP error status any more, so this is the best it can signal that
    over the wire.
    """
    run = plan.run
    # Resolved once, not per step: state_names.index() is an O(n) scan and
    # this loop can run thousands of times. A state the event removed is None.
    sim_idx = {n: i for i, n in enumerate(run.sim.state_names)}
    state_idx = [sim_idx.get(n) for n in plan.plot_states]

    # Pre-disturbance equilibrium samples (t_pre > 0), streamed first -- same
    # line shape as a solver step, so the client plots them like any other.
    t_pre = plan.req.t_pre
    meas_names = list(plan.req.plot_measurements) if plan.measurements is not None else []
    pre = None
    if t_pre > 0:
        pre = _pre_disturbance_samples(
            plan.model, t_pre, plan.plot_states, plan.plot_inputs, plan.plot_outputs,
            n=_pre_points(t_pre, plan.max_step, meas_names),
        )
        pre_t, pre_s, pre_i, pre_o = pre
        pre_m = plan.measurements.before_disturbance(meas_names, np.array(pre_t)) if meas_names else {}
        for k, t_k in enumerate(pre_t):
            yield json.dumps({
                "t": t_k,
                "states": {n: v[k] for n, v in pre_s.items()},
                "inputs": {n: v[k] for n, v in pre_i.items()},
                "outputs": {n: v[k] for n, v in pre_o.items()},
                "measurements": {n: float(v[k]) for n, v in pre_m.items()},
            }) + "\n"
    ms_post, present = _post_measurements(run, plan.measurements, meas_names)
    measure = ms_post.evaluator(present) if present else None
    smooth = ms_post.smoother(present) if present else None
    fill = {n: _removed_fill(n) for n in meas_names if n not in present}

    n_steps = 0
    truncated: str | None = None
    try:
        for step in simulate_steps(
            run.sim, (0.0, plan.t_final), x0=run.x0, u_exo_fn=run.u_exo_fn,
            max_step=plan.max_step, **_solver_kwargs(plan.req, with_max_step=False),
        ):
            if await request.is_disconnected():
                return  # client hit "Stop" / navigated away -- stop computing, not an error
            if n_steps >= EMT_LIVE_MAX_STEPS:
                # Not an error: everything traced so far is a valid trajectory.
                # End the stream normally and say where it stopped.
                truncated = (
                    f"stopped at t = {step.t:.4g} s of {plan.t_final:g} s: the trace reached the "
                    f"{EMT_LIVE_MAX_STEPS}-step limit. Everything up to there is plotted; run a shorter "
                    "duration, or the same run without live tracing, to see the rest."
                )
                break
            n_steps += 1
            meas = dict(fill)
            if measure is not None:
                meas.update(smooth(step.t, measure(step.x, step.z, step.u)))
            yield json.dumps({
                "t": step.t,
                "states": {n: (float(step.x[i]) if i is not None else None) for n, i in zip(plan.plot_states, state_idx)},
                "inputs": {n: step.inputs.get(n) for n in plan.plot_inputs},
                "outputs": {n: step.outputs.get(n) for n in plan.plot_outputs},
                "measurements": meas,
            }) + "\n"
    except RuntimeError as e:
        yield json.dumps({"error": f"nonlinear solver did not converge for this perturbation ({e}); try a smaller offset"}) + "\n"
        return
    done: dict = {"done": True, "perturbed": run.name, "perturb_kind": perturb_kind, "n_steps": n_steps}
    if truncated:
        done["truncated"] = truncated
    if plan.req.linear_overlay and truncated:
        done["linear_error"] = "not computed: the trace was cut short, so there is no full span to overlay it on"
    elif plan.req.linear_overlay:
        # Computed once the nonlinear run has finished, on an even grid over
        # the same span (the solver's own adaptive step times are too
        # irregular to be worth reproducing for a linear overlay).
        if run.dx0 is None:
            done["linear_error"] = run.linear_note
        else:
            try:
                grid = np.linspace(0.0, plan.t_final, EMT_DEFAULT_N_POINTS)
                lin = linear_overlay(plan.network, run, grid, plan.plot_states, plan.plot_inputs, plan.plot_outputs)
                if pre is not None:
                    lt, ls, li, lo = _prepend(pre, lin.t, lin.series, lin.inputs, lin.outputs)
                    lin = LinearOverlay(t=lt, series=ls, inputs=li, outputs=lo)
                done["linear"] = lin.model_dump()
            except HTTPException as e:
                done["linear_error"] = str(e.detail)
    yield json.dumps(done) + "\n"


# --- Model-order reduction ----------------------------------------------------


def model_levels_response() -> ModelLevelsResponse:
    """The reduction catalogue. Static -- it describes the models this
    build knows how to make, not anything about a particular network -- so
    the web UI fetches it once at boot and never hard-codes a level name or
    a state group of its own.
    """
    return ModelLevelsResponse(
        elements=[ElementModelInfo(**reduction.describe(kind)) for kind in reduction.ELEMENTS]
    )


def model_summary_response(network: Network) -> ModelSummaryResponse:
    """What the network's current settings add up to.

    The model class is derived, never chosen: a run is EMT when the passive
    network still integrates and every unit keeps its own fast electrical
    states, RMS when neither does, and Mixed in between. Naming it here
    rather than in the UI keeps the one definition in one place -- it is
    also what decides whether the time-domain page offers instantaneous
    3-phase waveforms, which mean nothing in a quasi-stationary model.
    """
    models = network.models
    net_group_modes = models.group_modes_for("network")
    network_dynamic = models.network_is_dynamic

    units: list[UnitModelRow] = []
    unit_fast_dynamic = False
    for der in network.der_units:
        kind = der.unit_type.value
        if kind not in reduction.ELEMENTS:
            continue
        group_modes = network.unit_group_modes(der)
        fast = "trafo_current" if kind == "sm" else "trafo_current"
        if group_modes.get(fast) == reduction.DYNAMIC:
            unit_fast_dynamic = True
        # This unit's own catalogue: a machine's AVR and PSS groups depend on
        # which regulator models it carries, so their state counts do too.
        element = network.unit_element(der)
        n_states = sum(
            len(element.group(gid).symbols)
            for gid, mode in group_modes.items() if mode == reduction.DYNAMIC
        )
        units.append(UnitModelRow(
            id=der.id,
            unit_type=kind,
            label=f"{TYPE_LABEL.get(kind, kind)} (unit {der.id})",
            level=network.unit_level(der),
            modes=group_modes,
            n_states=n_states,
            exciter=der.exciter_model if kind == "sm" else None,
            pss=der.pss_model if kind == "sm" else None,
            governor=der.governor_model if kind == "sm" else None,
            controller=der.controller_model if kind == "gfm" else None,
        ))

    if network_dynamic and unit_fast_dynamic:
        model_class = "EMT"
    elif not network_dynamic and not unit_fast_dynamic:
        model_class = "RMS"
    else:
        model_class = "Mixed"

    # Both counts come from the blocks that were actually built, never from
    # a sum done by hand here -- but from the *nonlinear* model, which is
    # cheap to build (a fraction of a second even at 118 buses) and already
    # cached for the states endpoint. Linearising a second time just to
    # count the full-order states would put a 20-second build behind a badge
    # that is drawn on every page view; the full-order model is built once,
    # deliberately, by the adequacy check.
    model = build_nonlinear_model_from_network(network)
    n_states = sum(b.comp.n_states for b in model.blocks)
    # A machine's full-order count depends on its regulators, so it is read
    # off that unit's own catalogue where the block can be traced back to one.
    unit_elements = {
        f"{der.unit_type.value.upper()}_{der.id}": network.unit_element(der)
        for der in network.der_units if der.unit_type.value in reduction.ELEMENTS
    }
    n_full = 0
    for block in model.blocks:
        kind = block.kind.removesuffix("_slack")
        if kind in reduction.ELEMENTS:
            element = unit_elements.get(block.name) or reduction.element(kind)
            n_full += sum(len(g.symbols) for g in element.groups)
        elif kind in ("node", "line", "load", "shunt"):
            n_full += 2  # every passive element is a dq pair at full order
        else:
            n_full += block.comp.n_states  # a frame or an infinite bus: nothing reducible

    return ModelSummaryResponse(
        model_class=model_class,
        network_level=reduction.element("network").matching_level(net_group_modes),
        network_modes=net_group_modes,
        network_frequency=models.network_frequency,
        units=units,
        n_states=n_states,
        n_states_full=n_full,
    )


def adequacy_response(network: Network, band_hz: float) -> AdequacyResponse:
    """Judge the chosen reduction against this network's own full-order
    model (g2elin_core.modal.adequacy). Builds two linearizations and two
    eigendecompositions, so it is the most expensive endpoint here -- it is
    meant to be run once when a level is chosen, not on every page view.
    """
    net = energized_or_422(network)
    result = run_power_flow(net)
    if not result.converged:
        raise HTTPException(
            status_code=422, detail="power flow did not converge; can't check the model against it"
        )
    try:
        report = check_adequacy(net, result, band_hz=band_hz)
    except NotImplementedError as e:
        raise HTTPException(status_code=501, detail=str(e)) from e
    except ValueError as e:
        raise HTTPException(status_code=422, detail=str(e)) from e
    data = report.as_dict()
    return AdequacyResponse(
        **{k: v for k, v in data.items() if k not in ("risks", "modes")},
        risks=[StateRiskRow(**r) for r in data["risks"]],
        modes=[ModePairRow(**m) for m in data["modes"]],
    )
