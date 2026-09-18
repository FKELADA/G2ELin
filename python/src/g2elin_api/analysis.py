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
from collections import OrderedDict
from dataclasses import asdict
from typing import AsyncIterator, Callable

import numpy as np
import scipy.signal
from fastapi import HTTPException, Request

from g2elin_core.interconnect import AssembledSystem
from g2elin_core.modal import ModalAnalysisResult, analyze, eigenvalue_sensitivity, free_response, mode_shape, step_response
from g2elin_core.network.schema import Network
from g2elin_core.network.topology import compute_topology_layout
from g2elin_core.network.validation import validate_network
from g2elin_core.pipeline import linearize_network
from g2elin_core.powerflow import run_power_flow
from g2elin_core.powerflow.pandapower_adapter import PowerFlowResult
from g2elin_core.timedomain import NonlinearNetworkModel, build_nonlinear_network, find_state_index, simulate, simulate_steps
from g2elin_core.timeseries import Snapshot, apply_snapshot, run_time_series, scale_loads

from .schemas import (
    BatchPowerFlowRequest,
    BatchPowerFlowResponse,
    BatchSnapshot,
    BusNodeRow,
    BusRow,
    EmtRequest,
    EmtResponse,
    FreeResponseRequest,
    FreeResponseResponse,
    LinearOverlay,
    ModalResponse,
    ModeRow,
    ModeShapeRequest,
    ModeShapeResponse,
    NetworkIssueRow,
    PowerFlowOptions,
    PowerFlowResponse,
    SensitivityEntryRow,
    SensitivityRequest,
    SensitivityResponse,
    StatesResponse,
    StepResponseRequest,
    StepResponseResponse,
    TimeSeriesResponse,
    TimeSeriesSnapshot,
    TopologyEdgeRow,
    TopologyResponse,
    ValidateResponse,
)

# Bounds on user-controlled simulation cost, enforced server-side rather than
# just documented -- these endpoints run a real nonlinear DAE integration
# synchronously inside one HTTP request, so an unbounded t_final directly
# translates to an unbounded request time.
EMT_MAX_T_FINAL = 3.0
EMT_MIN_N_POINTS = 10
EMT_MAX_N_POINTS = 2000  # each extra sample costs one more Newton solve when plot_inputs/plot_outputs is set
EMT_DEFAULT_N_POINTS = 200  # -> default dt = t_final/199 when req.dt isn't given
MODAL_MAX_T_FINAL = 20.0  # free/step response are cheap (linear algebra, no ODE solve) -- a looser bound
EMT_SIMULATE_KWARGS = dict(rtol=1e-4, atol=1e-6, first_step=1e-8)
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


def bus_rows(result: PowerFlowResult) -> list[BusRow]:
    table = result.bus_table()
    return [BusRow(**row) for row in table.to_dict(orient="records")]


def powerflow_tables(result: PowerFlowResult) -> dict:
    """Shared by the powerflow and timeseries responses below -- both wrap
    the exact same six other-pandapower-tables extraction, one per
    snapshot for the latter.
    """
    return dict(
        lines=result.line_table().to_dict(orient="records"),
        transformers=result.trafo_table().to_dict(orient="records"),
        loads=result.load_table().to_dict(orient="records"),
        generators=result.gen_table().to_dict(orient="records"),
        static_generators=result.sgen_table().to_dict(orient="records"),
        external_grid=result.ext_grid_table().to_dict(orient="records"),
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
    return ValidateResponse(
        ok=not any(i.severity == "error" for i in issues),
        issues=[NetworkIssueRow(severity=i.severity, message=i.message, affects=list(i.affects)) for i in issues],
    )


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
    return ModalResponse(
        n_states=system.A.shape[0],
        stable=bool((modal.eigenvalues.real < 1e-6).all()),
        max_real_part=float(modal.eigenvalues.real.max()),
        modes=[ModeRow(**row) for row in table.to_dict(orient="records")],
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
    return StepResponseResponse(t=t.tolist(), y=series[output_names[0]], series=series)


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
        state_names=model.state_names, input_names=model.input_names, output_names=model.output_names
    )


def _resolve_emt_perturbation(
    model: NonlinearNetworkModel, req: EmtRequest
) -> tuple[np.ndarray, Callable[[float], np.ndarray] | None, str]:
    """``(x0, u_exo_fn, perturbed_name)`` for either perturbation kind --
    shared by :func:`emt_response` and :func:`emt_live_response` so the
    two never resolve "state" vs. "input" differently.
    """
    x0 = model.initial_state()
    u_exo_fn = None

    if req.perturb_kind == "state":
        idx = find_state_or_422(model, req.perturb_name)
        x0 = x0.copy()
        x0[idx] += req.perturb_offset
        perturbed_name = model.state_names[idx]
    elif req.perturb_kind == "input":
        # A permanent step in one exogenous reference, held from t=0 for the
        # whole run (not just an initial-condition offset) -- the system
        # starts at its own equilibrium and that equilibrium itself moves,
        # the more standard "P_ref step" kind of disturbance test. See
        # timedomain/emt.py's simulate()'s own u_exo_fn parameter, built for
        # exactly this and unused until now.
        idx = named_index_or_422(model.input_names, req.perturb_name, "input")
        perturbed_u_exo = model.default_u_exo()
        perturbed_u_exo[idx] += req.perturb_offset
        u_exo_fn = lambda t: perturbed_u_exo  # noqa: E731 -- held constant at the stepped value for all t
        perturbed_name = model.input_names[idx]
    else:
        raise HTTPException(status_code=422, detail="perturb_kind must be 'state' or 'input'")

    return x0, u_exo_fn, perturbed_name


def _validate_emt_t_final_and_n_points(req: EmtRequest) -> int:
    """Shared ``t_final``/``dt`` bounds check -- returns the resolved
    sample count (one-shot endpoint) / step-count safety bound (live
    endpoint), both keyed off the same ``EMT_MIN/MAX_N_POINTS`` constants.
    """
    if not (0 < req.t_final <= EMT_MAX_T_FINAL):
        raise HTTPException(status_code=422, detail=f"t_final must be in (0, {EMT_MAX_T_FINAL}]")
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
    model: NonlinearNetworkModel, t_pre: float, plot_states: list[str], plot_inputs: list[str], plot_outputs: list[str]
) -> PreSamples:
    """Two equilibrium samples, at ``t=-t_pre`` and ``t=0`` (the instant just
    before the disturbance), for prepending to a trajectory: a plot then shows
    the flat pre-disturbance x0 and a vertical jump/step at t=0.
    """
    x_eq, inputs_eq, outputs_eq = _equilibrium_signals(model)
    t = [-t_pre, 0.0]
    states = {n: [float(x_eq[model.state_names.index(n)])] * 2 for n in plot_states}
    inputs = {n: [inputs_eq[n]] * 2 for n in plot_inputs}
    outputs = {n: [outputs_eq[n]] * 2 for n in plot_outputs}
    return t, states, inputs, outputs


def _prepend(pre: PreSamples, t: list[float], series: dict, inputs: dict, outputs: dict) -> PreSamples:
    pre_t, pre_s, pre_i, pre_o = pre
    return (
        pre_t + list(t),
        {n: pre_s.get(n, []) + list(v) for n, v in series.items()},
        {n: pre_i.get(n, []) + list(v) for n, v in inputs.items()},
        {n: pre_o.get(n, []) + list(v) for n, v in outputs.items()},
    )


def linear_overlay(
    network: Network, model: NonlinearNetworkModel, req: EmtRequest, t: np.ndarray,
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
    system, _ = build_modal_from_network(network)
    x_eq, inputs_eq, outputs_eq = _equilibrium_signals(model)
    n_x, n_u = system.A.shape[0], system.B.shape[1]
    if len(x_eq) != n_x or list(system.state_names) != list(model.state_names):
        raise HTTPException(
            status_code=501, detail="linearized and nonlinear models don't share a state vector for this network"
        )

    t = np.asarray(t, dtype=float)
    dx0 = np.zeros(n_x)
    du = np.zeros(n_u)
    if req.perturb_kind == "state":
        dx0[model.state_names.index(req.perturb_name)] = req.perturb_offset
    else:
        du[named_index_or_422(list(system.input_names), req.perturb_name, "input")] = req.perturb_offset

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


def emt_response(network: Network, req: EmtRequest) -> EmtResponse:
    n_points = _validate_emt_t_final_and_n_points(req)
    model = build_nonlinear_model_from_network(network)
    x0, u_exo_fn, perturbed_name = _resolve_emt_perturbation(model, req)

    t_eval = np.linspace(0.0, req.t_final, n_points)
    try:
        sim = simulate(model, (0.0, req.t_final), x0=x0, u_exo_fn=u_exo_fn, t_eval=t_eval, **EMT_SIMULATE_KWARGS)
    except RuntimeError as e:
        # The coupled Newton solve (see timedomain/emt.py) can fail to
        # converge for a large enough perturbation -- a trajectory that
        # fails mid-integration has no partial result worth returning, so
        # this surfaces as a clean 422 instead of an unhandled 500.
        raise HTTPException(
            status_code=422,
            detail=f"nonlinear solver did not converge for this perturbation ({e}); try a smaller offset",
        ) from e

    if req.plot_states:
        unknown = set(req.plot_states) - set(model.state_names)
        if unknown:
            raise HTTPException(status_code=422, detail=f"unknown state name(s): {sorted(unknown)}")
        plot_states = req.plot_states
    else:
        plot_states = [n for n in model.state_names if "dw_r" in n] or model.state_names
    series = {n: sim.x[model.state_names.index(n), :].tolist() for n in plot_states}

    # Inputs/outputs are opt-in (empty list = skip) -- recovering them costs
    # about as much again as the integration itself (see
    # recover_inputs_and_outputs's docstring), and inputs are held constant
    # (or, for perturb_kind="input", stepped-then-constant) throughout a
    # run, so there's no reason to pay for either by default.
    inputs: dict[str, list[float]] = {}
    outputs: dict[str, list[float]] = {}
    if req.plot_inputs or req.plot_outputs:
        # Reuse the *same* u_exo_fn the trajectory was actually integrated
        # with -- otherwise a perturb_kind="input" run would recover the
        # unperturbed exogenous inputs, silently contradicting the
        # trajectory that was actually simulated.
        all_inputs, all_outputs = model.recover_inputs_and_outputs(sim.t, sim.x, u_exo_fn=u_exo_fn)
        if req.plot_inputs:
            unknown = set(req.plot_inputs) - set(all_inputs)
            if unknown:
                raise HTTPException(status_code=422, detail=f"unknown input name(s): {sorted(unknown)}")
            inputs = {n: all_inputs[n].tolist() for n in req.plot_inputs}
        if req.plot_outputs:
            unknown = set(req.plot_outputs) - set(all_outputs)
            if unknown:
                raise HTTPException(status_code=422, detail=f"unknown output name(s): {sorted(unknown)}")
            outputs = {n: all_outputs[n].tolist() for n in req.plot_outputs}

    actual_dt = float(sim.t[1] - sim.t[0]) if len(sim.t) > 1 else req.t_final
    linear = None
    if req.linear_overlay:
        linear = linear_overlay(network, model, req, sim.t, list(series), list(inputs), list(outputs))
    t_out = sim.t.tolist()
    if req.t_pre > 0:
        pre = _pre_disturbance_samples(model, req.t_pre, list(series), list(inputs), list(outputs))
        t_out, series, inputs, outputs = _prepend(pre, t_out, series, inputs, outputs)
        if linear is not None:
            lt, ls, li, lo = _prepend(pre, linear.t, linear.series, linear.inputs, linear.outputs)
            linear = LinearOverlay(t=lt, series=ls, inputs=li, outputs=lo)
    return EmtResponse(
        perturbed=perturbed_name, perturb_kind=req.perturb_kind, dt=actual_dt,
        state_names=list(series.keys()), t=t_out, series=series, inputs=inputs, outputs=outputs, linear=linear,
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
        self, model: NonlinearNetworkModel, x0: np.ndarray, u_exo_fn: Callable[[float], np.ndarray] | None,
        perturbed_name: str, plot_states: list[str], plot_inputs: list[str], plot_outputs: list[str],
        max_step: float, t_final: float, network: Network | None = None, req: EmtRequest | None = None,
    ) -> None:
        self.model, self.x0, self.u_exo_fn = model, x0, u_exo_fn
        self.perturbed_name = perturbed_name
        self.plot_states, self.plot_inputs, self.plot_outputs = plot_states, plot_inputs, plot_outputs
        self.max_step, self.t_final = max_step, t_final
        self.network, self.req = network, req


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
    x0, u_exo_fn, perturbed_name = _resolve_emt_perturbation(model, req)

    if req.plot_states:
        unknown = set(req.plot_states) - set(model.state_names)
        if unknown:
            raise HTTPException(status_code=422, detail=f"unknown state name(s): {sorted(unknown)}")
        plot_states = req.plot_states
    else:
        plot_states = [n for n in model.state_names if "dw_r" in n] or model.state_names

    plot_inputs: list[str] = []
    if req.plot_inputs:
        unknown = set(req.plot_inputs) - set(model.input_names)
        if unknown:
            raise HTTPException(status_code=422, detail=f"unknown input name(s): {sorted(unknown)}")
        plot_inputs = req.plot_inputs

    plot_outputs: list[str] = []
    if req.plot_outputs:
        unknown = set(req.plot_outputs) - set(model.output_names)
        if unknown:
            raise HTTPException(status_code=422, detail=f"unknown output name(s): {sorted(unknown)}")
        plot_outputs = req.plot_outputs

    max_step = req.dt if req.dt else req.t_final / EMT_DEFAULT_N_POINTS
    return _EmtLivePlan(
        model, x0, u_exo_fn, perturbed_name, plot_states, plot_inputs, plot_outputs, max_step, req.t_final,
        network=network, req=req,
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
    # Resolved once, not per step: state_names.index() is an O(n) scan and
    # this loop can run thousands of times.
    state_idx = [plan.model.state_names.index(n) for n in plan.plot_states]

    # Pre-disturbance equilibrium samples (t_pre > 0), streamed first -- same
    # line shape as a solver step, so the client plots them like any other.
    t_pre = plan.req.t_pre if plan.req is not None else 0.0
    pre = None
    if t_pre > 0:
        pre = _pre_disturbance_samples(plan.model, t_pre, plan.plot_states, plan.plot_inputs, plan.plot_outputs)
        pre_t, pre_s, pre_i, pre_o = pre
        for k, t_k in enumerate(pre_t):
            yield json.dumps({
                "t": t_k,
                "states": {n: v[k] for n, v in pre_s.items()},
                "inputs": {n: v[k] for n, v in pre_i.items()},
                "outputs": {n: v[k] for n, v in pre_o.items()},
            }) + "\n"

    n_steps = 0
    try:
        for step in simulate_steps(
            plan.model, (0.0, plan.t_final), x0=plan.x0, u_exo_fn=plan.u_exo_fn,
            max_step=plan.max_step, **EMT_SIMULATE_KWARGS,
        ):
            if await request.is_disconnected():
                return  # client hit "Stop" / navigated away -- stop computing, not an error
            n_steps += 1
            if n_steps > EMT_MAX_N_POINTS:
                yield json.dumps({"error": f"exceeded {EMT_MAX_N_POINTS} steps; raise dt to see fewer, coarser updates"}) + "\n"
                return
            yield json.dumps({
                "t": step.t,
                "states": dict(zip(plan.plot_states, (float(step.x[i]) for i in state_idx))),
                "inputs": {n: step.inputs[n] for n in plan.plot_inputs},
                "outputs": {n: step.outputs[n] for n in plan.plot_outputs},
            }) + "\n"
    except RuntimeError as e:
        yield json.dumps({"error": f"nonlinear solver did not converge for this perturbation ({e}); try a smaller offset"}) + "\n"
        return
    done: dict = {"done": True, "perturbed": plan.perturbed_name, "perturb_kind": perturb_kind, "n_steps": n_steps}
    if plan.req is not None and plan.req.linear_overlay and plan.network is not None:
        # Computed once the nonlinear run has finished, on an even grid over
        # the same span (the solver's own adaptive step times are too
        # irregular to be worth reproducing for a linear overlay).
        try:
            grid = np.linspace(0.0, plan.t_final, EMT_DEFAULT_N_POINTS)
            lin = linear_overlay(
                plan.network, plan.model, plan.req, grid, plan.plot_states, plan.plot_inputs, plan.plot_outputs
            )
            if pre is not None:
                lt, ls, li, lo = _prepend(pre, lin.t, lin.series, lin.inputs, lin.outputs)
                lin = LinearOverlay(t=lt, series=ls, inputs=li, outputs=lo)
            done["linear"] = lin.model_dump()
        except HTTPException as e:
            done["linear_error"] = str(e.detail)
    yield json.dumps(done) + "\n"
