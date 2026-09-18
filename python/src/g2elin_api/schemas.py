"""HTTP response models. Kept separate from ``g2elin_core``'s own types
(``PowerFlowResult`` wraps a live pandapower net, ``ModalAnalysisResult``
holds numpy arrays — neither is JSON-serializable directly) rather than
serializing those internal types, so the wire format is a deliberate
choice, not an accident of what happened to be easy to `dict()`.
"""

from __future__ import annotations

from pydantic import BaseModel

from g2elin_core.network.schema import Network


class BusRow(BaseModel):
    bus: int
    vm_pu: float
    va_degree: float
    p_net_gen_mw: float
    q_net_gen_mvar: float


class PowerFlowOptions(BaseModel):
    """Solver settings forwarded to ``pandapower.runpp`` -- see
    ``analysis.POWERFLOW_ALGORITHMS`` for the accepted ``algorithm`` values.
    Every field's default is pandapower's own default, so an empty
    ``PowerFlowOptions()`` solves exactly as before these options existed.
    """

    algorithm: str = "nr"
    max_iteration: int | str = "auto"
    tolerance_mva: float = 1e-8
    init: str = "auto"


class PowerFlowRequest(BaseModel):
    options: PowerFlowOptions = PowerFlowOptions()


class PowerFlowResponse(BaseModel):
    converged: bool
    buses: list[BusRow]
    total_losses_mw: float
    # The rest of what pandapower computes, beyond the per-bus table above --
    # loosely typed (dict, not a row model per table) since each table's
    # column set genuinely differs (res_line has 13 columns, res_load has 2)
    # and this is an inspection/dropdown feature, not something else in this
    # codebase consumes programmatically. Empty (not omitted) when a preset
    # has no elements of that kind, e.g. no res_sgen without an sgen.
    lines: list[dict]
    transformers: list[dict]
    loads: list[dict]
    generators: list[dict]
    static_generators: list[dict]
    external_grid: list[dict]
    # Solver diagnostics, when pandapower reports them (None otherwise, e.g.
    # on non-convergence or an algorithm that doesn't record them).
    algorithm: str = "nr"
    iterations: int | None = None
    solve_time_s: float | None = None


class ModeRow(BaseModel):
    mode: int
    real: float
    imag: float
    undamped_hz: float
    damped_hz: float
    damping_pct: float
    state1: str
    part1_pct: float
    state2: str
    part2_pct: float
    state3: str
    part3_pct: float


class ModalResponse(BaseModel):
    n_states: int
    stable: bool
    max_real_part: float
    modes: list[ModeRow]
    state_names: list[str]
    input_names: list[str]
    output_names: list[str]
    participation: list[list[float]]  # [state][mode], each column sums to 1 -- see ModalAnalysisResult


class SensitivityRequest(BaseModel):
    mode: int


class SensitivityEntryRow(BaseModel):
    row_state: str
    col_state: str
    value: float


class SensitivityResponse(BaseModel):
    mode: int
    state_names: list[str]
    matrix: list[list[float]]
    top: list[SensitivityEntryRow]


class ModeShapeRequest(BaseModel):
    mode: int


class ModeShapeResponse(BaseModel):
    mode: int
    states: list[str]
    angles_deg: list[float]


class FreeResponseRequest(BaseModel):
    perturb_state: str
    offset: float = 0.05
    t_final: float = 2.0
    state_filter: str = "dw_r"
    # Exact state names to return -- when non-empty this replaces
    # state_filter (the web UI's per-subplot pickers send exact names).
    plot_states: list[str] = []


class FreeResponseResponse(BaseModel):
    perturb_state: str
    t: list[float]
    series: dict[str, list[float]]


class StepResponseRequest(BaseModel):
    input_name: str
    # One output (output_name) or several (output_names) -- the latter is a
    # single SIMO step response, one simulation for every output at once.
    output_name: str | None = None
    output_names: list[str] = []
    amplitude: float = 0.1
    t_final: float = 2.0


class StepResponseResponse(BaseModel):
    t: list[float]
    y: list[float]  # the first requested output, kept for single-output callers
    series: dict[str, list[float]] = {}  # every requested output, by name


class TimeSeriesSnapshot(BaseModel):
    label: str
    converged: bool
    buses: list[BusRow]
    total_losses_mw: float
    # Same other-pandapower-tables treatment as PowerFlowResponse, per
    # snapshot -- see that schema's own comment for why these are loosely
    # typed (list[dict], not a row model per table).
    lines: list[dict]
    transformers: list[dict]
    loads: list[dict]
    generators: list[dict]
    static_generators: list[dict]
    external_grid: list[dict]


class TimeSeriesResponse(BaseModel):
    snapshots: list[TimeSeriesSnapshot]


class BatchPowerFlowRequest(PowerFlowRequest):
    """A ramp of power-flow snapshots from the base operating point (every
    factor 1.0) to the target scale factors below, ``steps`` equal
    increments, base point included -- ``steps + 1`` snapshots in total.
    ``der_scale`` is keyed by DER unit type (``"sm"``/``"gfm"``/``"gfl"``)
    and scales that type's P *and* Q setpoints; a type left out stays at 1.0.
    """

    load_p_scale: float = 1.0
    load_q_scale: float = 1.0
    der_scale: dict[str, float] = {}
    steps: int = 10


class BatchSnapshot(TimeSeriesSnapshot):
    load_p_scale: float
    load_q_scale: float
    der_scale: dict[str, float]


class BatchPowerFlowResponse(BaseModel):
    snapshots: list[BatchSnapshot]


class StatesResponse(BaseModel):
    state_names: list[str]
    input_names: list[str]
    output_names: list[str]


class BusNodeRow(BaseModel):
    id: int
    name: str
    vn_kv: float
    x: float
    y: float
    unit_type: str | None
    load_p_mw: float | None
    der_info: dict | None = None


class TopologyEdgeRow(BaseModel):
    kind: str
    from_bus: int
    to_bus: int
    name: str
    params: dict = {}


class TopologyResponse(BaseModel):
    nodes: list[BusNodeRow]
    edges: list[TopologyEdgeRow]


class EmtRequest(BaseModel):
    # "state": an initial-condition offset (x0[idx] += perturb_offset) --
    # the system starts away from equilibrium and (maybe) settles back.
    # "input": a permanent step in one exogenous reference from t=0 onward
    # (u_exo[idx] += perturb_offset, held for the whole run) -- the system
    # starts *at* equilibrium and the equilibrium itself moves, the more
    # standard "P_ref step" kind of disturbance test. perturb_name is
    # looked up against state_names or input_names accordingly (exact name
    # from GET .../states, not a substring).
    perturb_kind: str = "state"  # "state" | "input"
    perturb_name: str
    perturb_offset: float = 0.02  # the perturbation's amplitude -- larger can push the coupled Newton solve past convergence, see main.py
    t_final: float = 1.0
    dt: float | None = None  # None = 200 samples over [0, t_final] (unchanged default); else t_final/dt (+1) samples, bounded server-side
    # Exact names from GET .../states, not substrings -- lets the frontend
    # use real multi-select pickers over the full list rather than a text
    # filter. Empty plot_states defaults (server-side) to every dw_r_* state
    # (speed deviations -- the most broadly interpretable signal); empty
    # plot_inputs/plot_outputs means "none" (they cost an extra Newton solve
    # per sample to recover -- see main.py -- so they're opt-in, not
    # defaulted the way plot_states is).
    plot_states: list[str] = []
    plot_inputs: list[str] = []
    plot_outputs: list[str] = []
    # Seconds of undisturbed equilibrium prepended before the disturbance at
    # t=0 (samples at t=-t_pre and t=0-), so a plot shows x0 before the
    # jump/step. 0 keeps the trajectory starting at t=0 exactly.
    t_pre: float = 0.0
    # Also return the linearized model's response to the same disturbance
    # (same equilibrium, same state/input/output names), for overlaying.
    linear_overlay: bool = False


class LinearOverlay(BaseModel):
    t: list[float]
    series: dict[str, list[float]]
    inputs: dict[str, list[float]]
    outputs: dict[str, list[float]]


class EmtResponse(BaseModel):
    perturbed: str  # the resolved exact name (== perturb_name)
    perturb_kind: str
    dt: float  # the actual sample spacing used (t[1] - t[0]) -- answers "what's the current timestep"
    state_names: list[str]  # names actually plotted (== plot_states, or the dw_r_* default if it was empty)
    t: list[float]
    series: dict[str, list[float]]
    inputs: dict[str, list[float]]
    outputs: dict[str, list[float]]
    linear: LinearOverlay | None = None


class PresetSummary(BaseModel):
    id: str
    name: str
    description: str
    n_buses: int
    n_lines: int
    n_transformers: int
    n_loads: int
    n_der_units: int
    der_unit_types: list[str]


class ErrorResponse(BaseModel):
    detail: str


# --- /api/network/* request bodies -----------------------------------------
# Each subclasses its preset-side counterpart, adding just the Network
# itself, rather than re-declaring every field -- see network_routes.py.


class NetworkRequest(BaseModel):
    """Body for the network-based endpoints that take no other input today
    (topology, powerflow, modal, timeseries, states)."""

    network: Network


class NetworkPowerFlowRequest(PowerFlowRequest):
    network: Network


class NetworkBatchPowerFlowRequest(BatchPowerFlowRequest):
    network: Network


class NetworkSensitivityRequest(SensitivityRequest):
    network: Network


class NetworkModeShapeRequest(ModeShapeRequest):
    network: Network


class NetworkFreeResponseRequest(FreeResponseRequest):
    network: Network


class NetworkStepResponseRequest(StepResponseRequest):
    network: Network


class NetworkEmtRequest(EmtRequest):
    network: Network


class NetworkIssueRow(BaseModel):
    severity: str  # "error" | "warning"
    message: str
    affects: list[str]  # which of "powerflow"/"modal"/"emt" this issue affects


class ValidateResponse(BaseModel):
    ok: bool  # true iff there are no "error"-severity issues (warnings don't block anything)
    issues: list[NetworkIssueRow]
