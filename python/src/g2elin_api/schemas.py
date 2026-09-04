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


class FreeResponseResponse(BaseModel):
    perturb_state: str
    t: list[float]
    series: dict[str, list[float]]


class StepResponseRequest(BaseModel):
    input_name: str
    output_name: str
    amplitude: float = 0.1
    t_final: float = 2.0


class StepResponseResponse(BaseModel):
    t: list[float]
    y: list[float]


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


class EmtResponse(BaseModel):
    perturbed: str  # the resolved exact name (== perturb_name)
    perturb_kind: str
    dt: float  # the actual sample spacing used (t[1] - t[0]) -- answers "what's the current timestep"
    state_names: list[str]  # names actually plotted (== plot_states, or the dw_r_* default if it was empty)
    t: list[float]
    series: dict[str, list[float]]
    inputs: dict[str, list[float]]
    outputs: dict[str, list[float]]


class RoaRequest(BaseModel):
    axis_x_state: str
    axis_x_range: float = 0.4
    axis_y_state: str
    axis_y_range: float = 0.02
    grid_n: int = 3
    t_final: float = 1.0
    t_early: float = 0.1


class RoaResponse(BaseModel):
    axis_x_label: str
    axis_x_offsets: list[float]
    axis_y_label: str
    axis_y_offsets: list[float]
    in_roa: list[list[bool]]
    failed: list[list[bool]]
    early_distance: list[list[float | None]]  # None where the Newton solve failed (not a real distance)
    late_distance: list[list[float | None]]


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


class NetworkRoaRequest(RoaRequest):
    network: Network


class NetworkIssueRow(BaseModel):
    severity: str  # "error" | "warning"
    message: str
    affects: list[str]  # which of "powerflow"/"modal"/"emt"/"roa" this issue affects


class ValidateResponse(BaseModel):
    ok: bool  # true iff there are no "error"-severity issues (warnings don't block anything)
    issues: list[NetworkIssueRow]
