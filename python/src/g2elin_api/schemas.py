"""HTTP response models. Kept separate from ``g2elin_core``'s own types
(``PowerFlowResult`` wraps a live pandapower net, ``ModalAnalysisResult``
holds numpy arrays — neither is JSON-serializable directly) rather than
serializing those internal types, so the wire format is a deliberate
choice, not an accident of what happened to be easy to `dict()`.
"""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field

from g2elin_core.network.schema import Network


class BusRow(BaseModel):
    bus: int
    # None on a bus de-energized by open breakers (see network.breakers)
    vm_pu: float | None
    va_degree: float | None
    p_net_gen_mw: float | None
    q_net_gen_mvar: float | None


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


class ModeCategoryInfo(BaseModel):
    """One of the kinds a mode can be (g2elin_core.reduction.CATEGORIES)."""

    id: str
    label: str
    note: str


class ModeRow(BaseModel):
    mode: int
    real: float
    imag: float
    undamped_hz: float
    damped_hz: float
    damping_pct: float
    # What kind of mode this is, from which states participate in it
    # (g2elin_core.modal.classify). "share" is how much of the mode that
    # category accounts for; "shares" gives every category's, so the map can
    # show the full breakdown rather than only the verdict.
    category: str = "mixed"
    category_share: float = 0.0
    category_shares: dict[str, float] = {}
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
    # The catalogue behind ModeRow.category, so the UI labels and explains
    # the groups without hard-coding either.
    categories: list[ModeCategoryInfo] = []
    # Modes that are only the model's free reference angle(s) (see
    # modal.reference_angle_modes): marginal by construction, so they are left
    # out of "stable"/"max_real_part".
    reference_modes: list[int] = []


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


class ParameterSensitivityRequest(BaseModel):
    mode: int
    # Narrow the scan; empty means every unit / every parameter.
    units: list[int] = []
    parameters: list[str] = []
    n_entries: int = Field(default=8, ge=1, le=40)


class ParameterEffectRow(BaseModel):
    unit: int
    unit_label: str
    parameter: str
    value: float
    d_lambda_real: float
    d_lambda_imag: float
    # What a 1% increase in the parameter does to the mode.
    d_freq_hz: float
    d_damping_pct: float
    magnitude: float


class EntryParametersRow(BaseModel):
    """One high-sensitivity entry of A, and the parameters it is built from."""

    row_state: str
    col_state: str
    sensitivity: float
    parameters: list[str]


class ParameterSensitivityResponse(BaseModel):
    mode: int
    eigenvalue_real: float
    eigenvalue_imag: float
    effects: list[ParameterEffectRow]
    entries: list[EntryParametersRow]
    notes: list[str] = []


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


class MeasurementInfo(BaseModel):
    name: str   # e.g. "V_{bus4}"
    group: str  # the element, e.g. "Bus 4 (bus4)"
    label: str  # e.g. "Voltage magnitude"
    unit: str   # e.g. "pu", "deg", "Hz"


class StatesResponse(BaseModel):
    state_names: list[str]
    input_names: list[str]
    output_names: list[str]
    # Measurement outputs (g2elin_core.timedomain.measurements) -- plottable
    # via EmtRequest.plot_measurements.
    measurements: list[MeasurementInfo] = []


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


class NetworkEventSpec(BaseModel):
    """A network event applied at t = 0 (timedomain.events): a breaker
    opening, a load step, or a phase jump."""

    kind: str  # "breaker" | "load_step" | "phase_jump"
    element: str = ""  # breaker: "line" | "transformer" | "load" | "unit"
    index: int = 0  # line/transformer/load index, or the unit's id
    dp_pct: float = 0.0  # load_step: active-power change, % of the load
    dq_pct: float = 0.0  # load_step: reactive-power change, %
    bus: int | None = None  # phase_jump: bus id (a network bus, or the infinite bus's own bus)
    angle_deg: float = 0.0  # phase_jump


class EmtRequest(BaseModel):
    # "state": an initial-condition offset (x0[idx] += perturb_offset) --
    # the system starts away from equilibrium and (maybe) settles back.
    # "input": a permanent step in one exogenous reference from t=0 onward
    # (u_exo[idx] += perturb_offset, held for the whole run) -- the system
    # starts *at* equilibrium and the equilibrium itself moves, the more
    # standard "P_ref step" kind of disturbance test. perturb_name is
    # looked up against state_names or input_names accordingly (exact name
    # from GET .../states, not a substring).
    # "event": a network event (``event``) -- breaker opening, load step or
    # phase jump -- at t=0; perturb_name/perturb_offset are then unused.
    perturb_kind: str = "state"  # "state" | "input" | "event"
    perturb_name: str = ""
    event: NetworkEventSpec | None = None
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
    # Measurement outputs by name (GET .../states lists them): power flows,
    # bus voltage/angle/frequency, 3-phase voltages, unit frequencies.
    plot_measurements: list[str] = []
    # Seconds of undisturbed equilibrium prepended before the disturbance at
    # t=0 (samples at t=-t_pre and t=0-), so a plot shows x0 before the
    # jump/step. 0 keeps the trajectory starting at t=0 exactly.
    t_pre: float = 0.0
    # Also return the linearized model's response to the same disturbance
    # (same equilibrium, same state/input/output names), for overlaying.
    linear_overlay: bool = False

    # --- solver ---------------------------------------------------------------
    # How the trajectory is integrated. The defaults are what this tool used
    # before any of these were settable, so leaving them alone changes
    # nothing; see g2elin_core.timedomain.emt.SOLVERS for what each one is
    # good for, and GET /api/solvers for the same catalogue as data.
    stepping: Literal["variable", "fixed"] = "variable"
    solver: str = "Radau"  # variable stepping only
    rtol: float = Field(default=1e-4, gt=0, le=1e-1)
    atol: float = Field(default=1e-6, gt=0, le=1e-1)
    # Largest step the adaptive solver may take. None = unbounded, which is
    # what it should normally be: capping it makes the solver take steps it
    # didn't need, and the sample spacing is `dt`'s job, not this one's.
    max_step: float | None = Field(default=None, gt=0)
    # Fixed stepping only: the step held for the whole run. None picks the
    # output spacing, so "fixed step" means "one step per plotted point".
    fixed_step: float | None = Field(default=None, gt=0)


class SolverInfo(BaseModel):
    id: str
    label: str
    note: str
    implicit: bool


class SolversResponse(BaseModel):
    """The integrators the time-domain page offers, and the defaults."""

    solvers: list[SolverInfo]
    default: str
    default_rtol: float
    default_atol: float


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
    series: dict[str, list[float | None]]
    inputs: dict[str, list[float | None]]
    outputs: dict[str, list[float | None]]
    # Signals of elements an event removes read null after t=0 (powers read 0).
    measurements: dict[str, list[float | None]] = {}
    linear: LinearOverlay | None = None
    linear_note: str | None = None  # why there's no linear overlay (events that change the model)


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


class NetworkParameterSensitivityRequest(ParameterSensitivityRequest):
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


class SweepTarget(BaseModel):
    """One network parameter to sweep. ``element`` is "network", "bus",
    "line", "transformer", "load" or "unit"; ``key`` identifies the element
    (bus/unit id, or list index for lines/transformers/loads; unused for
    "network"); ``field`` is the element's field name, or ``params.<name>``
    for a unit's control/electrical parameter (see DerUnit.params).
    """

    element: str
    key: int | None = None
    field: str


class SweepExtra(BaseModel):
    """A further parameter varied together with the main one: it goes from
    ``start`` to ``stop`` in lockstep with the main parameter's progress, so
    every step of the sweep is the combined effect of all of them."""

    target: SweepTarget
    start: float
    stop: float


class SweepRequest(BaseModel):
    target: SweepTarget
    start: float
    stop: float
    step: float
    extra: list[SweepExtra] = []
    # Also stream each mode's most participating states at every value (they
    # come out of the same eigendecomposition, so they cost nothing to
    # compute -- only to send).
    participation: bool = True


class NetworkSweepRequest(SweepRequest):
    network: Network


class UnitDefaultsRequest(BaseModel):
    """Everything a unit's default parameter set depends on -- no network
    needed, so the editor can show a unit's parameters before it's wired up.
    ``rt_pu``/``lt_pu`` are the network's first transformer (the models'
    shared-transformer convention, see operating_point.py)."""

    unit_type: str
    sn_mva: float
    f_hz: float
    un_kv: float
    rt_pu: float = 0.0
    lt_pu: float = 0.05


class UnitDefaultsResponse(BaseModel):
    params: dict[str, float]  # empty for a unit type without its own parameter set (infinite bus)


class NetworkIssueRow(BaseModel):
    severity: str  # "error" | "warning"
    message: str
    affects: list[str]  # which of "powerflow"/"modal"/"emt" this issue affects


class IslandRow(BaseModel):
    buses: list[int]
    reference: int | None  # the unit acting as this island's slack, None = blacked out


class ServiceInfo(BaseModel):
    """What open breakers leave in service (network.breakers.service_state)."""

    references: list[int]  # the units acting as a power-flow reference, one per energized island
    islands: list[IslandRow]
    energized_buses: list[int]
    lines: list[bool]
    transformers: list[bool]
    loads: list[bool]
    der_units: dict[int, bool]


class ValidateResponse(BaseModel):
    ok: bool  # true iff there are no "error"-severity issues (warnings don't block anything)
    issues: list[NetworkIssueRow]
    service: ServiceInfo | None = None


# --- Model-order reduction ----------------------------------------------------


class StateGroupInfo(BaseModel):
    """One switchable approximation of an element type (g2elin_core.reduction)."""

    id: str
    label: str
    states: list[str]  # the state names this group covers, for display
    allowed: list[str]  # which of dynamic/algebraic/frozen this group accepts
    default: str
    locked: bool  # true = dynamic only, the group can never be removed
    # Groups that must be algebraic before this one can be (a control loop
    # needs what it regulates to be an unknown). The UI carries these along
    # rather than letting the user pick a combination the server refuses.
    requires: list[str] = []
    note: str


class ModelLevelInfo(BaseModel):
    id: str
    label: str
    note: str
    modes: dict[str, str]  # group id -> mode


class ElementModelInfo(BaseModel):
    kind: str  # "network" | "sm" | "gfm" | "gfl"
    label: str
    default_level: str
    levels: list[ModelLevelInfo]
    groups: list[StateGroupInfo]


class ModelLevelsResponse(BaseModel):
    """The whole reduction catalogue -- what the web UI builds its level
    pickers and per-state-group controls from, so nothing about the
    available choices is hard-coded in the frontend."""

    elements: list[ElementModelInfo]


class UnitModelRow(BaseModel):
    """One unit's resolved model settings, after network defaults and the
    unit's own overrides."""

    id: int
    unit_type: str
    label: str
    level: str | None  # None = a custom combination matching no named level
    modes: dict[str, str]  # group id -> mode
    n_states: int


class ModelSummaryResponse(BaseModel):
    """What the current settings add up to: the model class the run will
    produce, and the size it will be."""

    model_class: str  # "EMT" | "RMS" | "Mixed"
    network_level: str | None
    network_modes: dict[str, str]
    network_frequency: str
    units: list[UnitModelRow]
    n_states: int
    n_states_full: int


class StateRiskRow(BaseModel):
    state: str
    group: str
    mode_hz: float
    mode_damping_pct: float
    participation: float


class ModePairRow(BaseModel):
    full_hz: float
    full_damping_pct: float
    reduced_hz: float | None
    reduced_damping_pct: float | None
    matched: bool
    # An unmatched mode that was made of the removed states themselves --
    # the reduction doing what was asked, not losing something.
    expected_loss: bool = False
    removed_share: float = 0.0
    d_hz: float | None
    d_damping_pct: float | None


class AdequacyRequest(BaseModel):
    network: Network
    band_hz: float = Field(
        default=5.0, gt=0, le=1000,
        description="The frequency band the reduced model is expected to reproduce. "
        "Electromechanical studies use a few Hz; converter interaction studies need more.",
    )


class AdequacyResponse(BaseModel):
    """Whether this network's chosen reduction is safe for this network --
    see g2elin_core.modal.adequacy."""

    verdict: str  # "safe" | "check" | "unsafe" | "full_order"
    band_hz: float
    n_states_full: int
    n_states_reduced: int
    removed_states: list[str]
    risks: list[StateRiskRow]
    modes: list[ModePairRow]
    max_d_hz: float
    max_d_damping_pct: float
    unmatched: int
    stiffness_full: float
    stiffness_reduced: float
    notes: list[str]
