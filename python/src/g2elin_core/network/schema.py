"""Typed network description for G2ELin.

Replaces the MATLAB toolbox's numeric-column matrix convention
(``Y_network`` / ``Y_line`` / ``Y_TR`` / ``Y_DER``) with named, validated
fields. See ``legacy.py`` for a converter from the original convention,
used to import existing G2ELin cases.

Per-unit convention: everything is in system per-unit on ``Network.sn_mva``,
matching the original tool (``Power_Fl.m`` / ``Load_flow.m``). Bus nominal
voltages (``Bus.vn_kv``) are physical, used only for display and for
transformer turns-ratio calculations.
"""

from __future__ import annotations

from enum import Enum
from typing import Literal

from pydantic import BaseModel, Field, PrivateAttr, model_validator

from g2elin_core import reduction


class BusType(str, Enum):
    """Load-flow bus type, matching ``Y_network`` column 4 / ``Y_DER`` column 4."""

    SLACK = "slack"  # type 1
    PV = "pv"  # type 2
    PQ = "pq"  # type 3


class UnitType(str, Enum):
    """DER unit type, matching ``Y_DER`` column 2 (``Y_network`` column 2 for raw nodes)."""

    NONE = "none"  # 0 - plain network node
    INFINITE_BUS = "infinite_bus"  # 1
    GFM = "gfm"  # 2 - grid-forming converter
    GFL = "gfl"  # 3 - grid-following converter
    SYNCHRONOUS_MACHINE = "sm"  # 4


class GfmController(str, Enum):
    """Outer power-control law for a GFM unit, matching ``GFM_P_control`` in script_generic.m."""

    DROOP = "droop"
    DROOP_FILTERED = "droop_filtered"
    DVOC = "dvoc"
    VSM = "vsm"
    MATCHING = "matching"


class ExciterModel(str, Enum):
    """Which AVR/exciter a synchronous machine carries.

    See ``docs/sphinx/design/controller-models.md`` for the equations and for
    how these relate to the IEEE 421.5 library.
    """

    #: The regulator this tool has always had: terminal-voltage transducer,
    #: first-order amplifier, first-order exciter and a rate feedback.
    G2ELIN = "g2elin"
    #: Kundur Fig. E12.9: a thyristor exciter -- transducer, gain and
    #: transient gain reduction, with no exciter lag. IEEE ST1A family.
    KUNDUR = "kundur"


class PssModel(str, Enum):
    """Which power system stabiliser a synchronous machine carries, if any."""

    #: The stabiliser this tool has always had: an input low-pass, a washout
    #: and two lead-lag stages.
    G2ELIN = "g2elin"
    #: Kundur Fig. E12.9: gain, washout and two lead-lag stages, taking the
    #: speed deviation unfiltered.
    KUNDUR = "kundur"
    #: No stabiliser fitted.
    NONE = "none"


class GovernorModel(str, Enum):
    """Whether a synchronous machine's prime mover responds to speed."""

    #: Droop into a first-order lag, the governor this tool has always had.
    G2ELIN = "g2elin"
    #: No governor: constant mechanical power, which is what most textbook
    #: small-signal examples assume.
    NONE = "none"


class Bus(BaseModel):
    """A network node. Corresponds to one row of ``Y_network``."""

    id: int
    name: str = ""
    vn_kv: float = Field(gt=0, description="Nominal line-to-line voltage (kV)")


class Line(BaseModel):
    """A pi-model transmission line, in per unit on the network's ``sn_mva`` base.

    Corresponds to one row of ``Y_line`` (excluding transformer rows, which
    ``network_form.m`` appends to the same matrix but are kept separate here
    as :class:`Transformer`).
    """

    from_bus: int
    to_bus: int
    r_pu: float = Field(ge=0)
    x_pu: float = Field(gt=0)
    b_pu: float = Field(ge=0, description="Total shunt susceptance (charging), per unit")
    length_km: float = Field(default=1.0, gt=0)
    name: str = ""
    from_closed: bool = Field(default=True, description="Breaker at the from-bus end (open = line out of service)")
    to_closed: bool = Field(default=True, description="Breaker at the to-bus end (open = line out of service)")


class Transformer(BaseModel):
    """A two-winding transformer with leakage impedance only (no magnetizing branch).

    Matches ``Y_TR``'s simple R+jX model, built from a DER unit's
    ``TR_R`` / ``TR_XL`` columns.

    Two quite different roles share this one element, told apart by whether a
    unit sits on the LV bus (``network/breakers.unit_transformers``):

    * a **unit step-up**, whose LV bus carries a DER. Its impedance is part of
      that unit's own model -- the classic machine-behind-transformer-
      impedance -- and its LV bus never appears in the dynamic model.
    * a **branch transformer** between two grid buses, which is a branch of
      the network like a line, and gets a block of its own.
    """

    hv_bus: int
    lv_bus: int
    r_pu: float = Field(ge=0)
    x_pu: float = Field(gt=0)
    sn_mva: float = Field(gt=0, description="Rating the r_pu/x_pu are referred to")
    tap_ratio: float = Field(
        default=1.0, gt=0,
        description="Off-nominal turns ratio on the HV side, per unit: v_hv = tap_ratio * v_lv at no "
        "load. 1.0 is the nominal ratio the two buses' vn_kv already imply. Only modelled for a branch "
        "transformer; a unit step-up's ratio is part of its unit's model.",
    )
    shift_degree: float = Field(
        default=0.0,
        description="Phase shift from HV to LV (a phase-shifting transformer). Rotates the dq frame "
        "between the two sides. Branch transformers only, as for tap_ratio.",
    )
    name: str = ""
    hv_closed: bool = Field(default=True, description="Breaker at the HV end (open = transformer out of service)")
    lv_closed: bool = Field(default=True, description="Breaker at the LV end (open = transformer out of service)")


class Load(BaseModel):
    """A constant-power load. Corresponds to the P_cons/Q_cons columns of ``Y_network``."""

    bus: int
    p_mw: float
    q_mvar: float
    name: str = ""
    closed: bool = Field(default=True, description="Breaker between the load and its bus")


class Shunt(BaseModel):
    """A shunt compensation device at a bus: a capacitor bank or a reactor.

    The sign convention is :class:`Load`'s (and pandapower's): ``q_mvar > 0``
    absorbs reactive power (a reactor), ``q_mvar < 0`` generates it (a
    capacitor bank).

    The two are modelled very differently, because a bus's own dynamic model
    (``components/node.py``) already *is* a capacitance to ground --
    ``C dv/dt = i - jwC v``. A capacitor bank is simply more of that
    capacitance, so it adds no states at all. A reactor is an inductor to
    ground and cannot be folded in: a negative capacitance would flip the
    sign of the node's dynamics, and netting the two against each other would
    only be right at exactly the nominal frequency. It therefore gets a
    branch of its own, the same RL branch a line uses with its far end at
    zero volts.
    """

    bus: int
    q_mvar: float = Field(
        description="Reactive power at 1.0 pu voltage: > 0 absorbs (reactor), < 0 generates (capacitor bank)"
    )
    r_pu: float | None = Field(
        default=None, ge=0,
        description="Series resistance of a *reactor*, per unit on the network base. None derives one from "
        "X/R = 50 (a typical shunt reactor); exactly 0 leaves its resonance with the bus capacitance "
        "undamped. Ignored for a capacitor bank, which has no branch of its own.",
    )
    closed: bool = Field(default=True, description="Breaker between the device and its bus")
    name: str = ""


class StateMode(str, Enum):
    """What becomes of one state group when the model order is lowered --
    see :func:`g2elin_core.components.base.apply_reduction`."""

    DYNAMIC = reduction.DYNAMIC
    ALGEBRAIC = reduction.ALGEBRAIC
    FROZEN = reduction.FROZEN


def _mode_values(states: dict) -> dict[str, str]:
    """``{group: mode string}`` from a ``dict[str, StateMode]``.

    Tolerates a plain string as well: ``model_copy(update=...)`` skips
    validation, so a field set that way holds whatever it was given, and a
    bare ``AttributeError`` here would say nothing about why.
    """
    return {k: (v.value if isinstance(v, StateMode) else str(v)) for k, v in states.items()}


def _check_requirements(kind: str, modes: dict[str, str]) -> None:
    """Refuse a combination whose algebraic groups depend on ones that are
    still dynamic.

    Without this the model still *builds* and then fails deep inside the
    linearisation with "singular algebraic Jacobian", which is true and
    tells the user nothing about which switch caused it.
    """
    unmet = reduction.element(kind).unmet_requirements(modes)
    if not unmet:
        return
    parts = [
        f"{group!r} needs {' and '.join(repr(m) for m in missing)} to be algebraic too"
        for group, missing in unmet
    ]
    raise ValueError(
        f"this {kind} model order isn't solvable: " + "; ".join(parts) + ". "
        "A control loop's integrator doesn't appear in its own equation, so it can only be made "
        "algebraic once what it regulates is an unknown as well."
    )


class ModelOptions(BaseModel):
    """Which dynamics this network's models keep.

    Lives on the ``Network`` rather than on each analysis request on
    purpose: the level is part of what a saved case *is*, it travels with a
    network that is exported and re-imported, and the analysis layer's model
    cache is keyed on the network's own JSON, so changing a level correctly
    invalidates every cached model without any extra plumbing.

    ``*_level`` names a preset from :mod:`g2elin_core.reduction`;
    ``*_states`` overrides individual state groups on top of it. A unit can
    override both again through ``DerUnit.level`` / ``DerUnit.states``.
    """

    network_level: str = Field(
        default="full",
        description="Passive-element dynamics: 'full' (EMT) or 'quasi_stationary' (RMS). "
        "See g2elin_core.reduction for the catalogue.",
    )
    network_states: dict[str, StateMode] = Field(
        default_factory=dict,
        description="Per-element-kind overrides (nodes/lines/loads/shunts/transformers) "
        "applied on top of network_level",
    )
    network_frequency: Literal["frame", "nominal"] = Field(
        default="frame",
        description="Whether the w*L / w*C speed terms in the passive elements follow the "
        "reference frame's own speed (the default, and what an EMT model does) or are pinned "
        "to nominal, which is the convention phasor tools use.",
    )
    sm_level: str = Field(default="full", description="Default synchronous-machine level")
    sm_states: dict[str, StateMode] = Field(default_factory=dict)
    gfm_level: str = Field(default="full", description="Default grid-forming converter level")
    gfm_states: dict[str, StateMode] = Field(default_factory=dict)
    gfl_level: str = Field(default="full", description="Default grid-following converter level")
    gfl_states: dict[str, StateMode] = Field(default_factory=dict)

    @model_validator(mode="after")
    def _levels_and_groups_exist(self) -> "ModelOptions":
        for kind in ("network", "sm", "gfm", "gfl"):
            level = getattr(self, f"{kind}_level")
            overrides = _mode_values(getattr(self, f"{kind}_states"))
            # Raises with the available levels / allowed modes named.
            modes = reduction.element(kind).modes_by_group(level, overrides)
            _check_requirements(kind, modes)
        return self

    def modes_for(self, kind: str) -> dict[str, str]:
        """``{symbol name: mode}`` for one element type."""
        overrides = _mode_values(getattr(self, f"{kind}_states"))
        return reduction.resolve_modes(kind, getattr(self, f"{kind}_level"), overrides)

    def group_modes_for(self, kind: str) -> dict[str, str]:
        """``{group id: mode}`` for one element type -- what the UI shows."""
        overrides = _mode_values(getattr(self, f"{kind}_states"))
        return reduction.element(kind).modes_by_group(getattr(self, f"{kind}_level"), overrides)

    @property
    def network_is_dynamic(self) -> bool:
        """True when any passive element still integrates something -- the
        property that decides whether a run is electromagnetic-transient."""
        return reduction.DYNAMIC in self.group_modes_for("network").values()

    @property
    def fixed_network_frequency(self) -> bool:
        return self.network_frequency == "nominal"


class DerUnit(BaseModel):
    """A dispatchable unit (infinite bus, GFM, GFL or synchronous machine).

    Corresponds to one row of ``Y_DER``. Sits on its own bus (``bus``),
    connected to the rest of the network through ``transformer`` — matching
    ``network_form.m``, which always places DER units behind a step-up
    transformer on a new node.
    """

    id: int
    bus: int
    unit_type: UnitType
    bus_type: BusType
    v_set_pu: float = Field(gt=0, description="Voltage setpoint (slack/PV) or initial guess (PQ)")
    p_set_mw: float = Field(description="Dispatched active power (generation positive)")
    q_set_mvar: float = 0.0
    p_cons_mw: float = Field(
        default=0.0, description="Small parallel load at the unit's own bus (Y_DER P_cons column)"
    )
    q_cons_mvar: float = 0.0
    controller: GfmController | None = Field(
        default=None, description="Outer power-control law, only meaningful for GFM units"
    )
    exciter: ExciterModel | None = Field(
        default=None,
        description="Which AVR/exciter model this machine carries; only meaningful for SM units. "
        "None = this tool's original regulator, which is also what every saved network predating "
        "the choice keeps.",
    )
    pss: PssModel | None = Field(
        default=None,
        description="Which power system stabiliser model this machine carries, or 'none' for a machine "
        "without one; only meaningful for SM units. None = this tool's original stabiliser.",
    )
    governor: GovernorModel | None = Field(
        default=None,
        description="Whether this machine has a governor, and which model; 'none' holds the mechanical "
        "power constant. Only meaningful for SM units. None = this tool's original governor.",
    )
    xd_pu: float | None = Field(
        default=None, description="Equivalent transient reactance, used for SCR calculations only"
    )
    sn_mva: float | None = Field(
        default=None, gt=0,
        description="The unit's own rating, i.e. the base this unit's `params` are given on. Published "
        "machine data is per unit of the machine's rating (Kundur's two-area machines are 900 MVA, while "
        "the network base is 100 MVA), so with this set those overrides are converted to the network base "
        "-- impedances x Sn_network/Sn_unit, inertia the other way. Leave it unset (the default) when the "
        "overrides are already on the network base. It rebases `params` only: the built-in defaults and "
        "the transformer impedance are on the network base either way.",
    )
    closed: bool = Field(default=True, description="Breaker between the unit and its bus")
    params: dict[str, float] = Field(
        default_factory=dict,
        description="Overrides of this unit's electrical/control parameters (e.g. KpCL, H, Ka), by the "
        "names sm_params()/gfm_params()/gfl_params() use; anything not listed keeps its default",
    )
    level: str | None = Field(
        default=None,
        description="This unit's own model level, overriding the network's default for its type "
        "(ModelOptions.sm_level / gfm_level / gfl_level). None = follow the network.",
    )
    states: dict[str, StateMode] = Field(
        default_factory=dict,
        description="Per-state-group overrides for this unit alone, applied on top of its level",
    )

    @model_validator(mode="after")
    def _controller_only_for_gfm(self) -> "DerUnit":
        if self.controller is not None and self.unit_type is not UnitType.GFM:
            raise ValueError("controller is only meaningful for GFM units")
        return self

    @model_validator(mode="after")
    def _regulators_only_for_sm(self) -> "DerUnit":
        if self.unit_type is not UnitType.SYNCHRONOUS_MACHINE:
            for name in ("exciter", "pss", "governor"):
                if getattr(self, name) is not None:
                    raise ValueError(f"{name} is only meaningful for synchronous machines")
        return self

    @property
    def exciter_model(self) -> str:
        """This machine's exciter, defaulted -- what the model builders take."""
        return (self.exciter or ExciterModel.G2ELIN).value

    @property
    def pss_model(self) -> str:
        """This machine's stabiliser, defaulted."""
        return (self.pss or PssModel.G2ELIN).value

    @property
    def governor_model(self) -> str:
        """This machine's governor, defaulted."""
        return (self.governor or GovernorModel.G2ELIN).value

    @model_validator(mode="after")
    def _level_and_groups_exist(self) -> "DerUnit":
        kind = self.unit_type.value
        if kind not in reduction.ELEMENTS:
            # An infinite bus has no reducible dynamics of its own.
            if self.level or self.states:
                raise ValueError(f"a {kind} unit has no model levels to choose from")
            return self
        if self.level is not None or self.states:
            overrides = _mode_values(self.states)
            modes = reduction.element(kind).modes_by_group(self.level, overrides)
            # Only checkable here against this unit's own settings; the
            # combination with the network's defaults is checked again in
            # Network.unit_group_modes' caller (validation.py), since a unit
            # can be valid on its own and not on top of the network's level.
            _check_requirements(kind, modes)
        return self


class Network(BaseModel):
    """A complete network case: base values, topology, loads and DER dispatch."""

    name: str
    f_hz: float = 60.0
    sn_mva: float = 100.0
    units_use_first_transformer: bool = Field(
        default=False,
        description="MATLAB-compatible mode: every unit's dynamic model uses the network's *first* "
        "transformer impedance (script_generic.m's Y_TR(1,:) convention) instead of its own. Off by "
        "default, so power flow and the dynamic models see the same transformer.",
    )
    nodes_share_first_line_b: bool = Field(
        default=False,
        description="MATLAB-compatible mode: every bus's dynamic model uses the *first* line's charging "
        "susceptance as its own capacitance, instead of the half-charging of the lines actually connected "
        "to it plus its capacitor banks. Off by default, so each bus gets its own; the ported presets set "
        "it so their numbers still match the MATLAB toolbox's (see network/breakers.node_capacitances).",
    )
    min_node_b_pu: float | None = Field(
        default=None, gt=0,
        description="The capacitance to give a bus whose in-service lines declare no charging at all, as "
        "distribution-feeder data often does (IEEE 33/69). Without it such a network is refused rather "
        "than run on an invented number, since a bus's dynamic model divides by its capacitance.",
    )
    frame_follows_slack: bool = Field(
        default=False,
        description="MATLAB-compatible mode: the dynamic models' common dq frame turns with the slack "
        "unit's own speed (script_generic.m's convention) instead of standing on its own at nominal speed. "
        "Off by default, so any unit -- the slack included -- can be disconnected and a network can be "
        "split into islands (see components/frame.py).",
    )
    models: ModelOptions = Field(
        default_factory=ModelOptions,
        description="Which dynamics the element models keep -- the model-order reduction "
        "settings. The default keeps everything, i.e. the full EMT model.",
    )
    buses: list[Bus]
    lines: list[Line] = Field(default_factory=list)
    transformers: list[Transformer] = Field(default_factory=list)
    loads: list[Load] = Field(default_factory=list)
    shunts: list[Shunt] = Field(default_factory=list)
    der_units: list[DerUnit] = Field(default_factory=list)
    # Set only on the reduced networks network.breakers.energized_network()
    # builds: the element numbering/labels of the network they came from, so
    # model names ("Ln_3", "SM_2", ...) don't shift when elements drop out.
    _labels: object = PrivateAttr(default=None)

    @model_validator(mode="after")
    def _bus_refs_exist(self) -> "Network":
        bus_ids = {b.id for b in self.buses}
        missing: set[int] = set()
        for ln in self.lines:
            missing |= {ln.from_bus, ln.to_bus} - bus_ids
        for tr in self.transformers:
            missing |= {tr.hv_bus, tr.lv_bus} - bus_ids
        for ld in self.loads:
            missing |= {ld.bus} - bus_ids
        for sh in self.shunts:
            missing |= {sh.bus} - bus_ids
        for der in self.der_units:
            missing |= {der.bus} - bus_ids
        if missing:
            raise ValueError(f"references to undefined bus ids: {sorted(missing)}")
        return self

    @model_validator(mode="after")
    def _exactly_one_slack(self) -> "Network":
        slacks = [d for d in self.der_units if d.bus_type is BusType.SLACK]
        if len(slacks) != 1:
            raise ValueError(f"expected exactly one slack DER unit, found {len(slacks)}")
        return self

    def bus(self, bus_id: int) -> Bus:
        for b in self.buses:
            if b.id == bus_id:
                return b
        raise KeyError(bus_id)

    def unit_group_modes(self, der: DerUnit) -> dict[str, str]:
        """``{group id: mode}`` for one unit: the network's default for its
        type, with the unit's own ``level``/``states`` applied on top."""
        kind = der.unit_type.value
        level = der.level or getattr(self.models, f"{kind}_level")
        overrides = _mode_values(getattr(self.models, f"{kind}_states"))
        overrides.update(_mode_values(der.states))
        modes = self.unit_element(der).modes_by_group(level, overrides)
        # A unit's own settings and the network's defaults are each valid on
        # their own; only their combination can be unsolvable, so this is
        # where that is caught.
        _check_requirements(kind, modes)
        return modes

    @staticmethod
    def unit_element(der: DerUnit) -> "reduction.ElementModel":
        """The reduction catalogue for one unit. A machine's AVR and PSS
        groups depend on which regulator models it carries, so the catalogue
        is per unit rather than per type."""
        if der.unit_type is UnitType.SYNCHRONOUS_MACHINE:
            return reduction.element("sm", exciter=der.exciter_model, pss=der.pss_model,
                                     governor=der.governor_model)
        return reduction.element(der.unit_type.value)

    def unit_modes(self, der: DerUnit) -> dict[str, str]:
        """``{symbol name: mode}`` for one unit -- what its ``*_dae()``
        builder takes."""
        e = self.unit_element(der)
        by_group = self.unit_group_modes(der)
        return {sym: by_group[g.id] for g in e.groups for sym in g.symbols}

    def unit_level(self, der: DerUnit) -> str | None:
        """The named level this unit's settings correspond to, or None when
        its per-group overrides don't match any."""
        return self.unit_element(der).matching_level(self.unit_group_modes(der))
