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

from pydantic import BaseModel, Field, PrivateAttr, model_validator


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
    """

    hv_bus: int
    lv_bus: int
    r_pu: float = Field(ge=0)
    x_pu: float = Field(gt=0)
    sn_mva: float = Field(gt=0, description="Rating the r_pu/x_pu are referred to")
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
    xd_pu: float | None = Field(
        default=None, description="Equivalent transient reactance, used for SCR calculations only"
    )
    closed: bool = Field(default=True, description="Breaker between the unit and its bus")
    params: dict[str, float] = Field(
        default_factory=dict,
        description="Overrides of this unit's electrical/control parameters (e.g. KpCL, H, Ka), by the "
        "names sm_params()/gfm_params()/gfl_params() use; anything not listed keeps its default",
    )

    @model_validator(mode="after")
    def _controller_only_for_gfm(self) -> "DerUnit":
        if self.controller is not None and self.unit_type is not UnitType.GFM:
            raise ValueError("controller is only meaningful for GFM units")
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
    buses: list[Bus]
    lines: list[Line] = Field(default_factory=list)
    transformers: list[Transformer] = Field(default_factory=list)
    loads: list[Load] = Field(default_factory=list)
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
