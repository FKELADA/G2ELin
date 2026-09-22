"""Builds a :class:`~g2elin_core.network.schema.Network` from a pandapower net.

This is how the IEEE test cases get in: pandapower ships them
(``pandapower.networks.case14`` and friends) as bus/line/trafo/gen tables,
which is a *power-flow* description. A G2ELin network is a dynamic one, and
the two differ in three ways that this module has to bridge.

**Every unit needs its own terminal bus.** A G2ELin unit sits behind a step-up
transformer whose LV side is the unit's alone, because that transformer's
impedance goes inside the unit's model (see ``schema.Transformer``). A
power-flow case instead puts its generators straight onto network buses, which
also carry loads and lines. So each generator gets a *new* bus and a *new*
step-up, leaving the original bus an ordinary node. The step-up is not in the
source data, so its reactance is a stated assumption
(``UNIT_STEP_UP_X_PU``, 0.15 pu on the machine's own rating -- a typical value
for a generator step-up).

**Machine ratings are usually missing.** The cases leave ``gen.sn_mva`` as NaN,
so a rating is derived from the dispatch at an assumed power factor. Anything
the case does give is used instead.

**Line charging is often left out.** Published line data frequently carries no
susceptance at all -- every line of ``case33bw``, nine of ``case14``'s fifteen
-- and a bus with no capacitance has no dynamic model, since its equation
divides by it (see ``breakers.node_capacitances``). ``estimate_line_charging``
fills those in from each line's own reactance and an assumed surge impedance,
which is a physical estimate rather than an invented constant.

Everything assumed rather than read is collected in the returned
:class:`ImportReport`, so a caller can see exactly what was not in the data.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from .machine_data import generic_machine_params
from .schema import Bus, BusType, DerUnit, Line, Load, Network, Shunt, Transformer, UnitType

# A generator step-up's reactance, per unit of the machine's own rating. Not in
# any power-flow case; 0.10-0.20 pu is the usual range for one.
UNIT_STEP_UP_X_PU = 0.15
# Rating assumed for a machine whose case gives none: its dispatch at this
# power factor. A synchronous condenser dispatches no active power at all, so
# its reactive capability is used instead, and failing that a share of the
# system load -- a rating near zero would put a huge step-up reactance in the
# way of a machine that has to hold a bus up.
ASSUMED_POWER_FACTOR = 0.85
MIN_RATING_SHARE_OF_LOAD = 0.05
# A load's dynamic model is a constant-impedance R-X equivalent, so a purely
# active load (Q = 0, which case118 has nine of) has no reactance to build it
# from. Such a load is given this much inductive Q as a fraction of its P --
# a power factor of 0.99995, far inside the rounding of published load data.
ZERO_Q_LOAD_TAN_PHI = 0.01
# Surge impedance used to estimate a line's charging when it declares none.
# Overhead transmission lines sit near 300 ohm.
DEFAULT_SURGE_IMPEDANCE_OHM = 300.0
# Bounds for regulate_grid_buses' search. A generator's terminal has to be able
# to sit well above its grid bus (it is behind a step-up), but not wander.
SECANT_MIN_SLOPE, SECANT_MAX_SLOPE = 0.15, 3.0
MAX_SETPOINT_STEP = 0.05
MIN_SETPOINT, MAX_SETPOINT = 0.80, 1.40
# A bus joined only by transformers has no line charging to take capacitance
# from, and a transformer is a series element that supplies none -- yet every
# bus needs some for its dynamic model to exist. What is really there is the
# substation: busbars, bushings, surge capacitors. This is a representative
# value for that, and it is what Network.min_node_b_pu is set from.
BUSBAR_CAPACITANCE_UF = 0.05
# Nominal voltage given to a generator's new terminal bus. It has no effect on
# any per-unit result -- the step-up's impedance is per unit of its own rating,
# and its ratio is nominal -- it only makes the diagram read sensibly.
TERMINAL_KV = 20.0


@dataclass
class ImportReport:
    """What had to be assumed, because the source case did not say."""

    name: str
    # unit id -> (the grid bus it was attached to, the voltage the case holds there)
    unit_grid_bus: dict[int, tuple[int, float]] = field(default_factory=dict)
    assumed_unit_ratings: dict[int, float] = field(default_factory=dict)
    estimated_line_charging: list[int] = field(default_factory=list)
    loads_given_reactive_power: list[int] = field(default_factory=list)
    step_up_x_pu: float = UNIT_STEP_UP_X_PU
    notes: list[str] = field(default_factory=list)

    def summary(self) -> str:
        bits = [f"{self.name}:"]
        if self.assumed_unit_ratings:
            bits.append(f"{len(self.assumed_unit_ratings)} unit rating(s) assumed at pf={ASSUMED_POWER_FACTOR}")
        if self.estimated_line_charging:
            bits.append(f"charging estimated for {len(self.estimated_line_charging)} line(s)")
        if self.loads_given_reactive_power:
            bits.append(f"{len(self.loads_given_reactive_power)} load(s) given a small inductive Q "
                        f"(tan(phi)={ZERO_Q_LOAD_TAN_PHI})")
        bits.append(f"step-up x={self.step_up_x_pu} pu on the machine base")
        bits.append("machines are the GENERIC one (no dynamic data in a power-flow case)")
        bits.extend(self.notes)
        return "; ".join(bits)


def _z_base(vn_kv: float, sn_mva: float) -> float:
    return vn_kv * vn_kv / sn_mva


def estimate_line_charging(x_pu: float, vn_kv: float, sn_mva: float,
                           surge_impedance_ohm: float = DEFAULT_SURGE_IMPEDANCE_OHM) -> float:
    """A line's total charging susceptance in pu, from its reactance.

    A line's surge impedance is ``Z0 = sqrt(L/C)``, so with
    ``X_pu = wL/Z_base`` and ``B_pu = wC*Z_base``::

        L/C = X_pu * Z_base^2 / B_pu = Z0^2   ->   B_pu = X_pu * (Z_base/Z0)^2

    which turns a reactance the case does give into the susceptance it does
    not, at the one stated assumption of ``Z0``.
    """
    return x_pu * (_z_base(vn_kv, sn_mva) / surge_impedance_ohm) ** 2


def _tap_ratio(net, row) -> float:
    """The off-nominal turns ratio of a pandapower transformer, HV side.

    pandapower encodes it twice: in ``vn_hv_kv``/``vn_lv_kv`` against the two
    buses' own nominal voltages, and in the tap changer's position. Both are
    folded into the single ratio ``schema.Transformer.tap_ratio`` carries.
    """
    hv_nom = float(net.bus.at[row.hv_bus, "vn_kv"])
    lv_nom = float(net.bus.at[row.lv_bus, "vn_kv"])
    ratio = (float(row.vn_hv_kv) / hv_nom) / (float(row.vn_lv_kv) / lv_nom)

    pos, neutral = row.get("tap_pos"), row.get("tap_neutral")
    step = row.get("tap_step_percent")
    if pos is not None and step is not None and not (_isnan(pos) or _isnan(step)):
        neutral = 0.0 if neutral is None or _isnan(neutral) else float(neutral)
        factor = 1.0 + (float(pos) - neutral) * float(step) / 100.0
        side = str(row.get("tap_side") or "hv").lower()
        ratio = ratio * factor if side == "hv" else ratio / factor
    return ratio


def _isnan(v) -> bool:
    try:
        return math.isnan(float(v))
    except (TypeError, ValueError):
        return False


def _finite(v, default=None):
    if v is None or _isnan(v):
        return default
    return float(v)


def from_pandapower(
    pp_net,
    *,
    name: str,
    slack_unit_type: UnitType = UnitType.SYNCHRONOUS_MACHINE,
    sgen_unit_type: UnitType = UnitType.GFL,
    step_up_x_pu: float = UNIT_STEP_UP_X_PU,
    surge_impedance_ohm: float = DEFAULT_SURGE_IMPEDANCE_OHM,
    f_hz: float = 60.0,
) -> tuple[Network, ImportReport]:
    """Convert a solved-or-unsolved pandapower net into a G2ELin network.

    Returns the network and an :class:`ImportReport` naming everything that was
    assumed rather than read.
    """
    sn_mva = float(pp_net.sn_mva) or 100.0
    report = ImportReport(name=name, step_up_x_pu=step_up_x_pu)
    omega = 2 * math.pi * f_hz

    buses = [
        Bus(id=int(i), name=str(row.get("name") or f"bus{int(i)}"), vn_kv=float(row.vn_kv))
        for i, row in pp_net.bus.iterrows()
    ]
    vn_of = {b.id: b.vn_kv for b in buses}
    next_bus_id = max(vn_of) + 1 if vn_of else 1

    lines: list[Line] = []
    for idx, row in pp_net.line.iterrows():
        if not bool(row.get("in_service", True)):
            continue
        km = float(row.length_km) or 1.0
        parallel = max(1, int(_finite(row.get("parallel"), 1) or 1))
        vn = vn_of[int(row.from_bus)]
        zb = _z_base(vn, sn_mva)
        r_pu = float(row.r_ohm_per_km) * km / zb / parallel
        x_pu = float(row.x_ohm_per_km) * km / zb / parallel
        c_nf = _finite(row.get("c_nf_per_km"), 0.0) or 0.0
        b_pu = omega * c_nf * 1e-9 * km * parallel * zb
        if b_pu <= 0.0:
            b_pu = estimate_line_charging(x_pu, vn, sn_mva, surge_impedance_ohm)
            report.estimated_line_charging.append(len(lines))
        lines.append(Line(
            from_bus=int(row.from_bus), to_bus=int(row.to_bus),
            r_pu=r_pu, x_pu=max(x_pu, 1e-9), b_pu=b_pu, length_km=km,
            name=str(row.get("name") or f"line {int(row.from_bus)}-{int(row.to_bus)}"),
        ))

    transformers: list[Transformer] = []
    for idx, row in pp_net.trafo.iterrows():
        if not bool(row.get("in_service", True)):
            continue
        rating = float(row.sn_mva)
        z_pct = float(row.vk_percent)
        r_pct = _finite(row.get("vkr_percent"), 0.0) or 0.0
        z_pu = z_pct / 100.0
        r_pu = r_pct / 100.0
        x_pu = math.sqrt(max(z_pu * z_pu - r_pu * r_pu, 0.0)) or z_pu
        transformers.append(Transformer(
            hv_bus=int(row.hv_bus), lv_bus=int(row.lv_bus),
            r_pu=r_pu, x_pu=max(x_pu, 1e-9), sn_mva=rating,
            tap_ratio=_tap_ratio(pp_net, row),
            shift_degree=_finite(row.get("shift_degree"), 0.0) or 0.0,
            name=str(row.get("name") or f"trafo {int(row.hv_bus)}-{int(row.lv_bus)}"),
        ))

    loads: list[Load] = []
    for _, row in pp_net.load.iterrows():
        if not bool(row.get("in_service", True)):
            continue
        p_mw = float(row.p_mw)
        q_mvar = float(_finite(row.get("q_mvar"), 0.0) or 0.0)
        if not p_mw and not q_mvar:
            continue
        if q_mvar == 0.0:
            q_mvar = abs(p_mw) * ZERO_Q_LOAD_TAN_PHI
            report.loads_given_reactive_power.append(len(loads))
        loads.append(Load(bus=int(row.bus), p_mw=p_mw, q_mvar=q_mvar,
                          name=str(row.get("name") or f"load at {int(row.bus)}")))

    shunts = [
        Shunt(bus=int(row.bus), q_mvar=float(row.q_mvar),
              name=str(row.get("name") or f"shunt at {int(row.bus)}"))
        for _, row in pp_net.shunt.iterrows()
        if bool(row.get("in_service", True)) and float(row.q_mvar)
    ]

    # --- units -----------------------------------------------------------------
    # Each one gets a terminal bus and a step-up of its own, so the bus it was
    # attached to stays an ordinary network node with its loads and lines.
    der_units: list[DerUnit] = []
    unit_id = 0

    total_load_mw = sum(abs(ld.p_mw) for ld in loads) or sn_mva

    def add_unit(grid_bus: int, p_mw: float, v_set: float, unit_type: UnitType,
                 bus_type: BusType, rating: float | None, label: str,
                 q_capability_mvar: float = 0.0) -> None:
        nonlocal next_bus_id, unit_id
        unit_id += 1
        terminal = next_bus_id
        next_bus_id += 1
        buses.append(Bus(id=terminal, name=f"{label}_terminal", vn_kv=TERMINAL_KV))
        if rating is None:
            rating = max(
                abs(p_mw) / ASSUMED_POWER_FACTOR,
                abs(q_capability_mvar),
                MIN_RATING_SHARE_OF_LOAD * total_load_mw,
            )
            report.assumed_unit_ratings[unit_id] = rating
        transformers.append(Transformer(
            hv_bus=grid_bus, lv_bus=terminal, r_pu=0.0, x_pu=step_up_x_pu,
            sn_mva=rating, name=f"{label}_xfmr",
        ))
        # A power-flow case carries no dynamic data, so the machine is the
        # generic one -- declared on its *own* rating, which is what makes it
        # physically sensible: left on the network base, a 500 MVA machine
        # would end up with a fifth of the inertia it should have.
        params = generic_machine_params(f_hz) if unit_type is UnitType.SYNCHRONOUS_MACHINE else {}
        der_units.append(DerUnit(
            id=unit_id, bus=terminal, unit_type=unit_type, bus_type=bus_type,
            v_set_pu=v_set, p_set_mw=p_mw, xd_pu=0.3,
            sn_mva=rating, params=params,
        ))
        report.unit_grid_bus[unit_id] = (grid_bus, v_set)

    if len(pp_net.ext_grid) != 1:
        raise ValueError(
            f"expected exactly one external grid to become the slack, found {len(pp_net.ext_grid)}"
        )
    eg = pp_net.ext_grid.iloc[0]
    # The slack's dispatch is whatever the power flow decides, so its rating
    # comes from what it would have to carry: the case's own solved value when
    # there is one, else the load the other units do not cover.
    slack_mw = abs(float(pp_net.res_ext_grid.p_mw.iloc[0])) if len(getattr(pp_net, "res_ext_grid", [])) else 0.0
    if not slack_mw:
        slack_mw = max(total_load_mw - sum(abs(float(r.p_mw)) for _, r in pp_net.gen.iterrows()), 0.0)
    add_unit(int(eg.bus), 0.0, float(_finite(eg.get("vm_pu"), 1.0) or 1.0),
             slack_unit_type, BusType.SLACK, None, "slack", q_capability_mvar=slack_mw / ASSUMED_POWER_FACTOR)

    for _, row in pp_net.gen.iterrows():
        if not bool(row.get("in_service", True)):
            continue
        q_cap = max(abs(_finite(row.get("max_q_mvar"), 0.0) or 0.0),
                    abs(_finite(row.get("min_q_mvar"), 0.0) or 0.0))
        add_unit(int(row.bus), float(row.p_mw), float(_finite(row.get("vm_pu"), 1.0) or 1.0),
                 UnitType.SYNCHRONOUS_MACHINE, BusType.PV, _finite(row.get("sn_mva")),
                 f"gen{int(row.bus)}", q_capability_mvar=q_cap)

    for _, row in pp_net.sgen.iterrows():
        if not bool(row.get("in_service", True)):
            continue
        add_unit(int(row.bus), float(row.p_mw), 1.0, sgen_unit_type, BusType.PQ,
                 _finite(row.get("sn_mva")), f"sgen{int(row.bus)}")

    # The floor for a bus that no line reaches (case14 has two, joined only by
    # transformers): the substation's own capacitance, at the network's highest
    # voltage, which is where such buses sit.
    top_kv = max(vn_of.values()) if vn_of else 1.0
    min_node_b_pu = omega * BUSBAR_CAPACITANCE_UF * 1e-6 * _z_base(top_kv, sn_mva)
    report.notes.append(
        f"buses reached by no line use min_node_b_pu={min_node_b_pu:.2e} "
        f"({BUSBAR_CAPACITANCE_UF} uF of substation capacitance at {top_kv:g} kV)"
    )

    network = Network(
        name=name, f_hz=f_hz, sn_mva=sn_mva,
        min_node_b_pu=min_node_b_pu,
        buses=buses, lines=lines, transformers=transformers,
        loads=loads, shunts=shunts, der_units=der_units,
    )
    return network, report


def regulate_grid_buses(network: Network, report: ImportReport, *, tol: float = 1e-5,
                        max_iter: int = 20) -> Network:
    """Retune each unit's voltage setpoint so the *grid* bus it feeds sits where
    the source case put it.

    A published case regulates the bus the generator is on. Giving each unit a
    terminal bus of its own moves the regulated point behind the step-up, so
    every voltage in the network comes out low by that transformer's drop --
    0.05 to 0.08 pu on ``case14``.

    Correcting each setpoint by its own error alone converges only slowly: at a
    strongly connected bus, raising a machine's terminal mostly pushes reactive
    power out into the network instead of lifting the bus, so the sensitivity
    is far below one (nearer 0.1 on ``case14``). Each unit therefore keeps its
    previous (setpoint, voltage) pair and takes a secant step, which needs no
    extra solves and converges in a few.

    Mutates and returns ``network``.
    """
    from g2elin_core.powerflow import run_power_flow  # lazy: powerflow imports this package

    by_id = {d.id: d for d in network.der_units}
    previous: dict[int, tuple[float, float]] = {}
    best = {uid: by_id[uid].v_set_pu for uid in report.unit_grid_bus}
    best_worst = float("inf")
    worst = float("inf")
    for _ in range(max_iter):
        result = run_power_flow(network)
        table = result.bus_table().set_index("bus")["vm_pu"]
        if not result.converged or not all(math.isfinite(float(table[g])) for g, _ in report.unit_grid_bus.values()):
            break                                  # keep the best setpoints found so far
        worst = max(abs(t - float(table[g])) for g, t in report.unit_grid_bus.values())
        if worst < best_worst:
            best_worst = worst
            best = {uid: by_id[uid].v_set_pu for uid in report.unit_grid_bus}
        if worst < tol:
            break
        for unit_id, (grid_bus, target) in report.unit_grid_bus.items():
            unit = by_id[unit_id]
            actual = float(table[grid_bus])
            error = target - actual
            # Secant on this unit's own setpoint. The slope is bounded well away
            # from zero: an unbounded one turns a small error into a huge step,
            # which walks the whole network into non-convergence.
            slope = 1.0
            if unit_id in previous:
                prev_set, prev_actual = previous[unit_id]
                d_set = unit.v_set_pu - prev_set
                if abs(d_set) > 1e-9:
                    slope = min(max((actual - prev_actual) / d_set, SECANT_MIN_SLOPE), SECANT_MAX_SLOPE)
            previous[unit_id] = (unit.v_set_pu, actual)
            step = max(min(error / slope, MAX_SETPOINT_STEP), -MAX_SETPOINT_STEP)
            unit.v_set_pu = min(max(unit.v_set_pu + step, MIN_SETPOINT), MAX_SETPOINT)
    for unit_id, v in best.items():
        by_id[unit_id].v_set_pu = v
    if best_worst >= tol:
        report.notes.append(f"grid-bus voltages matched to {best_worst:.2e} pu, not {tol:g}")
    return network
