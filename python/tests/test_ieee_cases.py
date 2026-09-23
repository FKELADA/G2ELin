"""The IEEE test cases, and the pandapower importer that produces them.

These start life as *power-flow* cases: bus/line/trafo/gen tables with no
dynamic data at all. Turning one into a G2ELin network means supplying things
the source does not contain -- a terminal bus and step-up per generator,
machine ratings, machine dynamics, line charging -- so what these tests pin is:

* the **power flow** reproduces the published case, which is what the source
  data actually determines;
* every assumption the importer makes is *recorded* and is the one intended;
* the machines are declared on their own rating, without which a 500 MVA unit
  would carry a fifth of the inertia it should.

What they deliberately do **not** assert is mode frequencies or stability. The
machine data is a generic placeholder (``machine_data.GENERIC_MACHINE``), so
dynamic results are qualitative until real data replaces it -- on the 39-bus
case one local mode is negatively damped with it.
"""

from __future__ import annotations

import math

import pytest

pp = pytest.importorskip("pandapower")
import pandapower.networks as nw  # noqa: E402

from g2elin_core.network.machine_data import flux_linkage_params, generic_machine_params  # noqa: E402
from g2elin_core.network.pandapower_import import (  # noqa: E402
    ASSUMED_POWER_FACTOR, estimate_line_charging, from_pandapower, regulate_grid_buses,
)
from g2elin_core.network.presets import ieee14, ieee39, ieee118  # noqa: E402
from g2elin_core.network.validation import validate_network  # noqa: E402
from g2elin_core.powerflow import run_power_flow  # noqa: E402

CASES = [(ieee14, nw.case14), (ieee39, nw.case39), (ieee118, nw.case118)]
IDS = ["ieee14", "ieee39", "ieee118"]


@pytest.fixture(scope="module")
def solved_sources():
    out = {}
    for _, factory in CASES:
        net = factory()
        pp.runpp(net)
        out[factory.__name__] = net
    return out


# --- the shipped presets ---------------------------------------------------------
@pytest.mark.parametrize("build,factory", CASES, ids=IDS)
def test_the_preset_is_a_valid_network(build, factory):
    assert [i for i in validate_network(build()) if i.severity == "error"] == []


@pytest.mark.parametrize("build,factory", CASES, ids=IDS)
def test_each_generator_got_its_own_terminal_bus_and_step_up(build, factory):
    """A G2ELin unit lives behind a step-up whose LV side is its own, while the
    source case puts generators straight onto network buses that carry loads
    and lines."""
    net, source = build(), factory()
    n_units = len(source.gen) + len(source.ext_grid) + len(source.sgen)
    assert len(net.der_units) == n_units
    assert len(net.buses) == len(source.bus) + n_units
    assert len(net.transformers) == len(source.trafo) + n_units
    # Each unit's own bus carries nothing else.
    unit_buses = {d.bus for d in net.der_units}
    assert len(unit_buses) == n_units
    assert not unit_buses & {ld.bus for ld in net.loads}
    assert not unit_buses & {b for ln in net.lines for b in (ln.from_bus, ln.to_bus)}


@pytest.mark.parametrize("build,factory", CASES, ids=IDS)
def test_the_power_flow_reproduces_the_published_case(build, factory, solved_sources):
    """The one thing the source data really determines. Tolerance is loose
    enough to absorb what had to be added -- estimated line charging, the
    step-ups -- and tight enough to catch a conversion error."""
    net = build()
    source = solved_sources[factory.__name__]
    result = run_power_flow(net)
    assert result.converged

    ours = result.bus_table().set_index("bus")
    worst_v = max(abs(float(ours["vm_pu"][b]) - float(source.res_bus.vm_pu[b])) for b in source.bus.index)
    assert worst_v < 0.01, f"worst bus voltage error {worst_v:.2e} pu"

    ref_ours = float(ours["va_degree"][source.ext_grid.bus.iloc[0]])
    ref_theirs = float(source.res_bus.va_degree[source.ext_grid.bus.iloc[0]])
    worst_a = max(
        abs((float(ours["va_degree"][b]) - ref_ours) - (float(source.res_bus.va_degree[b]) - ref_theirs))
        for b in source.bus.index
    )
    assert worst_a < 1.0, f"worst bus angle error {worst_a:.2e} deg"

    published_losses = float(source.res_line.pl_mw.sum() + source.res_trafo.pl_mw.sum())
    assert result.total_losses_mw() == pytest.approx(published_losses, rel=0.05)


def test_ieee14_keeps_its_published_transformer_taps():
    """0.978, 0.969 and 0.932 -- pandapower hides them in vn_hv_kv/vn_lv_kv and
    the tap position, which the importer has to fold into one ratio."""
    taps = sorted(round(t.tap_ratio, 3) for t in ieee14().transformers if t.tap_ratio != 1.0)
    assert taps == [0.932, 0.969, 0.978]


def test_ieee14_keeps_its_capacitor_bank():
    banks = [s for s in ieee14().shunts if s.q_mvar < 0]
    assert len(banks) == 1 and banks[0].q_mvar == pytest.approx(-19.0, abs=0.5)


@pytest.mark.parametrize("build,factory", CASES, ids=IDS)
def test_machines_are_declared_on_their_own_rating(build, factory):
    """Without this a 500 MVA machine would take the generic H = 4 s on the
    100 MVA network base, i.e. a fifth of the inertia it should have."""
    net = build()
    machines = [d for d in net.der_units if d.unit_type.value == "sm"]
    assert machines
    for unit in machines:
        assert unit.sn_mva and unit.sn_mva > 0
        assert unit.params["H"] == pytest.approx(4.0)     # on the machine's own base
    from g2elin_core.operating_point import unit_params
    big = max(machines, key=lambda d: d.sn_mva)
    effective_h = unit_params(net, big)["H"]
    assert effective_h == pytest.approx(4.0 * big.sn_mva / net.sn_mva)
    assert effective_h > 4.0                               # more inertia than the network base implies


# --- the importer's own conversions ------------------------------------------------
def test_tap_ratio_folds_in_both_of_pandapowers_encodings():
    source = nw.case14()
    net, _ = from_pandapower(source, name="t")
    branch = [t for t in net.transformers if t.name.startswith("trafo")]
    assert len(branch) == len(source.trafo)
    # Transformer 0: nominal winding voltages, one tap step of -2.2%.
    assert branch[0].tap_ratio == pytest.approx(1.0 - 0.022, abs=1e-9)


def test_charging_is_estimated_only_where_the_case_gives_none():
    source = nw.case14()
    net, report = from_pandapower(source, name="t")
    zero_in_source = int((source.line.c_nf_per_km == 0).sum())
    assert len(report.estimated_line_charging) == zero_in_source
    assert all(ln.b_pu > 0 for ln in net.lines)


def test_estimated_charging_follows_the_surge_impedance_relation():
    """B = X (Z_base/Z0)^2 -- a reactance the case does give turned into a
    susceptance it does not."""
    b = estimate_line_charging(x_pu=0.1, vn_kv=135.0, sn_mva=100.0, surge_impedance_ohm=300.0)
    assert b == pytest.approx(0.1 * ((135.0**2 / 100.0) / 300.0) ** 2)


def test_zero_reactive_loads_are_given_a_little_inductive_q():
    """A constant-impedance load model has no reactance to build from when
    Q = 0, which case118 has nine of."""
    source = nw.case118()
    net, report = from_pandapower(source, name="t")
    assert len(report.loads_given_reactive_power) == int((source.load.q_mvar == 0).sum())
    assert all(ld.q_mvar != 0 for ld in net.loads)


def test_a_condenser_is_rated_from_its_reactive_capability():
    """It dispatches no active power, so a rating from P alone would be zero --
    and a near-zero rating puts a huge step-up in front of a machine whose job
    is holding a bus up."""
    source = nw.case14()
    net, report = from_pandapower(source, name="t")
    condensers = [i for i, r in source.gen.iterrows() if r.p_mw == 0]
    assert condensers
    step_ups = {t.lv_bus: t for t in net.transformers}
    for unit in net.der_units:
        assert step_ups[unit.bus].sn_mva >= 20.0          # not a 1 MVA machine
    assert report.assumed_unit_ratings                     # and it was recorded as an assumption


def test_more_than_one_external_grid_is_refused():
    source = nw.case14()
    pp.create_ext_grid(source, bus=5, vm_pu=1.0)
    with pytest.raises(ValueError, match="exactly one external grid"):
        from_pandapower(source, name="t")


def test_regulating_the_grid_bus_is_what_matches_the_published_voltages(solved_sources):
    """Without it every voltage comes out low by the step-up's drop, because
    the machine regulates its new terminal instead of the bus the case did."""
    source = solved_sources["case14"]
    raw, report = from_pandapower(source, name="t")
    before = run_power_flow(raw).bus_table().set_index("bus")["vm_pu"]
    worst_before = max(abs(float(before[b]) - float(source.res_bus.vm_pu[b])) for b in source.bus.index)

    regulate_grid_buses(raw, report, max_iter=60)
    after = run_power_flow(raw).bus_table().set_index("bus")["vm_pu"]
    worst_after = max(abs(float(after[b]) - float(source.res_bus.vm_pu[b])) for b in source.bus.index)

    assert worst_before > 0.05                             # ~0.077 pu low
    assert worst_after < 0.01
    assert worst_after < worst_before / 20


# --- the machine-data conversion ------------------------------------------------------
def test_the_generic_machine_inverts_back_to_its_published_form():
    p = generic_machine_params(f_hz=60.0)
    wb = 2 * math.pi * 60.0
    par = lambda *xs: 1.0 / sum(1.0 / x for x in xs)  # noqa: E731
    xl = p["Ll"]
    assert xl + p["Lad"] == pytest.approx(1.8)
    assert xl + par(p["Lad"], p["Lfd"]) == pytest.approx(0.30)
    assert xl + par(p["Lad"], p["Lfd"], p["L1d"]) == pytest.approx(0.22)
    assert (p["Lad"] + p["Lfd"]) / (wb * p["Rfd"]) == pytest.approx(7.0)


def test_the_conversion_is_the_one_kundur_uses_too():
    """The same function serves the published Kundur machine and the generic
    one, so there is a single conversion to be right about."""
    from g2elin_core.network.presets import kundur_machine_params

    assert kundur_machine_params(6.5) == flux_linkage_params(
        xd=1.8, xq=1.7, xl=0.2, ra=0.0025, xdp=0.3, xqp=0.55, xdpp=0.25, xqpp=0.25,
        td0p=8.0, tq0p=0.4, td0pp=0.03, tq0pp=0.05, h=6.5, kd=0.0, f_hz=60.0,
    )


# --- the diagram at this size -----------------------------------------------------
@pytest.mark.parametrize("build,factory", CASES, ids=IDS)
def test_every_bus_still_gets_a_position(build, factory):
    from g2elin_core.network.topology import _unit_terminal_buses, compute_topology_layout

    net = build()
    assert len(_unit_terminal_buses(net)) == len(net.der_units)
    assert len({n.id for n in compute_topology_layout(net).nodes}) == len(net.buses)


def test_a_large_network_hangs_its_unit_terminals_off_the_bus_they_feed():
    """Laying each one out as a node of its own pulls the real topology apart
    and spends space on it -- 54 extra nodes among the 118-bus case's own,
    whose labels overlapped 152 times before this."""
    from g2elin_core.network.topology import _unit_terminal_buses, compute_topology_layout

    from g2elin_core.network.topology import _TERMINAL_OFFSET

    net = ieee118()
    terminals = _unit_terminal_buses(net)
    at = {n.id: (n.x, n.y) for n in compute_topology_layout(net).nodes}
    for terminal, grid_bus in terminals.items():
        distance = math.hypot(at[terminal][0] - at[grid_bus][0], at[terminal][1] - at[grid_bus][1])
        assert distance == pytest.approx(_TERMINAL_OFFSET, abs=1e-9)


@pytest.mark.parametrize("build", [ieee39, lambda: __import__(
    "g2elin_core.network.presets", fromlist=["x"]).wscc9_3sm()], ids=["ieee39", "wscc9_3sm"])
def test_a_small_network_keeps_the_spring_layout_for_its_terminals(build):
    """A fixed offset pushed the machines onto the very branches their bus
    connects to. Below the threshold the spring layout has room and does it
    better, so it keeps them as nodes."""
    from g2elin_core.network.topology import _unit_terminal_buses, compute_topology_layout

    from g2elin_core.network.topology import _TERMINAL_OFFSET

    net = build()
    terminals = _unit_terminal_buses(net)
    at = {n.id: (n.x, n.y) for n in compute_topology_layout(net).nodes}
    distances = [
        math.hypot(at[t][0] - at[g][0], at[t][1] - at[g][1]) for t, g in terminals.items()
    ]
    # Placed by the layout, so spaced by the topology -- not all sitting at the
    # one fixed offset a hung terminal would take.
    assert all(d != pytest.approx(_TERMINAL_OFFSET, abs=1e-6) for d in distances)
    assert all(d > 0 for d in distances)


def test_two_units_on_one_bus_do_not_land_on_top_of_each_other():
    """Only reachable on a network big enough to hang them."""
    from g2elin_core.network.schema import Bus, BusType, DerUnit, Transformer, UnitType
    from g2elin_core.network.topology import _unit_terminal_buses, compute_topology_layout

    net = ieee118()
    grid = next(t.hv_bus for t in net.transformers if t.lv_bus == net.der_units[1].bus)
    new_bus = max(b.id for b in net.buses) + 1
    net.buses.append(Bus(id=new_bus, name="extra_terminal", vn_kv=20.0))
    net.transformers.append(Transformer(hv_bus=grid, lv_bus=new_bus, r_pu=0.0, x_pu=0.05, sn_mva=100.0))
    net.der_units.append(DerUnit(id=9999, bus=new_bus, unit_type=UnitType.GFL, bus_type=BusType.PQ,
                                 v_set_pu=1.0, p_set_mw=1.0))
    sharing = [t for t, g in _unit_terminal_buses(net).items() if g == grid]
    assert len(sharing) == 2
    at = {n.id: (n.x, n.y) for n in compute_topology_layout(net).nodes}
    a, b = (at[t] for t in sharing)
    assert math.hypot(a[0] - b[0], a[1] - b[1]) > 1e-6
