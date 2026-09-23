"""Bus shunts, and each bus's own capacitance.

Both are the same piece of physics. A bus's dynamic model (components/node.py)
*is* a capacitance to ground -- ``C dv/dt = i - jwC v`` -- so a capacitor bank
is simply more of that capacitance and adds no states, while the charging of
the lines meeting at a bus is where that capacitance comes from in the first
place. A reactor cannot be folded in (a negative capacitance would flip the
sign of the node's dynamics) and gets an RL branch to ground of its own.

Before this, every bus borrowed the *first line's* charging, whatever was
actually connected to it -- the MATLAB toolbox's convention, which the ported
presets keep (``Network.nodes_share_first_line_b``) so their numbers still
match it.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from g2elin_core.network.breakers import energized_network, node_capacitances, shunt_reactor_pq_mw
from g2elin_core.network.presets import cigre_interconnected_1sm_1gfm_1gfl, wscc9_3sm
from g2elin_core.network.schema import Line, Shunt
from g2elin_core.pipeline import linearize_network
from g2elin_core.powerflow import run_power_flow
from g2elin_core.timedomain import build_nonlinear_network


def _fixed(build=wscc9_3sm):
    """A preset with each bus using its own capacitance -- the default since
    the presets stopped setting the MATLAB convention, so this is now just
    ``build()``; kept as a name because the tests below read better for it."""
    net = build()
    net.nodes_share_first_line_b = False
    return net


def _matlab(build=wscc9_3sm):
    """The same preset on the MATLAB toolbox's convention: every bus borrows
    the first line's charging. No longer the default, still supported, and
    still what the ported cases have to be able to reproduce."""
    net = build()
    net.nodes_share_first_line_b = True
    return net


def _node_residual(net) -> float:
    """max |dv/dt| over the bus voltages at the power flow's operating point.

    Zero means the point the model is linearized about really is an
    equilibrium of that model.
    """
    model = build_nonlinear_network(net, run_power_flow(net))
    xdot, _, _ = model.rhs(model.initial_state(), model.default_u_exo())
    nodes = [i for i, n in enumerate(model.state_names) if n.startswith("v_{g_")]
    return float(np.abs(xdot[nodes]).max())


# --- each bus's own capacitance ----------------------------------------------
def test_each_bus_gets_half_the_charging_of_its_own_lines():
    net = _fixed()
    cap = node_capacitances(net)
    for bus_id, value in cap.items():
        expected = sum(ln.b_pu / 2 for ln in net.lines if bus_id in (ln.from_bus, ln.to_bus))
        assert value == pytest.approx(expected)
    # WSCC's buses genuinely differ -- this is not one number in disguise.
    assert len(set(round(v, 6) for v in cap.values())) > 1


def test_parity_mode_gives_every_bus_the_first_lines_charging():
    net = _matlab()
    assert set(node_capacitances(net).values()) == {net.lines[0].b_pu}


def test_a_units_own_terminal_bus_has_no_node_of_its_own():
    """It lives inside the unit's model, behind the step-up transformer, so it
    has no node equation and needs no capacitance."""
    net = _fixed()
    assert {d.bus for d in net.der_units}.isdisjoint(node_capacitances(net))


def test_per_bus_capacitance_makes_the_operating_point_an_equilibrium():
    """The borrowed value leaves the model being linearized about a point that
    is not its own equilibrium: bus voltages drift with no disturbance.

    This is why the presets no longer borrow it. With a dynamic network the
    drift is a transient that decays in microseconds; with a quasi-stationary
    one it cannot decay at all, and moves the answer instead (see
    ``network/validation._model_order_issues``)."""
    borrowed = _node_residual(_matlab())
    fixed = _node_residual(wscc9_3sm())
    assert borrowed > 100.0                      # ~178 pu/s on this preset
    assert fixed < 1.0                           # ~0.09
    assert fixed < borrowed / 100


def test_what_is_left_over_is_only_the_machines_own_auxiliary_load():
    """On CIGRE the remainder is exactly the small load the SM model carries
    internally (SM.PL = Sb/RL_pu), which the power flow doesn't represent --
    inherited from the MATLAB toolbox, not from the capacitance."""
    net = _fixed(cigre_interconnected_1sm_1gfm_1gfl)
    with_aux = _node_residual(net)
    net.der_units[0].params = {"RL_pu": 1e9}     # effectively no auxiliary load
    assert _node_residual(net) < with_aux / 1000


def test_a_bus_with_no_charging_is_refused_rather_than_divided_by():
    """Distribution-feeder data (IEEE 33/69) often declares no line charging
    at all, which would divide by zero in every bus's equation."""
    net = _fixed()
    for ln in net.lines:
        ln.b_pu = 0.0
    with pytest.raises(ValueError, match="no shunt capacitance"):
        node_capacitances(net)


def test_min_node_b_pu_is_the_way_to_say_what_to_use_instead():
    net = _fixed()
    for ln in net.lines:
        ln.b_pu = 0.0
    net.min_node_b_pu = 1e-4
    assert set(node_capacitances(net).values()) == {1e-4}
    assert linearize_network(net, run_power_flow(net)).A.shape[0] == 88


def test_switching_a_line_out_takes_its_charging_with_it():
    """Via energized_network, which is what drops out-of-service elements --
    node_capacitances itself counts whatever the network still lists, so it
    always agrees with the blocks that were built."""
    net = _fixed()
    before = node_capacitances(net)
    ln = net.lines[0]
    net.lines[0].from_closed = False
    after = node_capacitances(energized_network(net))
    assert after[ln.from_bus] == pytest.approx(before[ln.from_bus] - ln.b_pu / 2)
    assert after[ln.to_bus] == pytest.approx(before[ln.to_bus] - ln.b_pu / 2)


# --- capacitor banks ----------------------------------------------------------
def test_a_capacitor_bank_is_capacitance_not_a_new_block():
    net = _fixed()
    plain = linearize_network(net, run_power_flow(net))
    banked = _fixed()
    banked.shunts.append(Shunt(bus=3, q_mvar=-19.0, name="bank"))
    with_bank = linearize_network(banked, run_power_flow(banked))

    assert with_bank.A.shape[0] == plain.A.shape[0]                 # no states added
    assert node_capacitances(banked)[3] == pytest.approx(node_capacitances(net)[3] + 0.19)
    assert not np.allclose(with_bank.A, plain.A)                    # but the model did change


def test_a_capacitor_bank_raises_its_buss_voltage():
    net = _fixed()
    before = run_power_flow(net).bus_table().set_index("bus").loc[3, "vm_pu"]
    net.shunts.append(Shunt(bus=3, q_mvar=-19.0))
    after = run_power_flow(net).bus_table().set_index("bus").loc[3, "vm_pu"]
    assert after > before


# --- reactors ------------------------------------------------------------------
def test_a_reactor_is_a_branch_of_its_own_and_lowers_the_voltage():
    net = _fixed()
    before = run_power_flow(net).bus_table().set_index("bus").loc[3, "vm_pu"]
    n_states = linearize_network(net, run_power_flow(net)).A.shape[0]

    net.shunts.append(Shunt(bus=3, q_mvar=25.0, name="reactor"))
    result = run_power_flow(net)
    sys = linearize_network(net, result)
    assert result.bus_table().set_index("bus").loc[3, "vm_pu"] < before
    assert sys.A.shape[0] == n_states + 2                            # its own d/q current
    assert [n for n in sys.state_names if "Sh_1" in n] == ["i_{l_d}_{Sh_1}", "i_{l_q}_{Sh_1}"]


def test_the_power_flow_and_the_dynamic_model_agree_on_what_a_reactor_draws():
    """The two used to disagree by the reactor's resistance: pandapower saw a
    pure susceptance while the branch had an X/R."""
    net = _fixed()
    shunt = Shunt(bus=3, q_mvar=25.0)
    net.shunts.append(shunt)
    result = run_power_flow(net)
    model = build_nonlinear_network(net, result)

    (ild, ilq), (vgd, vgq) = model.op.shunt_i0[0], model.op.node_vg[3]
    p_dyn = (vgd * ild + vgq * ilq) * net.sn_mva
    q_dyn = (vgq * ild - vgd * ilq) * net.sn_mva
    assert p_dyn == pytest.approx(float(result.net.res_shunt.p_mw.iloc[0]), rel=1e-3)
    assert q_dyn == pytest.approx(float(result.net.res_shunt.q_mvar.iloc[0]), rel=1e-3)
    # And both come from the one nameplate-to-branch conversion.
    p_nom, q_nom = shunt_reactor_pq_mw(shunt, net.sn_mva)
    assert (p_nom, q_nom) == pytest.approx((p_dyn, q_dyn), rel=1e-2)


def test_a_reactor_keeps_the_operating_point_an_equilibrium():
    net = _fixed()
    net.shunts.append(Shunt(bus=3, q_mvar=25.0))
    assert _node_residual(net) < 1.0


def test_a_reactors_reactance_comes_from_its_nameplate():
    """Q = V^2/X at nominal voltage, and X/R = 50 unless told otherwise."""
    net = _fixed()
    net.shunts.append(Shunt(bus=3, q_mvar=25.0))
    r_pu, x_pu = build_nonlinear_network(net, run_power_flow(net)).op.shunt_rx[0]
    assert x_pu == pytest.approx(net.sn_mva / 25.0)
    assert r_pu == pytest.approx(x_pu / 50.0)


def test_an_open_breaker_removes_a_shunt_from_both_sides():
    net = _fixed()
    banked = node_capacitances(net)[3]
    net.shunts.append(Shunt(bus=3, q_mvar=-19.0, closed=False))
    assert node_capacitances(energized_network(net))[3] == pytest.approx(banked)
    assert not bool(run_power_flow(net).net.shunt.in_service.iloc[0])


def test_an_open_reactor_breaker_drops_its_block():
    net = _fixed()
    n_states = linearize_network(net, run_power_flow(net)).A.shape[0]
    net.shunts.append(Shunt(bus=3, q_mvar=25.0, closed=False))
    reduced = energized_network(net)
    assert linearize_network(reduced, run_power_flow(reduced)).A.shape[0] == n_states


# --- the presets use their own capacitance; the MATLAB one is still there -----
@pytest.mark.parametrize("build", [wscc9_3sm, cigre_interconnected_1sm_1gfm_1gfl])
def test_presets_give_every_bus_its_own_capacitance(build):
    """The ported presets used to carry the toolbox's convention so their
    numbers matched it exactly. They no longer do: it is a convention rather
    than physics, and it is incompatible with a quasi-stationary network."""
    net = build()
    assert net.nodes_share_first_line_b is False
    caps = node_capacitances(net)
    assert len(set(caps.values())) > 1, "every bus borrowing one value is the old convention"


@pytest.mark.parametrize("build", [wscc9_3sm, cigre_interconnected_1sm_1gfm_1gfl])
def test_the_matlab_convention_is_still_available(build):
    """Dropping it as the default must not drop it as an option -- it is how
    a result gets compared against the original toolbox."""
    net = _matlab(build)
    assert set(node_capacitances(net).values()) == {net.lines[0].b_pu}
    linearize_network(net, run_power_flow(net))
