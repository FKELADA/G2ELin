"""Transformers between two grid buses.

One schema element covers two quite different things, told apart by whether a
unit sits on the LV bus:

* a **unit step-up**, whose impedance is part of that unit's own model (the
  machine-behind-transformer-impedance formulation) and whose LV bus never
  appears in the dynamic model;
* a **branch transformer** between two grid buses, which is a branch of the
  network like a line.

Only the first used to exist. A branch transformer carried power in the power
flow and was absent from the dynamic model entirely, so the point the model was
linearized about was not an equilibrium at the two buses it joined -- silently,
with no validation error. These tests pin the fix, including the turns ratio
and the phase shift, which the wiring applies as plain coefficients.
"""

from __future__ import annotations

import numpy as np
import pytest

from g2elin_core.network.breakers import branch_transformer_indices, unit_transformers
from g2elin_core.network.presets import wscc9_3sm
from g2elin_core.network.schema import Transformer
from g2elin_core.network.validation import validate_network
from g2elin_core.pipeline import linearize_network
from g2elin_core.powerflow import run_power_flow
from g2elin_core.timedomain import build_nonlinear_network

BRANCH = dict(hv_bus=4, lv_bus=5, r_pu=0.01, x_pu=0.085, sn_mva=100.0)


def _net(**transformer):
    net = wscc9_3sm()
    net.nodes_share_first_line_b = False       # so the only residual left is the machines' own
    if transformer:
        net.transformers.append(Transformer(name="branch", **{**BRANCH, **transformer}))
    return net


def _node_residuals(net) -> dict[int, float]:
    """max |dv/dt| per bus at the power flow's operating point."""
    model = build_nonlinear_network(net, run_power_flow(net))
    xdot, _, _ = model.rhs(model.initial_state(), model.default_u_exo())
    out: dict[int, float] = {}
    for i, name in enumerate(model.state_names):
        if name.startswith("v_{g_"):
            bus = int(name.split("{Nd_")[1].rstrip("}"))
            out[bus] = max(out.get(bus, 0.0), abs(float(xdot[i])))
    return out


# --- telling the two kinds apart ----------------------------------------------
def test_a_preset_has_only_unit_step_ups():
    net = wscc9_3sm()
    assert sorted(unit_transformers(net)) == [d.bus for d in net.der_units]
    assert branch_transformer_indices(net) == []


def test_a_transformer_between_grid_buses_is_a_branch():
    net = _net()
    net.transformers.append(Transformer(name="branch", **BRANCH))
    assert branch_transformer_indices(net) == [len(net.transformers) - 1]
    assert 5 not in unit_transformers(net)


# --- the bug this fixes --------------------------------------------------------
def test_a_branch_transformer_is_in_the_dynamic_model_at_all():
    plain = linearize_network(_net(), run_power_flow(_net()))
    net = _net(**BRANCH)
    with_branch = linearize_network(net, run_power_flow(net))
    assert with_branch.A.shape[0] == plain.A.shape[0] + 2       # its own d/q current
    assert [n for n in with_branch.state_names if "Tr_4" in n] == ["i_{l_d}_{Tr_4}", "i_{l_q}_{Tr_4}"]


def test_the_operating_point_stays_an_equilibrium_at_the_buses_it_joins():
    """It used to leave ~450 pu/s of drift at exactly those two buses, while
    the power flow pushed 24 MW through it."""
    before = _node_residuals(_net())
    net = _net(**BRANCH)
    after = _node_residuals(net)
    assert abs(float(run_power_flow(net).net.res_trafo.p_hv_mw.iloc[-1])) > 10.0   # it really carries power
    for bus in (4, 5):
        assert after[bus] < max(1.0, 2 * before[bus])


def test_it_actually_carries_the_power_the_power_flow_says():
    """The branch current at the operating point, against pandapower's own
    result for the same transformer."""
    net = _net(**BRANCH)
    result = run_power_flow(net)
    model = build_nonlinear_network(net, result)
    idx = branch_transformer_indices(net)[0]
    (ild, ilq) = model.op.transformer_i0[idx]
    vgd, vgq = model.op.node_vg[5]                      # the LV end
    # The branch current is what the transformer *delivers into* the LV bus,
    # while pandapower's p_lv_mw is what that bus draws *into* the transformer:
    # the same power, opposite sign conventions.
    delivered_to_lv = (vgd * ild + vgq * ilq) * net.sn_mva
    assert delivered_to_lv == pytest.approx(-float(result.net.res_trafo.p_lv_mw.iloc[-1]), rel=1e-2)
    assert abs(delivered_to_lv) > 10.0


# --- turns ratio and phase shift ------------------------------------------------
@pytest.mark.parametrize("tap,shift", [(1.0, 0.0), (1.05, 0.0), (0.95, 0.0), (1.0, 5.0), (1.0, -10.0), (1.025, 3.0)])
def test_tap_and_shift_agree_between_power_flow_and_dynamics(tap, shift):
    """The ratio is applied by wiring coefficients; getting them wrong shows up
    as the operating point no longer being an equilibrium."""
    baseline = max(_node_residuals(_net()).values())
    net = _net(tap_ratio=tap, shift_degree=shift)
    assert max(_node_residuals(net).values()) < max(1.0, 3 * baseline)


def test_a_higher_tap_pushes_more_power_through():
    low = float(run_power_flow(_net(tap_ratio=0.95)).net.res_trafo.p_lv_mw.iloc[-1])
    mid = float(run_power_flow(_net(tap_ratio=1.00)).net.res_trafo.p_lv_mw.iloc[-1])
    high = float(run_power_flow(_net(tap_ratio=1.05)).net.res_trafo.p_lv_mw.iloc[-1])
    assert low < mid < high


def test_a_phase_shift_moves_active_power_and_can_reverse_it():
    """What a phase-shifting transformer is for."""
    none = float(run_power_flow(_net(shift_degree=0.0)).net.res_trafo.p_lv_mw.iloc[-1])
    plus = float(run_power_flow(_net(shift_degree=5.0)).net.res_trafo.p_lv_mw.iloc[-1])
    minus = float(run_power_flow(_net(shift_degree=-10.0)).net.res_trafo.p_lv_mw.iloc[-1])
    assert plus > none > 0 > minus
    assert abs(plus - none) > abs(none) * 0.5          # a few degrees is a large effect


def test_the_ideal_ratio_conserves_power():
    """Across the ideal part; the difference is the series resistance's loss."""
    result = run_power_flow(_net(tap_ratio=1.05))
    row = result.net.res_trafo.iloc[-1]
    loss = abs(float(row.p_hv_mw) + float(row.p_lv_mw))
    assert loss < 0.02 * abs(float(row.p_lv_mw))


# --- switching ------------------------------------------------------------------
def test_opening_a_branch_transformer_drops_only_it():
    from g2elin_core.network.breakers import energized_network

    net = _net(**BRANCH)
    n_with = linearize_network(net, run_power_flow(net)).A.shape[0]
    net.transformers[-1].hv_closed = False
    reduced = energized_network(net)
    n_without = linearize_network(reduced, run_power_flow(reduced)).A.shape[0]
    assert n_without == n_with - 2
    assert len(reduced.der_units) == len(net.der_units)     # no unit went with it


# --- validation -------------------------------------------------------------------
def test_two_transformers_on_one_units_lv_bus_are_rejected():
    """This used to re-wire the unit to a different grid bus, silently, because
    the lookup is a dict and the last entry won."""
    net = wscc9_3sm()
    bus = net.der_units[1].bus
    net.transformers.append(Transformer(hv_bus=6, lv_bus=bus, r_pu=0.0, x_pu=0.06, sn_mva=100.0))
    errors = [i.message for i in validate_network(net) if i.severity == "error"]
    assert any("both have their LV side on bus" in m for m in errors)


def test_a_transformer_with_both_sides_on_one_bus_is_rejected():
    net = wscc9_3sm()
    net.transformers.append(Transformer(hv_bus=4, lv_bus=4, r_pu=0.0, x_pu=0.06, sn_mva=100.0))
    assert any("both sides on bus" in i.message for i in validate_network(net) if i.severity == "error")


def test_a_plain_preset_still_validates_clean():
    assert [i for i in validate_network(wscc9_3sm()) if i.severity == "error"] == []
