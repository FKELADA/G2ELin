"""Static power flow regression test for the CIGRE islanded 1SM+2GFM+1GFL preset.

Same status as the WSCC test: sanity bounds only, no MATLAB-exported
reference yet (see tests/golden/README.md). This case exercises what WSCC-9
doesn't: two voltage levels behind a mesh (not radial) MV feeder, an
islanded system with no external grid connection, and a heterogeneous DER
mix (SM slack + 2 PV-bus GFM + 1 PQ-dispatched GFL).
"""

from __future__ import annotations

import pytest

from g2elin_core.network.presets import cigre_islanded_1sm_2gfm_1gfl
from g2elin_core.powerflow import run_power_flow


@pytest.fixture(scope="module")
def result():
    network = cigre_islanded_1sm_2gfm_1gfl()
    return run_power_flow(network)


def test_converges(result):
    assert result.converged


def test_voltages_within_default_band(result):
    # Power_Fl.m defaults: Vmin = 0.95, Vmax = 1.05
    table = result.bus_table()
    assert (table["vm_pu"] >= 0.95).all()
    assert (table["vm_pu"] <= 1.05).all()


def test_pv_buses_hold_voltage_setpoint(result):
    table = result.bus_table()
    for der in result.network.der_units:
        if der.bus_type.value == "pv":
            vm = table.loc[table["bus"] == der.bus, "vm_pu"].iloc[0]
            assert vm == pytest.approx(der.v_set_pu, abs=1e-6)


def test_slack_covers_the_shortfall(result):
    # 2 GFM + 1 GFL dispatch 3 x 0.7 pu x 2.5 MVA = 5.25 MW; total load is
    # ~7.235 MW, so the islanded SM slack must supply the ~2 MW shortfall
    # plus losses rather than sitting near zero (its own Pgen dispatch is 0).
    table = result.bus_table()
    slack_bus = next(d.bus for d in result.network.der_units if d.bus_type.value == "slack")
    slack_p = table.loc[table["bus"] == slack_bus, "p_net_gen_mw"].iloc[0]
    assert slack_p > 1.0


def test_total_generation_covers_load_plus_losses(result):
    table = result.bus_table()
    total_gen = table["p_net_gen_mw"].clip(lower=0).sum()
    total_load = sum(l.p_mw for l in result.network.loads) + sum(
        d.p_cons_mw for d in result.network.der_units
    )
    losses = result.total_losses_mw()
    assert total_gen == pytest.approx(total_load + losses, abs=0.05)
