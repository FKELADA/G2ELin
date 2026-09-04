"""Static power flow regression test for the WSCC 9-bus / 3-SM preset.

There's no MATLAB-exported reference yet (see tests/golden/README.md) — until
one exists, this checks convergence and sanity bounds (voltages within
Power_Fl.m's default +/-5% band, slack picking up the right sign of power,
line/transformer losses positive and small relative to load) rather than
exact values. Tighten this to an exact comparison once bus_sol.json lands.
"""

from __future__ import annotations

import pytest

from g2elin_core.network.presets import wscc9_3sm
from g2elin_core.powerflow import run_power_flow


@pytest.fixture(scope="module")
def result():
    network = wscc9_3sm()
    return run_power_flow(network)


def test_converges(result):
    assert result.converged


def test_voltages_within_default_band(result):
    # Power_Fl.m defaults: Vmin = 0.95, Vmax = 1.05
    table = result.bus_table()
    assert (table["vm_pu"] >= 0.95).all()
    assert (table["vm_pu"] <= 1.05).all()


def test_slack_bus_is_generating(result):
    table = result.bus_table()
    slack_bus = next(d.bus for d in result.network.der_units if d.bus_type.value == "slack")
    slack_row = table.loc[table["bus"] == slack_bus].iloc[0]
    assert slack_row["p_net_gen_mw"] > 0


def test_total_generation_covers_load_plus_losses(result):
    table = result.bus_table()
    total_gen = table["p_net_gen_mw"].clip(lower=0).sum()
    total_load = sum(l.p_mw for l in result.network.loads) + sum(
        d.p_cons_mw for d in result.network.der_units
    )
    losses = result.total_losses_mw()
    assert total_gen == pytest.approx(total_load + losses, abs=0.5)
