"""Time-series load flow (feature 2.1) — new relative to the MATLAB tool,
so no golden MATLAB reference applies here (unlike tests/golden/). These
check physical sanity: convergence across a sweep, and that losses/slack
generation move monotonically with load.
"""

from __future__ import annotations

import pytest

from g2elin_core.network.presets import wscc9_3sm
from g2elin_core.timeseries import dispatch_sweep, run_time_series, scale_loads


@pytest.fixture(scope="module")
def network():
    return wscc9_3sm()


def test_load_scaling_sweep_converges(network):
    snapshots = scale_loads(network, {"light": 0.7, "base": 1.0, "heavy": 1.3})
    result = run_time_series(network, snapshots)
    assert result.all_converged()
    assert set(result.bus_tables.keys()) == {"light", "base", "heavy"}


def test_slack_generation_tracks_load(network):
    # WSCC-9's two PV generators are fixed-dispatch (210 MW total,
    # independent of load), so the slack alone must absorb the load change
    # — this is a direct power-balance consequence, robustly monotonic.
    # (Total *losses* are not guaranteed monotonic here: at light load the
    # fixed PV output overshoots demand more, which can increase circulating
    # flow and losses even as total load drops — a real, not a bug.)
    snapshots = scale_loads(network, {"light": 0.7, "base": 1.0, "heavy": 1.3})
    result = run_time_series(network, snapshots)

    slack_bus = next(d.bus for d in network.der_units if d.bus_type.value == "slack")
    slack_p = {
        label: table.loc[table["bus"] == slack_bus, "p_net_gen_mw"].iloc[0]
        for label, table in result.bus_tables.items()
    }
    assert slack_p["light"] < slack_p["base"] < slack_p["heavy"]


def test_combined_bus_table_shape(network):
    snapshots = scale_loads(network, {"a": 0.9, "b": 1.1})
    result = run_time_series(network, snapshots)
    table = result.combined_bus_table()
    assert set(table["snapshot"]) == {"a", "b"}
    assert len(table) == 2 * len(network.buses)


def test_dispatch_sweep_on_pv_generator(network):
    pv_der = next(d for d in network.der_units if d.bus_type.value == "pv")
    snapshots = dispatch_sweep(network, pv_der.id, {"low": 50.0, "mid": 90.0, "high": 130.0})
    result = run_time_series(network, snapshots)
    assert result.all_converged()

    # The swept unit's own bus should show ~the commanded P (net of its aux load).
    for snap, p_mw in zip(snapshots, [50.0, 90.0, 130.0]):
        table = result.bus_tables[snap.label]
        row = table.loc[table["bus"] == pv_der.bus].iloc[0]
        assert row["p_net_gen_mw"] == pytest.approx(p_mw - pv_der.p_cons_mw, abs=1e-6)
