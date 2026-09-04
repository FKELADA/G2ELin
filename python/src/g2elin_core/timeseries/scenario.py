"""Time-series load flow: multiple dispatch/load snapshots, no dynamics.

This is feature 2.1 from the migration plan — new relative to the MATLAB
tool. Each :class:`Snapshot` is a set of P/Q overrides (dispatch and/or
load) applied to a base :class:`~g2elin_core.network.schema.Network`, and
:func:`run_time_series` solves each one with the existing static
power-flow adapter (:mod:`g2elin_core.powerflow`) in a plain loop.

That's a deliberate simplicity choice over pandapower's native
``timeseries``/``controller`` batch API: it reuses the already-tested
:func:`~g2elin_core.powerflow.run_power_flow` path exactly, at the cost of
not getting pandapower's vectorized C-level speedup. Worth revisiting if
sweep size ever makes that speed difference matter (see the migration
plan's "pandapower solve speed at scale" risk note) — swap
:func:`run_time_series`'s inner loop for pandapower's batch
``update_data``/``OutputWriter`` machinery without changing the
``Snapshot``/``TimeSeriesResult`` interface.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

from g2elin_core.network.schema import Network
from g2elin_core.powerflow.pandapower_adapter import PowerFlowResult, run_power_flow


@dataclass(frozen=True)
class Snapshot:
    """P/Q overrides for one time step. Unlisted DER units/loads keep the
    base network's values. DER units are keyed by ``DerUnit.id``, loads by
    their index in ``Network.loads`` (loads have no id field).
    """

    label: str
    der_p_mw: dict[int, float] = field(default_factory=dict)
    der_q_mvar: dict[int, float] = field(default_factory=dict)
    load_p_mw: dict[int, float] = field(default_factory=dict)
    load_q_mvar: dict[int, float] = field(default_factory=dict)


def apply_snapshot(network: Network, snapshot: Snapshot) -> Network:
    """Returns a new :class:`Network` with the snapshot's overrides applied."""
    new_ders = [
        der.model_copy(
            update={
                "p_set_mw": snapshot.der_p_mw.get(der.id, der.p_set_mw),
                "q_set_mvar": snapshot.der_q_mvar.get(der.id, der.q_set_mvar),
            }
        )
        for der in network.der_units
    ]
    new_loads = [
        load.model_copy(
            update={
                "p_mw": snapshot.load_p_mw.get(idx, load.p_mw),
                "q_mvar": snapshot.load_q_mvar.get(idx, load.q_mvar),
            }
        )
        for idx, load in enumerate(network.loads)
    ]
    return network.model_copy(update={"der_units": new_ders, "loads": new_loads})


def scale_loads(network: Network, factors: dict[str, float]) -> list[Snapshot]:
    """One snapshot per ``(label, factor)`` pair, scaling every load's P and Q
    by ``factor`` relative to the base network — e.g. for a simple daily
    load curve: ``scale_loads(net, {"00:00": 0.6, "12:00": 1.0, "18:00": 1.15})``.
    """
    return [
        Snapshot(
            label=label,
            load_p_mw={i: ld.p_mw * factor for i, ld in enumerate(network.loads)},
            load_q_mvar={i: ld.q_mvar * factor for i, ld in enumerate(network.loads)},
        )
        for label, factor in factors.items()
    ]


def dispatch_sweep(network: Network, der_id: int, p_mw_values: dict[str, float]) -> list[Snapshot]:
    """One snapshot per ``(label, P)`` pair, overriding a single DER unit's
    dispatch — e.g. a generator ramp/PV-curve sweep.
    """
    return [Snapshot(label=label, der_p_mw={der_id: p_mw}) for label, p_mw in p_mw_values.items()]


@dataclass
class TimeSeriesResult:
    network: Network
    snapshots: list[Snapshot]
    bus_tables: dict[str, pd.DataFrame]  # snapshot label -> bus_table()
    converged: dict[str, bool]
    # Full per-snapshot result (the pandapower net + everything PowerFlowResult
    # exposes -- res_line/res_trafo/etc.), for callers that want more than the
    # bus table `bus_tables` already gives (e.g. g2elin_api mirroring the
    # Power Flow tab's other-result-table dropdown per snapshot).
    results: dict[str, PowerFlowResult] = field(default_factory=dict)

    def combined_bus_table(self) -> pd.DataFrame:
        """All snapshots' bus tables stacked, with a ``snapshot`` column."""
        frames = []
        for label, table in self.bus_tables.items():
            t = table.copy()
            t.insert(0, "snapshot", label)
            frames.append(t)
        return pd.concat(frames, ignore_index=True)

    def all_converged(self) -> bool:
        return all(self.converged.values())


def run_time_series(network: Network, snapshots: list[Snapshot]) -> TimeSeriesResult:
    bus_tables: dict[str, pd.DataFrame] = {}
    converged: dict[str, bool] = {}
    results: dict[str, PowerFlowResult] = {}
    for snap in snapshots:
        result = run_power_flow(apply_snapshot(network, snap))
        converged[snap.label] = result.converged
        bus_tables[snap.label] = result.bus_table() if result.converged else pd.DataFrame()
        results[snap.label] = result
    return TimeSeriesResult(
        network=network, snapshots=snapshots, bus_tables=bus_tables, converged=converged, results=results
    )
