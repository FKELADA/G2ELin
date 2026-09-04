# `timeseries` — multi-scenario load flow, no dynamics

This is feature 2.1 from the migration plan: **new** relative to the
MATLAB tool, which had no equivalent. Each `Snapshot` is a set of P/Q
overrides (dispatch and/or load) applied to a base `Network`;
`run_time_series()` solves each one with the existing static power-flow
adapter in a plain loop — no dynamics, no coupling between snapshots.

```{mermaid}
sequenceDiagram
    participant U as caller
    participant S as scale_loads() / dispatch_sweep()
    participant A as apply_snapshot()
    participant PF as run_power_flow()
    participant R as TimeSeriesResult

    U->>S: base Network, {label: factor, ...}
    S-->>U: list[Snapshot]
    loop each Snapshot
        U->>A: Network, Snapshot
        A-->>U: new Network (der_units/loads overridden,<br/>everything else unchanged)
        U->>PF: run_power_flow(overridden Network)
        PF-->>R: bus_table() or empty if not converged
    end
    R-->>U: bus_tables[label], converged[label]
```

This is a deliberate simplicity choice over pandapower's native
`timeseries`/`controller` batch API: it reuses the already-tested
`run_power_flow()` path exactly, at the cost of not getting pandapower's
vectorized C-level speedup. Worth revisiting if sweep size ever makes that
difference matter — swap the inner loop for pandapower's batch
`update_data`/`OutputWriter` machinery without changing the
`Snapshot`/`TimeSeriesResult` interface.

`Network` is a frozen pydantic model, so `apply_snapshot()` never mutates
the base network — `model_copy(update=...)` produces an independent copy
per snapshot, which is what makes running snapshots in any order (or in
parallel) safe.

## Reference

```{eval-rst}
.. automodule:: g2elin_core.timeseries.scenario
```
