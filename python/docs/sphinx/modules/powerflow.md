# `powerflow` — static power flow

Replaces `Functions/Power_Fl.m` / `Functions/Load_flow.m` (a hand-rolled
sparse Newton-Raphson solver) with [pandapower](https://www.pandapower.org/).
pandapower's `gen`/`ext_grid` elements give native PV-bus (voltage-controlled)
and slack-bus behavior, so — unlike
[power-grid-model](https://github.com/PowerGridModel/power-grid-model), the
other candidate considered in the migration plan — no outer Q-iteration loop
is needed to hold PV buses at their voltage setpoint.

## Build → solve → read back

```{mermaid}
sequenceDiagram
    participant N as Network
    participant B as build_pandapower_net()
    participant PP as pandapower (pp.runpp)
    participant R as PowerFlowResult

    N->>B: buses, lines, transformers, loads, der_units
    Note over B: per-unit -> physical (Ohm/Siemens/MVA)<br/>conversion at each bus's own voltage base
    B->>PP: pp.create_bus / create_line_from_parameters /<br/>create_transformer_from_parameters / create_load /<br/>create_ext_grid / create_gen / create_sgen
    B-->>R: net, bus_index (Network id -> pandapower internal index)
    R->>PP: pp.runpp(calculate_voltage_angles=True, ...)
    alt converged
        PP-->>R: res_bus, res_line, res_trafo
        R->>R: bus_table() re-indexes back to Network bus ids
    else LoadflowNotConverged
        PP-->>R: converged = False
    end
```

Each `DerUnit`'s `bus_type` selects which pandapower element it becomes:
`SLACK` → `ext_grid` (fixed voltage & angle), `PV` → `gen` (fixed P and |V|,
Q free within a wide default limit matching `Power_Fl.m`'s unspecified
`qg_max`/`qg_min`), anything else → `sgen` (fixed P and Q — a
PQ-dispatched unit, e.g. a GFL not under voltage control).

## Physics: the AC power-flow problem

pandapower solves the standard nonlinear AC power-flow problem via
Newton-Raphson internally (not reimplemented here — this project's own
contribution is the `Network` → pandapower translation and reading
results back, not the solver itself). For completeness: given the bus
admittance matrix $Y_{bus}$ built from every line/transformer's
$R+jX$ (and line charging $B$), each bus $i$'s complex power injection
must satisfy

$$
S_i = P_i + jQ_i = \bar V_i \sum_{k} Y_{ik}^{*}\,\bar V_k^{\,*}
$$

with three bus-type behaviors selected by each `DerUnit.bus_type`
(exactly mirroring `Power_Fl.m`'s own `bus_type` convention, now
pandapower's instead of a hand-rolled solver):

- **slack** (`ext_grid`): $V_i$ and $\theta_i$ fixed (here $\theta=0$,
  the network's own angle reference), $P_i$/$Q_i$ solved for — the slack
  absorbs whatever active/reactive mismatch every other bus's fixed
  dispatch leaves.
- **PV** (`gen`): $P_i$ and $|V_i|$ fixed, $Q_i$ and $\theta_i$ solved for
  (within the wide default $Q$ limits $\pm 9999$ MVAr, matching
  `Power_Fl.m`'s unspecified `qg_max`/`qg_min`) — a voltage-controlled
  unit.
- **PQ** (`sgen`): $P_i$ and $Q_i$ both fixed, $|V_i|$ and $\theta_i$
  solved for — a unit dispatched at a fixed power factor, not under
  voltage control (e.g. a GFL).

Newton-Raphson linearizes the real/reactive mismatch equations
$\Delta P = P_{spec}-P_{calc}(V,\theta)$, $\Delta Q =
Q_{spec}-Q_{calc}(V,\theta)$ around the current iterate via the power-flow
Jacobian $J = \partial(\Delta P,\Delta Q)/\partial(\theta,V)$ and solves
$J\,\delta = -[\Delta P;\Delta Q]$ repeatedly until $\delta$ is smaller
than a tolerance — `pp.runpp`'s own default is 10 iterations at $10^{-8}$
MVA mismatch, which every preset/generated network in this codebase
converges well inside.

**Per-unit → physical conversion**, since pandapower's own line/
transformer model wants physical (Ω, S, MVA) parameters rather than this
project's per-unit convention: given a line's own base impedance
$Z_b = U_n^2/S_n$ at its bus's nominal voltage,

$$
R_{\Omega/km} = \frac{r_{pu}\,Z_b}{\ell}
\qquad
X_{\Omega/km} = \frac{x_{pu}\,Z_b}{\ell}
\qquad
C_{nF/km} = \frac{b_{pu}}{Z_b\,\omega_n\,\ell}\times10^{9}
$$

($\ell=$ `length_km`; note $r_{pu}$/$x_{pu}$/$b_{pu}$ are each the
line's *total* per-unit impedance/susceptance regardless of $\ell$ — the
$\ell$ in the numerator and pandapower's own per-km model cancel exactly,
so `length_km` is bookkeeping metadata here, not a second independent
electrical parameter). A transformer's short-circuit voltage percentages
follow directly from its own $r_{pu}$/$x_{pu}$: $v_{k\%} =
100\sqrt{r_{pu}^2+x_{pu}^2}$, $v_{kr\%}=100\,r_{pu}$.

## Reference

```{eval-rst}
.. automodule:: g2elin_core.powerflow.pandapower_adapter
```
