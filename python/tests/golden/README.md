# Golden fixtures

These tests are meant to compare the Python port's output against the
original MATLAB tool's output, bus for bus. That comparison isn't wired up
yet because it needs a `bus_sol` export from a real MATLAB run.

## How to add the real golden reference

1. In MATLAB, run `WSCC/script_WSCC.m` with `model.name = 'WSCC_3SM'`.
2. After `Load_flow` runs, export `bus_sol` and `P_loss`, e.g.:

   ```matlab
   result.bus_sol = bus_sol;
   result.P_loss = P_loss;
   savejson('', result, fullfile(baseDirectory, 'python', 'tests', 'golden', 'wscc9_3sm_bus_sol.json'));
   % or, without the JSON toolbox:
   writematrix(bus_sol, fullfile(baseDirectory, 'python', 'tests', 'golden', 'wscc9_3sm_bus_sol.csv'));
   ```

3. Replace the sanity-bound assertions in `test_wscc9_powerflow.py` with an
   exact comparison (voltage magnitude/angle to ~1e-4, P/Q flows to
   ~1e-3 pu) against the exported values, accounting for the bus
   renumbering documented in `g2elin_core.network.presets.wscc9_3sm`.

The same applies to `test_cigre_powerflow.py`: run
`CIGRE/script_CIGRE_Islanded.m` with `model.name = 'CIGRE_Islanded_1SM_2GFM_1GFL'`
and export `bus_sol` the same way. The CIGRE interconnected case (`S0=1`,
with an infinite-bus upstream connection) isn't ported yet and would be a
good third fixture — it exercises a slack that's an infinite bus rather
than a synchronous machine.

## Modal analysis (test_wscc9_modal_analysis.py)

Same situation: run `WSCC/script_WSCC.m` with `model.name = 'WSCC_3SM'` and
`ss_analysis = 1`, then after `modal_analysis(...)` runs, export `eigA`
(or `d`/`eigenVal` inside `modal_analysis.m`) and `stateNames` to compare
eigenvalues and mode participation against `analyze(system.A, system.state_names)`.
Expect the state *ordering* to differ from the MATLAB run (this port
orders states slack-DG, other-DGs, nodes, lines, loads — cosmetically the
same as `script_generic.m`, but component-internal state order should be
double-checked against `states_SG`) — compare the eigenvalue *set*
(sorted) first, then match states by name rather than by position.
