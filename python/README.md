# g2elin-core

Python core for **G2ELin**, porting the original MATLAB/Simulink toolbox to
Python + a web front end. See the migration plan for the full picture
(architecture, feature scope, phased roadmap). `g2elin_core` is the compute
core — no web dependencies, usable standalone. `g2elin_api` (optional,
`pip install -e ".[api]"`) is a thin HTTP layer over it, serving a static
web frontend from `web/`.


## Selectable unit models

A unit's regulators are part of its description, not baked into the
component. Each model brings its *own* parameters under its own names, so
what is typed is what that model's diagram shows, and an override under a
name the fitted model does not have is refused rather than ignored.

| Unit | Field | Models |
|---|---|---|
| Synchronous machine | `exciter` | `g2elin` (transducer, amplifier, exciter, rate feedback), `kundur` (thyristor + transient gain reduction, Fig. E12.9) |
| | `pss` | `g2elin` (input filter, washout, two lead-lags), `kundur` (washout, two lead-lags), `none` |
| | `governor` | `g2elin` (droop into a first-order lag), `none` |
| Grid-forming converter | `controller` | `droop`, `droop_filtered`, `dvoc`, `vsm`, `matching` |

Every default is the model the component always had, so any saved network
linearises to exactly what it did before.

`none` is equipment that is not fitted rather than a gain turned down: the
states, the parameters and the reduction group all go with it. Two
consequences worth knowing, both of which the tools handle explicitly:

- A machine's `P_ref` is its *governor's* setpoint. With no governor its
  column of `B` is exactly zero, so a step on it does nothing — and the step
  response says so rather than drawing an unexplained flat line.
- A fleet with no governor anywhere has no frequency regulation, so its
  common frequency is a free integrator. That is a marginal direction of the
  model, not an instability, and it is set aside from the stability verdict
  beside the reference angles (`Network.frequency_is_regulated()`).

The grid-forming laws are `symGFM_types.m`'s, with `script_generic.m`'s
gains — each written in terms of the droop tuning (`eta = mp`, `J =
1/(mp*wf)`, `K_theta = mp*Kpdc`), so the laws are *comparable* and a study
that swaps one for another sees what the law changes rather than what a
different tuning changes. VSM and matching carry no power filters at all.

See `docs/sphinx/design/controller-models.md` for the equations, the block
diagrams, and the plan for the rest of the IEEE 421.5 library.

## Status

Phase P0 (foundations), P1 (static power flow), P2 (linear small-signal /
modal analysis — SM/GFM/GFL/IB slack, plus the full toolbox:
sensitivity, mode shape, free/step response), P3 (time-series load flow),
P4 (EMT/nonlinear time-domain simulation — now with an analytic Newton
Jacobian, ~9x faster than the first working slice), and P7 (web UI — result
tabs, always-light theme, a network editor and drag-and-drop from-scratch
builder over an arbitrary `Network` JSON, not just the fixed presets):

- `g2elin_core.network` — typed network schema (buses, lines, transformers,
  DER units), replacing the MATLAB `Y_network`/`Y_line`/`Y_TR`/`Y_DER`
  numeric-column convention with named, validated fields. Being a pydantic
  model (not a dataclass) is what lets it round-trip to/from JSON for
  free — the basis for the web UI's network editor/builder, further down.
- `g2elin_core.network.presets` — 15 preset cases in 4 families (5 in the
  original slice, described below; 10 more DER-mix variants on the WSCC-9
  and CIGRE-islanded topologies, added later — see "Web UI (P7)" further
  down):
  - `wscc9_3sm()` — WSCC 9-bus, 3 synchronous machines (`Functions/WSCC_raw.m`
    + `WSCC/script_WSCC.m`, case `'WSCC_3SM'`). Transcribed line-by-line.
  - `cigre_islanded_1sm_2gfm_1gfl()` — CIGRE MV benchmark feeder, islanded,
    1 SM + 2 GFM + 1 GFL (`Functions/CIGRE_raw.m` +
    `CIGRE/script_CIGRE_Islanded.m`, case `'CIGRE_Islanded_1SM_2GFM_1GFL'`
    — matches the `.slx` model actually present in this repo). Transcribed
    line-by-line; used for power flow, modal analysis, and EMT simulation
    alike.
  - `sm_smib()` / `gfm_smib()` / `gfl_smib()` — single-machine-infinite-bus
    presets: one SM/GFM/GFL against an infinite bus. **Not** a straight
    port like the two above — `Functions/preset_networks.m` names exactly
    these three cases but their driving scripts were never finished
    (empty `switch` bodies), so the topology/line/load values are real
    (from `Functions/SMIB_raw.m`, confirmed to share CIGRE's own
    base-value convention numerically) while the DER attachment is a
    documented reconstruction — see `sm_smib()`'s own docstring. These are
    also what finally exercises `g2elin_core.components.ib` (the
    infinite-bus component, built early in this project but never wired
    into the interconnection until these presets needed a real slack for
    it — see the interconnect bullet below).
- `g2elin_core.network.topology` — a 2D layout (`networkx`'s Kamada-Kawai,
  over the bus graph) for the web UI's Network and Power Flow tabs. A DER
  unit or load isn't a separate geometric node in this schema (a DER sits
  on its own real `Bus`, a load just references an existing bus id), so
  the graph to lay out is exactly the buses + lines/transformers.
- `g2elin_core.powerflow` — static power flow via **pandapower** (see the
  plan for why pandapower was chosen over power-grid-model: native PV-bus
  voltage control, no outer Q-iteration loop needed).
- `g2elin_core.components` — symbolic (sympy) nonlinear-then-linearized
  small-signal models, ported from `Functions/sym*.m`: synchronous machine
  (`sm.py`, from `symSG.m`/`SG_subs.m`), grid-forming converter with Droop
  control (`gfm.py`, from `symGFM_types.m`/`GFM_subs.m`), grid-following
  converter (`gfl.py`, from `symGFL.m`/`GFL_subs.m`), infinite bus
  (`ib.py`), line (`line.py`), node (`node.py`), constant-impedance load
  (`load.py`). GFM's other 5 outer-control variants
  (`Droop+filter`/`Droop+PLL`/`dVOC`/`VSM`/`Matching`) aren't ported yet —
  only `Droop`, which is what both presets use.
  `ComponentDAE.nonlinear_funcs()` (base.py) now also exposes the
  *nonlinear* `f`/`g`/`h` as lambdified numpy callables, alongside the
  linearized `A`/`B`/`C`/`D` — both come from the same symbolic equations
  (see the migration plan's "the nonlinear models already exist" finding).
  This is the foundation P4 (EMT simulation) needs; it isn't wired into a
  time-domain integrator yet, but it's already doing real work as an
  independent correctness check (see below).
- `g2elin_core.pu_base` — per-unit AC/DC base-value calculator, ported from
  `Functions/PU_calc.m` (needed for GFM/GFL's DC-link parameters).
- `g2elin_core.operating_point` — computes each component's linearization
  point from a solved power flow (ports `script_generic.m`'s "SM/GFM/GFL
  Initializations" sections + the node/line/load operating-point loops).
- `g2elin_core.interconnect` — assembles per-component state-spaces into one
  closed-loop system. Functionally ports `Functions/assoc_matrices.m`, but
  via explicit named-port wiring rules instead of replicating its ~300
  lines of index arithmetic — see the module docstring for the mapping.
  `assemble_network` handles any mix of SM/GFM/GFL DER units, with the
  slack being either a synchronous machine or an infinite bus
  (`_SLACK_KIND_BY_UNIT_TYPE` picks the right port kind for either).
  **Block names are per-type counters**, not the DER's raw id — `"SM_1"`,
  `"GFM_1"`, `"GFM_2"`, `"GFL_1"` rather than `"DG_1"`/`"DG_2"`/`"DG_3"`/
  `"DG_4"` — so every state/input/output name built by interpolating a
  block's name (`state_names`, `input_names`, `output_names`, in both the
  linear and nonlinear models, since both are built by the same
  `build_blocks_and_wiring()`) says what *kind* of unit it is without
  cross-referencing the network definition, e.g. `dw_r_{SM_2}` rather than
  `dw_r_{DG_2}`.
- `g2elin_core.modal` — eigenvalues + participation-factor matrix from the
  numeric core of `Functions/modal_analysis.m`, **plus the rest of that
  file's toolbox** (`modal.toolbox`, read in full to get this right, not
  guessed): `eigenvalue_sensitivity()` (masked to A's structurally-nonzero
  entries), `mode_shape()` (top-5 participating states' phase angle),
  `free_response()` (closed-form modal-expansion response to an
  initial-condition perturbation — no ODE solve needed), and
  `step_response()` (MIMO step response via a SISO reduction +
  `scipy.signal.step`). Each cross-checked against an independent
  reference in `tests/test_modal_toolbox.py` (`free_response` against
  `scipy.linalg.expm(A*t) @ x0` directly, `step_response` against
  `scipy.signal.step` on a hand-built SISO system), not just "runs and is
  finite."
- `g2elin_core.pipeline.linearize_network` — ties the above together:
  solved power flow -> operating point -> linearized components -> closed
  loop.

Two MATLAB quirks are preserved rather than silently fixed (see the
`g2elin_core.components.sm` module docstring): every unit's transformer
impedance and every node's shunt susceptance come from the *first*
DG's/line's values in `script_generic.m`, not each one's own.

Validated on both presets — WSCC-9/3SM (SM-only) and CIGRE islanded
1SM+2GFM+1GFL (the first to exercise GFM/GFL in the interconnection): state
counts match expectations, no NaN/Inf, every eigenvalue stable except one
at the expected zero (the unreferenced global angle — no infinite bus in
either case), and the mode structure is exactly what each model should
produce: SM electromechanical modes (0.7-2.2 Hz), a GFM-GFM oscillation
mode (~6 Hz on CIGRE), AVR modes (5-8 Hz), a GFL PLL mode (~13.5 Hz on
CIGRE), GFM voltage-loop modes matching their tuned bandwidth (~45 Hz), and
network/electrical modes clustered at each case's own synchronous frequency
(+/-60 Hz on WSCC, +/-50 Hz on CIGRE). No MATLAB-exported reference to diff
against exactly yet (see `tests/golden/README.md`).

While building the CIGRE case, a real bug surfaced and was fixed: the
participation-factor code computed right eigenvectors from `eig(A)` and left
eigenvectors from a *separate* `eig(A.conj().T)` call, then assumed they
came back in matching order — nothing guarantees that, and it silently
produced NaN participation factors under exactly the conditions CIGRE has
(two near-identical GFM units). Fixed by computing both from a single
`scipy.linalg.eig(A, left=True, right=True)` call.

- `g2elin_core.timeseries` — feature 2.1 (time-series load flow, new
  relative to the MATLAB tool): a `Snapshot` of P/Q overrides applied to a
  base `Network`, solved in a plain loop over `powerflow.run_power_flow`.
  Deliberately simpler than pandapower's native batch API for now — see the
  module docstring for the tradeoff and how to swap in the faster path
  later if sweep size ever makes it matter.

### A second, more consequential bug — found by cross-validating linear against nonlinear

`tests/test_nonlinear_cross_validation.py` evaluates each component's
*nonlinear* `f`/`g` at its own linearization operating point (should be
~0), and separately checks that a finite-difference Jacobian of that same
nonlinear `f`/`g` reconstructs the symbolic `A` matrix `linearize_*()`
produces — two different code paths (symbolic differentiation +
substitution vs. lambdified functions + finite differences) that have to
agree if the implementation is right.

This caught a real bug in `SmOperatingPoint`: the machine's own dq-frame
voltage was computed as `v0 * exp(j*(angle_terminal_rad - theta0))`, which
double-counts `angle_terminal_rad` — `v0` already equals
`v_terminal_pu * exp(j*angle_terminal_rad)`, so it should have been just
`v0 * exp(-j*theta0)`. It silently corrupted `ved0`/`veq0` (and everything
derived from them: flux linkages, AVR/PSS initial states, the linearized
`A`/`B`/`C`/`D`) for **every non-slack synchronous machine** — WSCC's two
PV-dispatched SMs, in particular. It was invisible on the slack unit only
because the power-flow slack reference pins that bus's angle to exactly 0,
making the bug's contribution vanish by coincidence. Every WSCC-9 modal
analysis result reported earlier in this project was computed with this
bug present; it's now fixed in `components/sm.py`, and both golden modal
tests were re-run clean afterward (their structural checks — stability,
electromechanical modes present — still hold, but the exact eigenvalues
shifted).

Chasing this down also surfaced three legitimate (non-bug) reasons a
nonlinear residual isn't exactly zero at these operating points — a
rotating machine's absolute angle state has nonzero derivative by
definition, floating-point rotation round-off gets amplified by a large
`wb/Lt` coefficient in SM's current-state equations, and the "every unit
uses DG#1's transformer impedance" quirk (see below) creates a small real
mismatch for every other unit — all documented in detail in that test
file's `assert_nonlinear_matches_linear` docstring, since future
maintainers will hit the same "why isn't this exactly zero" question.

### P4: a first working EMT (nonlinear) time-domain integrator

A naming note, since this changed after the fact: this integrates the
network's actual nonlinear differential-algebraic equations (not a
linearized approximation) — a dq-frame, averaged-converter, single-
rotating-reference-frame formulation, so no switching-level converter
physics, no abc/unbalanced-fault representation, and lines/transformers
are lumped RL branches rather than distributed/traveling-wave models. That
still integrates real nonlinear time-domain dynamics, which is what
"EMT simulation" means for this project's purposes, so it's named and
scoped that way here rather than as a separate "RMS" fidelity tier.

- `g2elin_core.timedomain` — integrates the network's *nonlinear* DAE
  forward in time, reusing the same nonlinear `f`/`g`/`h` callables the
  cross-validation above already exercises, and the same interconnection
  topology (`interconnect.build_blocks_and_wiring`/`compute_topology`) the
  linear path uses — refactored out of `interconnect/assemble.py` so both
  paths share it rather than the nonlinear side re-deriving the wiring.
  `build_nonlinear_network(network, result)` mirrors
  `pipeline.linearize_network`; `simulate(model, t_span, ...)` integrates
  via `scipy.integrate.solve_ivp` (`Radau` — the stiffest modes reach
  ~1e6-1e7 rad/s, see the SM operating-point note above).

  The real complexity here: a component's inputs are wired to *other*
  components' outputs, which can themselves depend on those same inputs
  (`u = F@u_exo + G@h(x, z, u)`), so every component's algebraic variables
  and the full interconnection have to be solved *together* by one
  coupled Newton system at every state the integrator asks about — not
  component-by-component. See `timedomain/emt.py`'s module docstring for
  the exact formulation.

  Building this immediately re-surfaced the transformer-impedance quirk
  from a new angle: solving that Newton system at the "operating point"
  used for linearization does *not* return that same point back — `u_sol`
  differs from the per-component operating-point guess by ~0.08 (small
  relative to per-unit scale, but not solver noise), because the "every
  node uses line #1's susceptance" quirk means the per-component operating
  points were never quite a network-wide KCL fixed point to begin with.
  That's invisible to linear modal analysis (a Jacobian is valid at any
  point) but real for nonlinear simulation, which starts with an actual
  transient to settle rather than a clean equilibrium. Documented (not
  asserted away) in `test_emt_simulation.py`.

  One implementation-robustness finding worth flagging for whoever extends
  this: SciPy's default derivative-free Newton method (`"hybr"`) for the
  coupled algebraic solve converges fine exactly at the operating point but
  reliably fails for the nearby points an ODE integrator's step-size and
  collocation probing evaluate. Levenberg-Marquardt (`"lm"`) was robust in
  the same situations, at the cost of estimating its own Jacobian
  numerically every iteration.

  **Update:** that Jacobian is no longer numeric. `Gz`/`Gx`/`Gu`/`Hz`/`Hx`/
  `Hu` — the same per-component symbolic Jacobians `linearize()` already
  used, just not yet lambdified into numeric callables — now are
  (`ComponentDAE.nonlinear_jacobians()`), and the coupled residual's own
  Jacobian is assembled from them the same way the linear A/B/C/D
  elimination combines per-component Jacobians with the interconnection's
  `G` matrix (`NonlinearNetworkModel._residual_jacobian`). With an exact
  Jacobian, `"hybr"` became viable again — and roughly **9x faster** than
  `"lm"`'s numeric one: a single algebraic solve went from ~0.05s to
  ~0.006s. Simulating WSCC-9 for a full second (previously not attempted —
  even 1e-3s wasn't validated) now takes ~13s wall time, no NaNs or
  blow-up. `theta`/`theta_pll` legitimately grow to ~377 (≈ `wb`*t) over
  that second — not instability, just the absolute rotor-angle states
  doing what they're defined to do (see `test_emt_simulation.py`'s
  `_non_rotating_mask` note) — everything else stays bounded.

### P7: a first slice of the web UI

- `g2elin_api` (`src/g2elin_api/`) — a FastAPI app wrapping the real
  `g2elin_core` functions (no mocking) over the five validated presets:
  `POST /api/presets/{id}/powerflow`, `/modal` (+ `/modal/sensitivity`,
  `/modal/mode_shape`, `/modal/free_response`, `/modal/step_response`),
  `/timeseries`, `/emt`, plus `GET /api/presets`, `GET
  /api/presets/{id}/topology`, `GET /api/presets/{id}/states` for the
  picker/diagram/pickers. Response models (`g2elin_api/schemas.py`) are
  deliberately separate from `g2elin_core`'s own result types
  (`PowerFlowResult` wraps a live pandapower net, `ModalAnalysisResult`
  holds numpy arrays — neither is JSON-serializable), so the wire format
  is a considered choice, not whatever `dict()` happened to produce.
- `web/index.html` — a single-page frontend: a preset picker, and (now)
  six tabs (Power Flow, Modal Analysis, Time Series, Network, EMT
  Simulation, Documentation) that call the API and
  render results, including an eigenvalue map (SVG, symlog-scaled real
  axis and frequency axis since eigenvalues here span ~15 orders of
  magnitude — see the Rg-penalty-parameter note above — with points
  colored by stability status per the dataviz skill's status-palette
  convention, plus the constant-damping-ratio guide lines
  `Functions/modal_analysis.m` also draws) alongside the mode table.
  Always light — the dark-mode CSS block that used to shadow-follow
  `prefers-color-scheme` was removed outright, so every chart's color
  logic only reasons about one palette.

  **Deliberately plain JS, not React** — the migration plan's original
  target stack is React/TypeScript via Vite, but this environment has no
  Node.js/npm, so there's no way to install or build one. Rather than
  block on tooling this session doesn't have, this ships as a build-free
  page FastAPI serves directly as a static file — real, working, and
  swappable for a proper Vite+React+TS app later without touching the API
  contract it calls. Treat this as scaffolding to replace, not a style
  decision to preserve.

  **Not personally visually verified.** This environment also has no
  browser-automation tooling (no `chromium-cli`, no Playwright/Node), so
  unlike everything else in this project, this wasn't verified by running
  it and checking the result — only by starting the real server and
  confirming correct HTTP-level behavior (`tests/test_api.py`: every
  endpoint, status codes, response shapes; the static file serves with
  the right content-type and contains the expected markup) and by careful
  manual review of the JS against those confirmed response shapes. Opened
  in the user's own browser via `Start-Process` for them to check what I
  couldn't.

  No network editor or scenario builder yet — this is "run orchestration +
  results dashboards" over fixed presets (P7's original scope also
  includes building/editing arbitrary networks interactively, which is
  substantially more work and deferred).

- **EMT Simulation tab**, wrapping P4 (`g2elin_core.timedomain`) over
  HTTP: two new endpoints (`GET /api/presets/{id}/states`, `POST .../emt`).
  `/states` exists because the nonlinear model's state names are entirely
  preset-dependent (which DER id is the slack, which unit types exist), so
  the tab's state pickers are populated from it rather than hard-coded —
  the same problem `find_state_index`'s substring matching solves
  internally, which `/emt` reuses server-side so a typo'd/ambiguous state
  name comes back as a 422 with that function's own error message, not a
  silently wrong answer. `t_final` is bounded server-side (not just
  documented) since this endpoint runs a real nonlinear DAE integration
  synchronously inside one HTTP request.

  The **EMT Simulation** tab perturbs a chosen state and plots the
  resulting trajectory as a multi-series line chart — built directly via
  the SVG DOM API (not the `innerHTML` string templates the other tabs
  use), since it needed live `mousemove` listeners for a real
  crosshair+tooltip (per the dataviz skill: not optional for a line
  chart).

  Verified the same way as the rest of P7: `tests/test_api.py` exercises
  the new endpoints (happy path, an ambiguous/unknown `perturb_state`, an
  out-of-bounds `t_final`) over HTTP, plus the frontend was reviewed by
  hand against those confirmed response shapes; not personally visually
  verified in a browser, same no-browser-automation constraint as the
  rest of this tab set.

- **Sphinx documentation, integrated as a tab.** `docs/sphinx/`
  builds full API-reference + architecture documentation via
  [Sphinx](https://www.sphinx-doc.org/) autodoc — one page per module/
  subpackage, each with a hand-written [Mermaid](https://mermaid.js.org/)
  diagram (class diagrams for the data structures, flowcharts/sequence
  diagrams for the non-obvious control flow — the coupled Newton solve in
  `timedomain`, the EMT subsystem hierarchy)
  plus the real docstrings pulled straight from the code, so the diagrams
  and the reference can't drift apart silently. Build it with
  `pip install -e ".[docs]"` then `python tools/build_docs.py`; `g2elin_api`
  mounts the output at `/manual` (deliberately not `/docs`, which is
  FastAPI's own Swagger UI) if it's been built, and `web/index.html`'s new
  **Documentation** tab lazily loads it in an iframe on first click, with a
  plain "run this command" message instead of a broken iframe if it hasn't
  been built yet. Mermaid diagrams render client-side via a CDN-hosted
  `mermaid.js` (loaded by `sphinxcontrib-mermaid`), so no Node.js/
  mermaid-cli is needed at build time — consistent with the rest of this
  environment's no-Node.js constraint. Verified the same way as the rest of
  P7: `tests/test_api.py` checks the `/manual` mount and the tab's markup
  over HTTP, and the build itself was run and inspected (0 Sphinx warnings,
  autodoc content spot-checked in the generated HTML) rather than just
  assumed to work; not personally visually verified in a browser for the
  same no-browser-automation reason as the rest of the web UI, per the note
  above — opened via `Start-Process` for the user to check.

- **Network tab, and a network diagram on Power Flow.** New
  `g2elin_core.network.topology` (Kamada-Kawai layout over the bus graph)
  + `GET .../topology`, rendered by a shared `renderNetworkDiagram()`
  used two ways: on the **Network** tab, colored by unit type (SM/GFM/
  GFL/IB, fixed categorical order); on **Power Flow**, colored by voltage
  magnitude (sequential blue, fixed [0.9, 1.1] pu domain so a color means
  the same thing across presets, not just within one) above the results
  table. Power Flow also gained a dropdown over pandapower's other result
  tables (`res_line`/`res_trafo`/`res_load`/`res_gen`/`res_sgen`/
  `res_ext_grid` — `PowerFlowResult` gained one accessor per table,
  joined back to this project's own bus ids).

  **Time Series got the same diagram-plus-dropdown treatment**, one
  instance per snapshot: `g2elin_core.timeseries.TimeSeriesResult` gained
  a `results: dict[str, PowerFlowResult]` field (the full per-snapshot
  result, not just the bus table `bus_tables` already held — a purely
  additive change, existing fields/tests untouched) so the API layer can
  pull the same six other-result tables per snapshot. The load-scaling
  sweep's network *topology* doesn't change across snapshots (only P/Q
  values do), so the diagram layout is fetched once and just recolored
  per snapshot's own bus results, not refetched three times.
- **Modal Analysis toolbox.** Covered above (`g2elin_core.modal.toolbox`)
  — the web UI side is a "View" dropdown exposing all of it: eigenvalue
  map, participation heatmap, single-mode participation bar chart,
  sensitivity heatmap (a 2D heatmap here, not the original's 3D bar plot
  — same information, more readable; capped to state counts ≤30, since an
  87×87 HTML table isn't a heatmap anyone can read — the API's top-8 list
  is the useful output past that size), mode shape (polar plot), free
  motion response and step response (both reuse the EMT tab's line-chart
  renderer).
- **Eigenvalue map: zoomable/pannable, full grid; mode table: top-3
  states; network diagrams: click for parameters.** The eigenvalue map now
  has a full x-y gridline set (`eigenGridLines()`, not just the origin
  crosshair) with tick labels inverse-transformed back to physical
  real-part/Hz units, plus mouse-wheel zoom and drag-pan via a small
  reusable `attachSvgZoomPan()` helper (SVG `viewBox` manipulation, no
  external library) and a "Reset zoom" button. The mode summary table now
  lists each mode's top **three** participating states, not just one — the
  backend already computed `state2`/`state3`/`part2_pct`/`part3_pct` (this
  was a frontend-only gap). Every bus/line/transformer in a network
  diagram (Network, Power Flow, and Time Series tabs, which all share
  `renderNetworkDiagram()`) is now clickable, not just hoverable, opening
  a detail panel with that element's actual parameters — a line's `r_pu`/
  `x_pu`/`b_pu`/`length_km`, a transformer's `r_pu`/`x_pu`/`sn_mva`, and
  for a DER-hosting bus its dispatch fields plus (for SM/GFM/GFL) the
  *derived* control-parameter dict from `sm_params()`/`gfm_params()`/
  `gfl_params()` — the actual values used to build that unit's dynamic
  model, not a re-statement of the static schema. `network.topology`'s
  `BusNode`/`TopologyEdge` and the `/topology` response schema were
  extended to carry this; it needs no power flow (those three functions
  are pure functions of static network data), so the endpoint stays cheap.
- **EMT tab: full state/input/output exposure.** The perturb-state picker
  and the three new "plot states/inputs/outputs" pickers are genuine
  `<select multiple>` elements over the *complete* list from `/states`
  (which now also returns `input_names`/`output_names`), not a
  substring text filter — the explicit ask this round. Recovering
  input/output *trajectories* needed new machinery:
  `NonlinearNetworkModel.recover_inputs_and_outputs()` re-solves the
  coupled algebraic system once per already-integrated sample (`solve_ivp`
  only returns state trajectories; the `z`/`u` it solves for internally
  live at the integrator's own adaptive step points, not the returned
  sample times) to recover a consistent `(z, u)` at each sample, then
  evaluates each block's own named outputs there. Costs about as much
  again as the integration itself, so it's opt-in (empty `plot_inputs`/
  `plot_outputs` skips it) rather than always paid for. States/inputs/
  outputs render as up to three *separate* charts, not one shared axis —
  they're physically different kinds of quantity (per-unit dynamics,
  fixed exogenous references, named physical outputs like P/Q/V), so one
  y-axis across all of them would risk exactly the "two measures of
  different scale" dual-axis problem the dataviz skill rules out.

  **Follow-up: perturb an input, not just a state; control the output
  timestep.** What can be perturbed was states-only at first; a state
  perturbation is an initial-condition offset (`x0[idx] += offset`) — the
  system starts away from equilibrium and (maybe) settles back. The more
  standard EMT disturbance test is a *step in an exogenous reference*
  (e.g. "what happens if `P_ref` steps up by 0.1 pu") — the system starts
  *at* equilibrium and the equilibrium itself moves. `timedomain.simulate()`'s
  `u_exo_fn` parameter already existed for exactly this ("lets a caller
  drive a disturbance... by returning a modified exogenous-input vector at
  each t") but nothing had ever called it with a non-default function —
  `run_emt` now builds one when `perturb_kind="input"`, and threads the
  *same* closure into `recover_inputs_and_outputs()` too (using the
  unperturbed default there instead would silently contradict the
  trajectory that was actually integrated). The web UI's perturb-target
  picker is now one merged `<select>` (states and inputs grouped via
  `<optgroup>`, kind read off the chosen option's `data-kind`) rather than
  a states-only dropdown; the "Offset" field is now labeled "Amplitude"
  since that's what it always was.

  The **output sampling timestep** (`dt`) is now a request parameter too
  — previously always exactly 200 samples over `[0, t_final]` (so `dt =
  t_final/199`, e.g. ~5.03 ms at the default `t_final=1.0`). This is the
  spacing of the *returned/plotted* samples, not the ODE solver's own
  adaptive internal step size (`Radau`, `rtol=1e-4`, `first_step=1e-8`,
  unexposed and untouched — that's an accuracy control, tuned once for
  this problem, not something to casually change per-request) — a smaller
  `dt` gives a finer plot, not a more accurate integration. Bounded
  server-side to 10–2000 samples (`dt=None` keeps the original 200-sample
  default), since each sample costs a full Newton solve when
  `plot_inputs`/`plot_outputs` is also set.

- **10 more WSCC/CIGRE-islanded preset variants (15 presets total).**
  `script_WSCC.m` and `CIGRE/script_CIGRE_Islanded.m` each define several
  DER-mix cases on one shared topology (only `Y_DER` differs between
  cases) — `network/presets.py` was refactored to match: a shared
  `_wscc9_topology()`/`_cigre_raw_topology()` per family, plus a thin
  per-variant function transcribing just its own DER rows
  (`_WsccDer`/`_CigreDer` `NamedTuple`s). Two of the eight WSCC variants
  (`wscc9_1gfm_2gfl`, `wscc9_3gfm`) have a **GFM as the slack** in the
  MATLAB source — power flow works (pandapower's slack handling doesn't
  care about unit type), but `interconnect/network_assembly.py` only wires
  up an SM/IB slack so far, so modal/EMT raise a documented
  `NotImplementedError` (HTTP 501) for those two, not a silent wrong
  answer. `g2elin_api/presets.py`'s `PRESETS` dict is now built from a
  `(id, name, description, build)` table rather than 5 (now 15)
  hand-written entries; the web UI's preset `<select>` groups them by
  family via `<optgroup>`.
- **Network editor and drag-and-drop builder.** `Network`
  (`g2elin_core/network/schema.py`) is already a pydantic model, so it
  round-trips to/from JSON for free — an edited or from-scratch-built
  network is just another `Network`, validated the same way a preset is.
  New `/api/network/*` endpoints (`network_routes.py`) mirror every
  preset-id-keyed endpoint but take a `Network` JSON body instead;
  `analysis.py` holds the actual analysis logic both sets of endpoints
  call into, so nothing is duplicated between them. Nothing is persisted
  server-side (confirmed with the user: ephemeral, client-held edits
  only) — the web UI clones a preset via `GET .../network`, holds it in
  `state.customNetwork`, and two dispatcher functions
  (`runAnalysis()`/`fetchTopologyActive()`) transparently redirect every
  existing tab's own run button to `/api/network/*` while it's set. A
  slide-down `#network-editor` panel (not a new tab, since edits must
  affect every tab identically) gives a table-per-element-kind editor
  (inline-editable buses/lines/transformers/loads/DER units, add/delete,
  a debounced live preview) plus an "Open canvas builder" mode: a palette
  of draggable element chips (Bus/Load/SM/GFM/GFL/IB), an explicit
  Select/Wire-Line/Wire-Transformer mode toggle for connecting buses
  (disambiguates "drag a bus to move it" from "drag a bus to wire it" on
  the same circle), and a click-to-edit parameter panel — all pointer-event
  driven (`setPointerCapture`, matching the eigenvalue map's own
  `attachSvgZoomPan()` precedent), not native HTML5 drag-and-drop. A
  dragged bus updates only the specific SVG attributes of what's visually
  connected to it rather than re-rendering the canvas on every
  `pointermove` — a full re-render would recreate the dragged circle
  itself and silently drop its pointer capture mid-gesture. Not yet
  supported (said explicitly, not silently dropped): multi-select,
  undo/redo, copy/paste, edge auto-routing, cross-reload persistence.

- **`network.validate_network()`: one comprehensive structural pre-flight
  check.** A hand-built network can violate invariants every preset
  satisfies by construction (`network_form.m`'s "each DER sits behind its
  own transformer" convention, primarily) in several independent ways --
  originally each surfaced one at a time as a crash (a bare `KeyError`/
  `IndexError`, deep inside `operating_point.py` or
  `interconnect/network_assembly.py`) the moment a user happened to hit it,
  requiring a fix-one/hit-the-next cycle to find them all.
  `network/validation.py`'s `validate_network(network) -> list[NetworkIssue]`
  instead checks everything in one pass -- duplicate bus/DER ids,
  self-loop lines/transformers, a DER not behind its own transformer, that
  transformer landing on another DER's bus, a Line/Load on a DER's own
  bus, graph connectivity, and an unsupported slack unit type -- and
  returns *every* issue found, each tagged `"error"` (blocks the affected
  capability) or `"warning"` (e.g. a GFM/GFL slack: power flow still
  works, modal/EMT don't support it yet). `compute_operating_point()`
  calls this itself and raises one combined `ValueError` listing every
  error (-> HTTP 422, not a 500) regardless of whether a caller checks
  proactively; a new `POST /api/network/validate` endpoint also exposes it
  directly, and the web UI's network editor calls it on the same debounce
  as its live preview, showing every problem in one panel instead of
  discovering them one crash at a time.

- **Full physics documentation, and per-tab pedagogy in the web UI.** The
  Sphinx docs previously described every module's *methodology* in prose
  but never wrote out an actual equation — added `sphinx.ext.mathjax` +
  MyST's `dollarmath`/`amsmath` (`docs/sphinx/conf.py`) and, module by
  module, transcribed every component's real nonlinear DAE (SM, GFM,
  GFL, IB, line, node, load — all 7, straight from each `*_dae()`
  function's own sympy expressions, not paraphrased) into LaTeX in
  `modules/components.md`; the per-unit base formulas, default
  SM/AVR/PSS parameter tables, and every GFM/GFL loop's rise-time
  pole-placement formula in `modules/operating_point.md`; the AC
  power-flow bus-injection equations and per-unit-to-physical conversion
  in `modules/powerflow.md`; the closed-loop elimination algebra and
  wiring-rule KCL/KVL statements in `modules/interconnect.md`; the
  eigenvalue/damping/participation-factor and sensitivity/mode-shape/
  free-response/step-response formulas in `modules/pipeline.md`; the
  coupled-Newton DAE residual and its analytic Jacobian in
  `modules/timedomain.md`. Every equation was
  transcribed directly from the source it documents (Sphinx builds
  clean, `-W --keep-going`, 0 warnings) rather than re-derived from
  memory. Separately, `web/index.html` gained a short pedagogical
  `.panel-intro` card at the top of every tab (what it computes, the
  method in one line, and a concrete "useful for" case) plus a small
  hoverable/tappable "(i)" bubble next to each tab label in the nav bar
  itself with a one-sentence version of the same — both new,
  self-contained CSS/JS additions (`.info-icon`/`.tooltip`,
  `.panel-intro`), no new dependencies.
- **A single worked example, run end to end.** New
  `notebooks/cigre_1sm_1gfm_1gfl_walkthrough.ipynb` follows one preset
  (`cigre_islanded_1sm_1gfm_1gfl()` — the smallest network in this
  codebase with one of each major DER type, fully islanded) through
  every module in order, interpreting the *actual* computed results, not
  just calling each function. Executed, not hand-typed (`jupyter
  nbconvert --execute`, 0 cell errors) — which caught a real mistake
  before it reached the docs: an earlier draft picked "the least-damped
  mode" mechanically for its free-response/EMT-comparison sections, and
  that mode turned out to be the network's numerically-near-zero
  reference-angle eigenvalue (the "theta problem" — an absolute angle has
  no restoring force, so every closed-loop network here has exactly one
  structurally-zero eigenvalue), not a real dynamic mode; comparing a nonlinear
  trajectory against a frozen baseline for that kind of state diverges
  by construction and has nothing to do with instability. Fixed by
  searching explicitly for the genuine SM electromechanical swing mode
  (`dw_r`-dominated, found at ~2.2 Hz / ~22% damping) instead — the
  corrected notebook's linear-vs-nonlinear comparison for that mode now
  tracks to within ~5% of the initial perturbation over a full damped
  swing. A companion `docs/sphinx/cigre_walkthrough.md` page narrates
  the same walkthrough with the real results and five of the notebook's
  own exported plots embedded (`docs/sphinx/_static/cigre_walkthrough/`),
  linked from `index.md`'s toctree and its own architecture overview.
- **The notebook itself, rendered directly in the docs.** Added
  `myst-nb` (new `docs` extra in `pyproject.toml`) alongside
  `myst-parser` — it supersedes it as the registered Markdown/`.ipynb`
  parser (`source_suffix` now maps both to `"myst-nb"`; `conf.py`'s
  `nb_execution_mode = "off"` renders the notebook's own already-executed,
  already-verified outputs instead of re-running a multi-minute EMT
  simulation on every docs build). `notebooks/` stays the single
  git-tracked source of truth (same as `tour.ipynb`/`random_network.ipynb`)
  — `tools/build_docs.py` stages a copy into a new, gitignored
  `docs/sphinx/_notebooks/` immediately before invoking `sphinx-build`,
  the same way `_build/` itself is a disposable artifact, not a second
  copy to keep in sync by hand. One real fix needed along the way: a
  markdown cell inside the notebook had used a Sphinx-only `` {doc}` ``
  cross-reference role, which is meaningless (and would render as
  literal broken syntax) if that same notebook is opened directly in
  plain Jupyter — replaced with plain text so the notebook stays fully
  portable/standalone, embedded-in-Sphinx or not. Linked from
  `index.md`'s toctree and from `cigre_walkthrough.md`'s own intro as
  `_notebooks/cigre_1sm_1gfm_1gfl_walkthrough`.
- **`notebooks/random_network.ipynb` embedded the same way**, added to
  `EMBEDDED_NOTEBOOKS` in `tools/build_docs.py` and `index.md`'s
  toctree. This one already had no Sphinx-only syntax in its own
  markdown (checked before embedding, having just found that exact bug
  in the CIGRE notebook), and re-executed clean (0 cell errors) to
  confirm it's still current against today's codebase before staging it
  in. Cross-linked from `modules/network.md`'s own "Validation" section
  — the natural home, since the notebook is `validate_network()`
  exercised hands-on. Fixed that section's own stale "Five networks"
  opening line while there (the preset-expansion work earlier in this
  session brought the real count to 15 across 4 families, but never
  updated this one paragraph) — flagged to the user as a drive-by fix
  found via the file, not silently absorbed into this notebook's own scope.
- **`random_network()` retargeted to a specific, requested topology**:
  radial (its grid *is* a random spanning tree now, no chords added on
  top — was meshed with a few redundant lines before), at least 15 buses
  (was 3-10), and exactly 2 SM + 2 GFM + 2 GFL (was 1 forced-slack SM +
  3 more of a randomly-chosen type) each dispatched to a *distinct* grid
  bus via `rng.choice(..., replace=False)` (dispatch buses could
  previously collide). Two real bugs surfaced and were fixed while
  re-verifying results at the new scale, not assumed to still hold from
  the smaller version: (1) with 15+ load buses now splitting one
  Dirichlet draw, a share could round to an exact 0.00 MVAr, making that
  load's equivalent reactance exactly zero — a division by zero deep in
  the linearizer, surfacing as a cryptic "cannot convert complex to
  float" (sympy's `zoo` for `1/0`); fixed with a strictly-positive floor
  on both P and Q, not just finer rounding. (2) The line-impedance range
  was tuned for a handful of buses with redundant paths (WSCC/
  transmission scale); on a purely radial chain of 15+ buses with no
  redundant path to share current, that impedance accumulated hop after
  hop and dropped the batch loop's power-flow convergence to 18/30 —
  rescaled to CIGRE's own per-line order of magnitude (a real MV
  distribution feeder), confirmed 30/30 across the same 30 seeds
  afterward. Also fixed the notebook's own "stable" checks (a leftover
  `< 0` reads `False` for *every* network in this codebase, not a
  genuine instability — the same structural near-zero "theta problem"
  eigenvalue the CIGRE walkthrough bullet above already documents) to a
  `< 1e-6` tolerance instead, now correctly reporting 0 unstable modes.
- **`numba` added as a real dependency, not left implicit.** pandapower
  optionally JIT-compiles its Newton-Raphson internals with numba;
  without it, every single power flow call in this project — every test,
  every notebook cell, every API request — silently fell back to a much
  slower pure-Python path and printed a warning saying so on first use
  (visible dozens of times over in this session's own notebook output).
  Installed it and added `numba>=0.60` to `pyproject.toml`'s core
  `dependencies` (not an extra — every consumer of `g2elin_core` benefits,
  the same way pandapower/scipy/sympy already aren't optional). No
  behavior change (28/28 power-flow tests still pass); confirmed via a
  direct timing check the warning is gone and the steady-state (JIT
  already warmed up) WSCC-9 solve takes ~0.17s.
- **Every `modal.toolbox` and `timedomain.emt.simulate()` option now
  exercised in both the CIGRE and random-network notebooks**, plus a
  `RUN_EMT_SECTION` toggle (default `False`, matching `tour.ipynb`'s own
  established convention) gating the slow nonlinear cells in each.
  `eigenvalue_sensitivity()` was the one `modal.toolbox` function neither
  notebook had used before (`mode_shape`/`free_response`/`step_response`
  already were) — added to both. `simulate()`'s full parameter surface
  (`t_span`, `x0`, `u_exo_fn`, `t_eval`, `method`, `rtol`, `atol`,
  `first_step`) is now documented in a table and exercised via two
  demonstrations per notebook: perturbing a **state** (the existing
  linear-vs-nonlinear comparison) and perturbing an **input** via
  `u_exo_fn` (a permanent reference step, held from $t=0$ — the same
  `perturb_kind="input"` path `g2elin_api`'s own EMT endpoint already
  supported but neither notebook had demonstrated). Every new *gated*
  code path (state perturbation, input perturbation) was verified
  standalone with the toggle forced `True` before shipping with it
  defaulted `False` — a toggle guard only proves the *unguarded* path
  still runs; it says nothing about whether the code inside the guard is
  even correct. That check caught a real problem: a leftover dead-code
  artifact (`... if False else None`) in the CIGRE notebook's own
  step-response cell, predating this change, cleaned up while there.
  Also fixed a `myst-nb` structural-linter warning (H1 straight to H3,
  "Non-consecutive header level increase") both new toggle-explanation
  cells introduced now that these notebooks render through Sphinx, not
  just plain Jupyter — `### Run-time toggles` -> `## Run-time toggles` in
  both.
- **`random_network.ipynb`'s Section 6a split into a live-traced popup
  window + a separate final static comparison plot** (still gated by
  `RUN_EMT_SECTION`, unchanged elsewhere). Bypasses `simulate()`'s
  one-shot `solve_ivp()` call in favor of driving `scipy.integrate.Radau`'s
  own low-level OOP stepper directly (`.step()` in a loop), replicating
  `simulate()`'s own warm-started coupled-algebraic-solve pattern by hand
  since a one-shot call can't yield control back between steps. A real
  limitation surfaced verifying this: **this repository's own tools have
  no way to observe an external GUI window** — my first attempt to verify
  the live cell via this project's usual `nbconvert --execute` pipeline
  produced no output for 90+ seconds, which I (wrongly) read as a hang
  and killed; the user, watching their own screen, confirmed it was
  genuinely tracing live the whole time. A separate GUI-free check
  (same stepping logic, no window) settled the "correct or not"
  question independently either way: it matches `simulate()`'s own
  reference trajectory to within 4e-5 pu — but also took 341s of pure
  computation (2051 steps at `max_step=0.01`) before any redraw
  overhead, confirming the long wait was real work, not a stuck
  process. Left `max_step=0.01` as originally written once the user
  confirmed the pacing was already what they wanted, rather than
  "fixing" a slowdown that wasn't a complaint. Net effect: this specific
  cell's own correctness can no longer be fully confirmed by this
  project's own automated notebook-verification discipline (every other
  cell in every notebook here is) — only by a human actually running it
  and watching, which is inherent to what a live GUI window is, not a
  gap to close.

Not yet started: legacy `Y_*` importer, other GFM controller variants,
a GFM slack (see the two WSCC variants above), and (once Node.js is
available) migrating the frontend to the planned React/TypeScript/Vite
stack.

## Setup

```powershell
py -m venv .venv
.venv\Scripts\Activate.ps1
pip install -e ".[dev]"          # core only
pip install -e ".[dev,api]"      # core + the FastAPI backend
pip install -e ".[dev,api,docs]" # + Sphinx docs (see "Sphinx documentation" under P7 above)
pip install -e ".[notebook]"     # + Jupyter/matplotlib, for notebooks/tour.ipynb below
```

## Run the tests

```powershell
pytest -v
```

`tests/golden/` holds regression tests meant to compare against the
original MATLAB tool's output; see `tests/golden/README.md` for how to
generate the reference values that test currently lacks (no MATLAB
license in this environment).

## Build the documentation

```powershell
python tools/build_docs.py
```

Builds `docs/sphinx/` to `docs/sphinx/_build/html` (0 warnings on a clean
build). Run the API (`uvicorn g2elin_api.main:app`) afterward and it's
served at `/manual`, and from the **Documentation** tab in `web/index.html`.

## Run the API + web UI

```powershell
uvicorn g2elin_api.main:app --reload
```

Then open `http://127.0.0.1:8000/` — FastAPI serves both the JSON API and
the static frontend from one process. The first modal-analysis request per
preset takes a few seconds (symbolic model derivation, cached in-process
after that); everything else is fast.

## Notebook: a tour of every module

`notebooks/tour.ipynb` exercises every `g2elin_core` module directly (no
web UI, no API layer) in the order data actually flows through the tool —
`network` -> `pu_base` -> `powerflow` -> `network.topology` ->
`operating_point` -> `components` -> `interconnect` -> `pipeline` ->
`modal` (`analysis` + `toolbox`) -> `timedomain` (EMT) -> `timeseries` —
plotting each result with matplotlib. Run it with:

```powershell
pip install -e ".[notebook]"
python -m ipykernel install --user --name g2elin --display-name "G2ELin (.venv)"
jupyter lab notebooks/tour.ipynb
```

Already executed top-to-bottom (outputs committed) so it's readable
without running it first. The EMT section is the slowest cell (full
nonlinear coupled-Newton time integration, not just linear algebra) and
defaults to **off** via the `RUN_EMT_SECTION` toggle in a cell near the
top, so "Restart Kernel and Run All" finishes quickly; flip it to `True`
to exercise it.

A "9b" section traces a root locus of a GFM's inner current-loop tuning
(`operating_point.gfm_params(..., tr_cl=...)`, transcribed from
`script_generic.m`'s pole-placement formula) against
`cigre_islanded_1sm_2gfm_1gfl()`: first the one current-loop pole in
isolation, then two more views — every mode's trajectory in the complex
plane as `tr_cl` sweeps (colored by `tr_cl`, arrow on the traced pole,
zoomed inset for the near-origin cluster; a greedy nearest-neighbor
eigenvalue trace, since `scipy.linalg.eig` doesn't return eigenvalues in a
parameter-continuous order — this surfaces other modes that also cross
into instability, not just the traced one), and a dual-axis plot of that
pole's real part and damping ratio against `tr_cl` with the interpolated
stability-boundary crossing marked.

## Layout

```
src/g2elin_core/
  network/          typed schema + preset cases + topology.py (2D layout) + validation.py (structural pre-flight checks)
  powerflow/        pandapower adapter
  components/       symbolic linear + nonlinear component models (SM, GFM, GFL, IB, line, node, load)
  pu_base.py        AC/DC per-unit base-value calculator
  operating_point.py  solved power flow -> per-component linearization points
  interconnect/     shared topology (Block/Wiring/F/G) + linear closed-loop assembly
  modal/            eigenvalues + participation factors + toolbox.py (sensitivity, mode shape, free/step response)
  timeseries/       time-series load flow (dispatch/load snapshots)
  timedomain/       EMT (nonlinear) time-domain integration
  pipeline.py       end-to-end: power flow -> linearized closed-loop model
src/g2elin_api/     FastAPI backend over g2elin_core (optional `[api]` extra)
  presets.py          the fixed, named preset catalog (15 entries)
  analysis.py          analysis logic shared by the preset-id and Network-JSON endpoints
  network_routes.py    /api/network/* -- run any analysis on an arbitrary Network JSON body
  main.py               preset-id-keyed endpoints (thin wrappers over analysis.py) + app wiring
web/                static frontend (plain JS — see the P7 README note on why)
notebooks/          tour.ipynb -- every module exercised directly (optional `[notebook]` extra)
tests/
  golden/                             regression fixtures against the MATLAB original
  test_nonlinear_cross_validation.py  linear vs. nonlinear cross-checks (no MATLAB needed)
  test_timeseries.py                  time-series load flow sanity checks
  test_emt_simulation.py              EMT integration smoke tests
  test_modal_toolbox.py               modal toolbox tests (cross-checked against scipy directly)
  test_smib_presets.py                SMIB preset + IB-slack wiring tests
  test_wscc_cigre_variants.py         the 10 additional WSCC/CIGRE DER-mix preset variants
  test_api.py                         API smoke tests (needs the `api` extra)
  test_network_endpoints.py           /api/network/* smoke tests, incl. the pydantic-422 path
```
