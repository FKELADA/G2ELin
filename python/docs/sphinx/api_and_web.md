# `g2elin_api` and the web UI

`g2elin_api` is a thin FastAPI layer wrapping the real `g2elin_core`
functions (no mocking). `web/` is a build-free, vanilla-JS single-page app
that calls it (`index.html` + `app.css` + one script per page under
`web/js/`). This page (the one you're reading, if you got here through the
**Documentation** page) is itself served as part of that same app — see
[Mount layout](#mount-layout) below.

## Endpoints → pages

The UI loads a preset's `Network` once (`GET /api/presets/{id}/network`)
and from then on calls only the `/api/network/*` endpoints with that
(possibly edited) network. The preset-id endpoints below are kept for
scripts and tests; both sets share the same `analysis.py` functions.

```{mermaid}
flowchart LR
    subgraph API["g2elin_api (FastAPI)"]
        EP1["GET /api/presets<br/>GET /api/presets/{id}/network"]
        EPT["POST /api/network/topology<br/>POST /api/network/validate"]
        EP2["POST .../powerflow<br/>(solver options)"]
        EP2B["POST .../powerflow/batch"]
        EP3["POST .../modal"]
        EP3X["POST .../modal/sensitivity, mode_shape,<br/>free_response, step_response"]
        EP5["POST .../states"]
        EP6["POST .../emt<br/>(t_pre, linear_overlay)"]
        EP6L["POST .../emt/live<br/>(streamed NDJSON)"]
    end
    subgraph UI["web/ pages (left navigation)"]
        T0["Home"]
        TN["Network"]
        T1["Power Flow"]
        T2["Modal Analysis<br/>(7 sub-pages)"]
        T4["EMT Simulation"]
        T6["Documentation"]
    end
    EP1 --> T0
    EP1 --> TN
    EPT --> TN
    EP2 --> T1
    EP2B --> T1
    EP3 --> T2
    EP3X --> T2
    EP5 --> T4
    EP6 --> T4
    EP6L --> T4

    EPT -.-> TOPO["network.compute_topology_layout()<br/>network.validate_network()"]
    EP2 -.-> PF["powerflow.run_power_flow(**runpp options)"]
    EP2B -.-> TS["timeseries.apply_snapshot() per ramp step"]
    EP3 -.-> PIPE["pipeline.linearize_network()<br/>+ modal.analyze()"]
    EP3X -.-> TOOL["modal toolbox"]
    EP5 -.-> BNM["timedomain.build_nonlinear_network()"]
    EP6 -.-> EMT["timedomain.simulate()<br/>+ scipy.signal.lsim on (A, B, C, D) for the overlay"]
    EP6L -.-> EMTLIVE["timedomain.simulate_steps()"]
```

## `/api/network/*` — arbitrary networks, not just presets

Every endpoint above also has a counterpart under `/api/network/*`
(`network_routes.py`) taking a full `Network` JSON body instead of a
`preset_id` path param — `POST /api/network/powerflow`, `/modal`,
`/modal/sensitivity`, `/modal/mode_shape`, `/modal/free_response`,
`/modal/step_response`, `/timeseries`, `/emt`, `/emt/live`, plus `/topology` and
`/states` as POST (their preset-side counterparts are GET, but a `Network`
body can't ride a GET request — an intentional divergence). There's no
preset-side equivalent of `POST /api/network/validate` (a preset is always
already known-valid) — it runs `network.validate_network()` (see
{doc}`network <modules/network>`) and returns every structural issue
found, not just the first, without needing power flow to have run first;
the web UI's network editor calls it on a debounce as the user edits.
Since `Network`
(`g2elin_core.network.schema`) is already a pydantic model, it round-trips
to/from JSON for free — a client-edited or from-scratch-built network is
just another `Network`, validated by pydantic before any handler runs (bus
references resolve, exactly one slack, etc.), same as any preset.

Both sets of endpoints delegate to the same functions in `analysis.py`
(each takes a `Network`, returns the exact response model its endpoint
returns) so the actual analysis logic exists in exactly one place —
`main.py`'s preset handlers are just "resolve `preset_id` -> `Network`
(404 on unknown id), then call into `analysis.py`," and
`network_routes.py`'s handlers skip straight to the same call with the
POSTed `Network`. Nothing here is persisted server-side: the web UI's
network editor holds the edited `Network` client-side and resends it whole
on every call (see "The network editor," below) — this is a deliberate
scope choice, not a missing feature; see the module docstring in
`g2elin_api/presets.py`.

`GET /api/presets/{id}/network` returns a preset's own `Network`
definition verbatim — how the web UI clones a preset into an editable copy.

Response models (`g2elin_api/schemas.py`) are deliberately separate from
`g2elin_core`'s own result types — `PowerFlowResult` wraps a live
pandapower net, `ModalAnalysisResult` holds numpy arrays, neither is
JSON-serializable — so the wire format is a considered choice, not
whatever `dict()` happened to produce.
`PowerFlowResponse`'s `lines`/`transformers`/`loads`/`generators`/
`static_generators`/`external_grid` fields are loosely typed (`list[dict]`,
not a row model per table) since each pandapower result table's column set
genuinely differs (13 columns for lines, 2 for loads) and this is an
inspection/dropdown feature, not something else in this codebase consumes
programmatically.

`/states` exists because the nonlinear model's state/input/output names
are entirely preset-dependent — which DER id is the slack, which unit
types exist, whether a `dw_r_*` (synchronous-machine speed deviation)
state exists at all — so the EMT tab's pickers are populated from
this endpoint at preset-selection time rather than hard-coded, the same
problem `find_state_index()`'s substring matching solves (which the EMT
tab's *perturb*-state selection reuses server-side). `/emt`'s
`plot_states`/`plot_inputs`/`plot_outputs`, by contrast, take **exact**
names from that same list (not substrings) — real multi-select pickers
over the full list, not a text filter — with an unknown name coming back
as a 422 listing what *is* known, not a silent empty result.

**`/emt/live`** is the same request body, streamed instead of returned in
one piece — newline-delimited JSON, one line per accepted solver step
(`{"t": ..., "states": {...}, "inputs": {...}, "outputs": {...}}`), a
final `{"done": true, ...}` line on success, or a single `{"error": "..."}`
line if the coupled Newton solve fails partway through (a stream that's
already started can't turn into an HTTP error status any more, unlike the
one-shot endpoint's clean 422 on the same failure). Backed by
`timedomain.simulate_steps()` — the manually-stepped sibling of
`simulate()`, driving `scipy`'s low-level OOP stepper class directly
instead of the one-shot `solve_ivp()` wrapper, since only the low-level
interface can yield control back between steps. All validation
(`analysis.prepare_emt_live()`) runs synchronously in the route handler
*before* `StreamingResponse` is constructed, so a bad request still 422s
normally rather than producing a broken stream. The web UI's own EMT tab
gates this behind an opt-in **"Trace live"** checkbox next to the regular
"Run EMT simulation" button (default off) rather than a second button —
tracing step by step is a real speed trade-off, not just a display
preference: a trajectory the one-shot endpoint answers in a couple of
seconds can take several minutes traced live (same integration, same
tolerance; every step now separately pays a JSON-serialize + HTTP-write +
DOM update it wouldn't otherwise pay), so it needs to be something a user
opts into, not something that silently replaces the fast path.

## Mount layout

```{mermaid}
flowchart TD
    ROOT["FastAPI app"] --> API["/api/*<br/>JSON endpoints (defined in main.py)"]
    ROOT --> MANUAL["/manual<br/>StaticFiles: docs/sphinx/_build/html<br/>(this Sphinx site — build it first, see below)"]
    ROOT --> WEB["/ (catch-all, mounted last)<br/>StaticFiles: web/, html=True"]
```

Mount order matters: FastAPI/Starlette match mounts in registration order,
so `/api/*` routes and the `/manual` static mount are registered **before**
the catch-all `/` mount for `web/` — otherwise the catch-all would shadow
them. `/manual` (not `/docs`) was chosen specifically to not collide with
FastAPI's own built-in Swagger UI, which lives at `/docs` by default and is
left alone.

The Sphinx site has to be built before `/manual` serves anything real —
`app.py` checks whether `docs/sphinx/_build/html` exists and only mounts it
if so, so a fresh checkout without a docs build doesn't crash the app, it
just doesn't have a Documentation tab content yet. Build it with:

```bash
pip install -e ".[docs]"
python tools/build_docs.py
```

## The web UI itself

A fixed left navigation panel with six pages, routed by URL hash
(`#/home`, `#/docs/<page>`, `#/network`, `#/powerflow`,
`#/modal/<view>`, `#/emt`). Every page works on the same client-held
`Network` (`state.network` in `web/js/core.js`); editing it bumps a version
counter that invalidates every cached result, and nothing is ever
persisted server-side. Always light-themed.

- **Home** — what the tool is, its features and workflow, the component
  models, and the preset list (click one to load it).
- **Documentation** — this Sphinx site in an `<iframe>`. Its table of
  contents is read from the built `index.html` and listed as sub-items
  of *Documentation* in the app's own navigation (with the active page's
  sections underneath); Furo's own sidebars are hidden inside the frame so
  the pages read as part of the app. In-frame navigation is mirrored into
  the app URL.
- **Network** — preset picker (the network is plotted immediately), a
  drag-and-drop canvas (palette of Bus / Load / SM / GFM / GFL / IB; a unit
  dropped on or near a bus gets its own terminal bus and transformer;
  *Draw line* / *Draw transformer* modes), and an inspector that opens on
  the right when an element is clicked, listing **every** editable
  `Network` field for it (bus, its loads and its unit; line; transformer;
  network-level settings) plus the unit's derived control parameters
  (read-only). Live validation (`/validate`), a bulk table view, and JSON
  import/export (network + diagram positions).
  A unit's control & electrical parameters are editable: changed values are
  stored as overrides in `DerUnit.params` (only what differs from the
  default) and applied on top of `sm_params()`/`gfm_params()`/`gfl_params()`
  by `operating_point.compute_operating_point`, so modal analysis and EMT
  both use them; an unknown parameter name is a validation error. Defaults
  come from `POST /api/units/defaults` (unit type + base values only), so
  they're available before the unit is connected. For GFM/GFL units a
  **loop tuner** converts response time t_r and damping ζ to Kp/Ki and back
  with the same pole-placement formulas as the defaults (ωn = 3/(ζ·t_r)).
  Changing a plant value a loop was tuned against (Rf, Lf, Cf, Cdc, Gdc)
  does not retune its gains; the tuner flags the affected loops and can
  re-apply the default t_r/ζ to the new plant.
  Every per-unit value is shown next to its SI equivalent, and either can
  be edited (per unit stays the stored value; SI is derived from it with
  the current base values, so changing a bus voltage or the base power
  changes the SI figures, not the per-unit ones). Lines show R/X in Ω and
  B in µS, both for the whole line and per km; transformers in Ω referred
  to the HV side on their rating; a unit's R/L/C in Ω/mH/µF on its own
  base (DC link in mF and µS on the DC base, `pu_base.py`). A line's
  R, X and B are whole-line values, so its length has no effect on the
  results by itself: editing the length rescales them at constant per-km
  values (the root-locus sweep of a length does the same).
  A unit's Rt/Lt are its own transformer -- the element the power flow
  uses -- converted from the transformer's rating to the system base
  (`operating_point.unit_transformer_rx`), and the editor links the two:
  editing Rt/Lt edits the transformer and vice versa. The network setting
  *MATLAB-compatible transformers* (`Network.units_use_first_transformer`)
  restores the MATLAB tool's convention of every unit using the first
  transformer, to reproduce its results.
  **Breakers** (`g2elin_core.network.breakers`): every line and transformer
  has one at each end (`from_closed`/`to_closed`, `hv_closed`/`lv_closed`),
  every load and unit one at its bus (`closed`); all closed by default. They
  are drawn on the diagram (filled square = closed, hollow red = open) and
  toggle on click — on the Network and Power Flow pages — or from the
  inspector/table checkboxes. A line or transformer with either breaker
  open, and a load or unit with its own open, is out of service; so is
  anything cut off from the slack's bus (a de-energized island). The slack
  unit's breaker can't be opened. Out-of-service elements are drawn dashed /
  faded, validation lists them, and every analysis leaves them out: the
  power flow keeps them as pandapower `in_service=False` elements (so result
  tables keep their numbering; de-energized buses read null), and modal
  analysis, root locus and EMT are built from
  `breakers.energized_network()` — the network without them — whose model
  names keep the full network's numbers (line #3 is always `Ln_4`, the
  second SM always `SM_2`; `breakers.BlockLabels`). A unit that is out of
  service takes its own transformer and terminal bus with it.
  **Any breaker may be opened, the slack unit's included.** Open breakers
  split the network into islands, and each island that still holds a unit
  able to set a voltage and a frequency — a synchronous machine, a
  grid-forming converter or an infinite bus — is solved against its own
  reference (its own `ext_grid`), as every load-flow tool requires: the
  designated slack keeps the role in its island, otherwise the largest grid
  former takes it. An island left with only grid-following converters and
  loads is blacked out, since a grid-following converter can only follow a
  voltage, never start one (and anti-islanding protection would trip it).
  Only a network with no grid former left at all is an error. In the dynamic
  models this works because the reference frame is a block of its own
  ({doc}`frame.py <modules/components>`), one per island, so no unit is
  load-bearing for the model — EMT can trip the slack and watch the rest
  island and drift to its own frequency.
- **Power Flow** — solver (`nr`, `iwamoto_nr`, `fdbx`, `fdxb`, `gs`,
  `bfsw`), max iterations, tolerance and initialisation, passed to
  `pandapower.runpp`. The network diagram is shown before any run;
  results colour it as a heatmap (buses on a blue→red diverging scale
  centred on 1 pu or 0°, lines/transformers light→dark by |P|, |Q|,
  losses or current) and every element's own results appear on hover.
  **Batch power flow** opens sliders (0–200 %, default 100 %) for load
  P, load Q and each unit type's P/Q setpoints; it solves a ramp of
  independent power flows from the base point to those targets, with a
  snapshot scrubber (and *Play*) over the diagram and overview charts.
  The **Show** selector sits directly above the results table; clicking
  an element on the diagram jumps to its row.
Every figure and table can be taken out of the page (`web/js/exports.js`):
each one gets a small toolbar in its corner with **PNG** and **SVG** for the
figure and **CSV** for the numbers behind it -- a chart's own series, a
result table's rows, or, on the root-locus page, the selected mode's
participation across the sweep. A `MutationObserver` decorates figures as
they are drawn, so no page has to opt in. The exported SVG carries its
styles on the elements themselves (read off the live page with
`getComputedStyle`), since the app's rules are descendant selectors that
stop matching once the figure stands on its own; the PNG is that SVG drawn
onto a canvas at 2x.

In a chart's legend, clicking an entry **hides that signal** and clicking it
again brings it back; **double-clicking shows only it**, and double-clicking
again restores the rest. A signal hidden this way also leaves the hover
tooltip, and the PNG/SVG export follows what is on screen (the CSV always
holds every series). Entries drive every chart in their own card, so a scope
split over stacked subplots switches together.

Signals, everywhere they are picked, are organised by the **element** they
belong to: a name's trailing block (`..._{SM_2}`, `..._{Nd_4}`, `V_{bus4}`,
`P_from_{line3}`) is mapped back to the unit, bus, line or load it names, so
a bus's states and its measurements sit under the same entry (`core.js`'s
`elementCatalog()`). A perturbation is chosen as **element type → element →
signal**, and the "+ add signal" pickers ask the same way before listing
anything. Buses, lines and loads have two states each and are rarely what a
study is about, so a tick-box (on by default) keeps their *states* out of the
lists — their measurements are unaffected, and whatever is already plotted or
perturbed is never hidden away.

- **Modal Analysis** — one linearisation per network version (cached
  server-side too, see `analysis._cached`), seven sub-pages: eigenvalue
  map (+ mode table, click to select a mode), participation heatmap,
  single-mode participation, sensitivity heatmap, mode shape, free-motion
  response and step response. The last two hold *channels* — one
  perturbation (free motion) or one input step (step response) each —
  and every channel has any number of *subplots*, each with its own set of
  signals; crosshairs are synced across a channel's subplots.
  **Root locus** sweeps any one parameter (network, bus, line,
  transformer, load, unit setpoint or any `params.*` value) over a range
  and step: `POST /api/network/modal/sweep` re-solves the power flow and the
  eigenvalues at each value (at most 201 values), streaming one NDJSON line
  per value. Eigenvalues are reordered step to step by a minimum-cost
  assignment (`g2elin_api/sweep.py`), so each mode draws as one locus,
  coloured by the parameter value, on symlog or linear axes; a table ranks
  the modes the parameter moves most. A value that doesn't solve is
  reported for that step and the sweep continues.
  For GFM/GFL units the sweepable quantities also include each control
  loop's response time t_r and damping ζ (and the droop loop's inertia H
  and power-filter time constant), as `tune.<loop>.<quantity>` targets:
  at each value `g2elin_core.tuning` re-derives that loop's Kp/Ki with the
  same pole-placement formulas as the defaults, holding the loop's other
  quantity at its current value.
  Several parameters can be swept together (`extra` in the request):
  each moves from its own start to its stop in lockstep with the main
  parameter, so every step is their combined effect.
  Selecting a mode (on the plot or in the table) opens its **participation
  beside the locus, moving**: the same reading as the single-mode
  participation page, played through the sweep, so a mode's composition can
  be watched changing with the parameter (one machine handing a swing mode
  over to another, an inner loop taking over as a gain rises). The bars keep
  a fixed order — each state's largest participation over the whole sweep —
  so only their lengths move, a ring on the locus marks the value being
  shown, and hovering a point of that locus jumps the bars to it. The factors
  come from the sweep itself: every value's eigendecomposition already
  computes them, so each solved line carries each mode's top
  `sweep.PARTICIPATION_TOP` states (`participation: false` turns this off if
  the stream should stay small; about 14 kB per value for a 78-state model).
  **Video** records the sweep as it plays: one frame per solved value, the
  loci growing as the parameter moves, with a caption strip carrying the
  parameter's value and — when a mode is selected — that mode's frequency and
  damping, in red once it turns negative. The frames are the page's own plot,
  serialized exactly as the PNG export is and painted onto a canvas that
  `MediaRecorder` captures: MP4 where the browser can write it (Chrome,
  Edge), WebM otherwise. The canvas is repainted while each value is held,
  since a canvas left untouched emits no frames and the video would come out
  a fraction of its intended length. The participation panel records the same
  way (**Video** beside its CSV): its bars are HTML rather than a figure, so
  those frames are drawn straight onto the canvas -- the state names, the
  bars and their values, over the parameter, frequency and damping of that
  value.
- **EMT Simulation** — state offset, input step or **network event** at
  T0 = 0, integrated
  with the nonlinear model; plots start 0.01 s before T0 (`t_pre`) so the
  initial point x0 is visible. *Overlay the equivalent linearised
  response* adds the linear model's response to the same disturbance
  (`linear_overlay`: `lsim` of the full `(A, B, C, D)` around the same
  operating point, in absolute units) as dotted lines in the same colours.
  The network diagram is folded into the page's network strip (closed by
  default, and breakers can be operated from it), and *Trace live* is on by
  default.
  The nonlinear model's initial point is not an exact equilibrium (see
  {doc}`modules/timedomain`), so the EMT traces include a small initial
  transient the linear response doesn't have — the difference
  EMT(disturbed) − EMT(undisturbed) is what matches the linear response
  (checked in `tests/test_analysis_options.py`). *Trace live* streams the
  solver steps as before.
  **Limits, and what happens at them.** A one-shot run is limited to
  `analysis.EMT_MAX_N_POINTS` = 10 000 samples (each costs a Newton solve
  when inputs, outputs or measurements are plotted), and a live stream to
  `EMT_LIVE_MAX_STEPS` = 50 000 solver steps -- a bound on what one stream
  may send, not a physical limit: how many steps an adaptive stiff solver
  takes depends on the network, and a stiff one (CIGRE) takes far more than a
  small one. Reaching the live limit is a normal end of the stream, not an
  error: the final line carries a `truncated` message saying where it
  stopped, and the page keeps every plot drawn so far with that message above
  it. The same holds when the solver fails part-way through: the partial trace
  stays, with the reason, and only a run that fails before its first step is
  shown as a bare error.
  **Measurements** (`g2elin_core.timedomain.measurements`, requested with
  `plot_measurements`, listed by `.../states`) are computed from each
  element's own model variables at every sample: bus voltage magnitude
  (pu), angle (deg, equal to the power flow's at t = 0), frequency (Hz --
  instantaneous, from the analytic derivative of the voltage angle, and
  "measured", through a one-cycle first-order filter) and instantaneous
  3-phase voltages (pu); P and Q (pu on the network base) at both ends of
  each line (π-line terminal flows, charging included), consumed by each
  load, and injected by each unit at its grid-side bus (which is also its
  transformer's HV-side flow); and each unit's own frequency (rotor, droop
  or PLL). At the operating point they reproduce the power flow
  (`tests/test_measurements.py`). A unit's LV terminal bus is not a node of
  the dynamic model, so it has no bus measurements. When 3-phase voltages
  are plotted the page picks a timestep of about 40 samples per cycle.
  **Network events** (`perturb_kind: "event"`, `event: {...}`;
  `g2elin_core.timedomain.events`): a *breaker opening* on a line, a unit
  transformer (trips its unit), a load or a unit; a *load step* (ΔP/ΔQ in %
  of the load, as a constant impedance at the operating-point voltage); or
  a *phase jump* at a bus (rotates that bus's voltage phasor) or of the
  infinite-bus source (the grid-code phase-jump test). A breaker opening or
  load step integrates a *post-event model* built from the pre-event
  model's own component instances — same parameters and setpoints, no new
  power flow — rewired without the element or with the new load impedance,
  starting from the pre-event state block by block; whatever the opening
  islands keeps running on its own. Signals of removed elements read null
  after T0 (power flows read 0). The slack can't be tripped. Those two
  events change the model, so they have no linear overlay (the response
  says why in `linear_note`); a phase jump only moves the initial state, so
  it has one. Tested in `tests/test_breakers.py`.

## Reference

```{eval-rst}
.. automodule:: g2elin_api.main
```

```{eval-rst}
.. automodule:: g2elin_api.analysis
```

```{eval-rst}
.. automodule:: g2elin_api.network_routes
```

```{eval-rst}
.. automodule:: g2elin_api.presets
```

```{eval-rst}
.. automodule:: g2elin_api.schemas
```
