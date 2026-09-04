# `g2elin_api` and the web UI

`g2elin_api` is a thin FastAPI layer wrapping the real `g2elin_core`
functions (no mocking) over the two validated presets. `web/index.html`
is a build-free, vanilla-JS single-page app that calls it. This page (the
one you're reading, if you got here through the **Documentation** tab) is
itself served as part of that same app — see [Mount layout](#mount-layout)
below.

## Endpoints → tabs

```{mermaid}
flowchart LR
    subgraph API["g2elin_api (FastAPI)"]
        EP1["GET /api/presets"]
        EPT["GET /api/presets/{id}/topology"]
        EP2["POST /api/presets/{id}/powerflow"]
        EP3["POST /api/presets/{id}/modal"]
        EP3S["POST .../modal/sensitivity"]
        EP3M["POST .../modal/mode_shape"]
        EP3F["POST .../modal/free_response"]
        EP3R["POST .../modal/step_response"]
        EP4["POST /api/presets/{id}/timeseries"]
        EP5["GET /api/presets/{id}/states"]
        EP6["POST /api/presets/{id}/emt"]
        EP6L["POST .../emt/live<br/>(streamed NDJSON)"]
        EP7["POST /api/presets/{id}/roa"]
    end
    subgraph UI["web/index.html tabs"]
        T0["Preset picker"]
        T1["Power Flow"]
        T2["Modal Analysis"]
        T3["Time Series"]
        TN["Network"]
        T4["EMT Simulation"]
        T5["Region of Attraction"]
        T6["Documentation"]
    end
    EP1 --> T0
    EPT --> TN
    EPT --> T1
    EP2 --> T1
    EP3 --> T2
    EP3S --> T2
    EP3M --> T2
    EP3F --> T2
    EP3R --> T2
    EP4 --> T3
    EP5 --> T4
    EP5 --> T5
    EP6 --> T4
    EP6L --> T4
    EP7 --> T5

    EPT -.-> TOPO["network.compute_topology_layout()"]
    EP2 -.-> PF["powerflow.run_power_flow()"]
    EP3 -.-> PIPE["pipeline.linearize_network()<br/>+ modal.analyze()"]
    EP3S -.-> SENS["modal.eigenvalue_sensitivity()"]
    EP3M -.-> MSH["modal.mode_shape()"]
    EP3F -.-> FREE["modal.free_response()"]
    EP3R -.-> STEP["modal.step_response()"]
    EP4 -.-> TS["timeseries.run_time_series()"]
    EP5 -.-> BNM["timedomain.build_nonlinear_network()"]
    EP6 -.-> EMT["timedomain.simulate()<br/>+ recover_inputs_and_outputs() (opt-in)"]
    EP6L -.-> EMTLIVE["timedomain.simulate_steps()<br/>one solver step per streamed line"]
    EP7 -.-> ROA["stability.trace_roa_grid()<br/>(a grid of trajectories)"]
```

## `/api/network/*` — arbitrary networks, not just presets

Every endpoint above also has a counterpart under `/api/network/*`
(`network_routes.py`) taking a full `Network` JSON body instead of a
`preset_id` path param — `POST /api/network/powerflow`, `/modal`,
`/modal/sensitivity`, `/modal/mode_shape`, `/modal/free_response`,
`/modal/step_response`, `/timeseries`, `/emt`, `/emt/live`, `/roa`, plus `/topology` and
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
whatever `dict()` happened to produce. `RoaResponse`'s distance grids use
`float | None` rather than `float`, specifically because JSON has no `NaN`
literal — the underlying `RoaGridResult.early_distance`/`late_distance`
arrays are `NaN` at any grid point where the coupled Newton solve failed
(see {doc}`stability <modules/stability>`), and `None`/`null` is what
survives the round trip to the browser instead of either crashing
`JSON.stringify` server-side (Python's `json` module *would* emit a bare
`NaN` token, which is invalid JSON and fails `fetch().json()` in the
browser) or silently coercing a "we don't know" into a number.
`PowerFlowResponse`'s `lines`/`transformers`/`loads`/`generators`/
`static_generators`/`external_grid` fields are loosely typed (`list[dict]`,
not a row model per table) since each pandapower result table's column set
genuinely differs (13 columns for lines, 2 for loads) and this is an
inspection/dropdown feature, not something else in this codebase consumes
programmatically.

`/states` exists because the nonlinear model's state/input/output names
are entirely preset-dependent — which DER id is the slack, which unit
types exist, whether a `dw_r_*` (synchronous-machine speed deviation)
state exists at all — so the EMT/ROA tabs' pickers are populated from
this endpoint at preset-selection time rather than hard-coded, the same
problem {doc}`the ROA module itself <modules/stability>` solves with
`find_state_index()`'s substring matching (which `/roa` and the EMT
tab's *perturb*-state selection both reuse server-side). `/emt`'s
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

Seven tabs: **Power Flow**, **Modal Analysis**, **Time Series**,
**Network**, **EMT Simulation**, **Region of Attraction**, and
**Documentation** (this Sphinx site, lazily loaded in an `<iframe>` on
first click so the page's initial load doesn't pay for it). Always light —
the dark-mode CSS block that used to shadow-follow `prefers-color-scheme`
was removed outright rather than kept as an unused branch, so every chart's
color logic only has to reason about one palette.

Every tab now carries a short `.panel-intro` card right below its
controls — what the tab computes, the underlying method in one line, and
a concrete "useful for" case, written for someone learning the physics
rather than already knowing it — and a small hoverable/tappable "(i)"
bubble next to each tab's own label in the nav bar with a one-sentence
version of the same, so the summary is visible before a tab is even
opened. Both are self-contained CSS/JS (`.info-icon`/`.tooltip` for the
bubble, `.panel-intro` for the card); the bubble's tooltip opens on
`:hover`/`:focus-visible` and also toggles on click/tap (`.open` class,
closed by a document-level click listener) so it's reachable without a
mouse. The full underlying physics each summary is a one-line version of
lives in this Sphinx site itself (mostly {doc}`modules/components`,
{doc}`modules/operating_point`, {doc}`modules/pipeline`,
{doc}`modules/timedomain`, and {doc}`modules/stability` — see each for
the real LaTeX equations), which is exactly what the Documentation tab
embeds.

The eigenvalue map (SVG, symlog-scaled since eigenvalues here span ~15
orders of magnitude) uses the dataviz skill's status-palette convention for
stable/marginal/unstable points, now with the constant-damping-ratio guide
lines (5%, 70.7%) `Functions/modal_analysis.m` also draws — straight lines
in linear (real, freq) space, sampled log-spaced in `|real|` and
symlog-transformed like every other point on the chart, so they render as
curves here. It's also zoomable/pannable (`attachSvgZoomPan()`, a small
reusable viewBox-manipulation helper — mouse wheel zooms centered on the
cursor, drag pans, a "Reset zoom" button restores the original view) with a
full x-y gridline set (`eigenGridLines()`, not just the origin crosshair),
tick labels inverse-transformed back to physical real-part/Hz units since
that's what a reader actually wants to read off the axis. The mode summary
table lists each mode's top **three** participating states (`state1..3`/
`part1_pct..3`, already computed by `modal.summary_table()`), matching
`modal_analysis.m`'s own table layout. Every other chart on these tabs
follows the same validated categorical/sequential/status palette
(`references/palette.md`, already used as this file's `--series-*`/
`--good`/`--critical`/sequential-blue CSS custom properties):

- **Network** (and the **Power Flow** tab, reusing the same function)
  render `network.compute_topology_layout()`'s bus/line/transformer graph
  as an SVG diagram — buses as circles (colored by unit type on the
  Network tab: fixed categorical order SM/GFM/GFL/IB; colored by voltage
  magnitude on the Power Flow tab: sequential blue, fixed [0.9, 1.1] pu
  domain so a color is comparable across presets, not just within one),
  transformers as a small square marker on their edge, loads as a
  triangle marker. Every bus/line/transformer is also clickable (not just
  hoverable) — the topology response now carries each line's `r_pu`/
  `x_pu`/`b_pu`/`length_km`, each transformer's `r_pu`/`x_pu`/`sn_mva`, and
  each DER-hosting bus's own dispatch fields plus its *derived*
  electrical/control-parameter dict (`sm_params()`/`gfm_params()`/
  `gfl_params()`'s actual output — the same values used to build that
  unit's dynamic model, not just the static schema fields), rendered into
  a detail panel below the diagram on click. Power Flow also gained a
  dropdown over pandapower's other result tables (lines/transformers/
  loads/generators/static generators/external grid), reusing a
  generalized version of the existing bus-table renderer rather than one
  function per table shape. **Time Series** got the identical
  diagram-plus-dropdown-plus-click-detail treatment, once per snapshot —
  the network's topology doesn't change across load-scaling snapshots
  (only P/Q values do), so the layout is fetched once and just recolored
  per snapshot's own bus results, not refetched.
- **Modal Analysis** gained a "View" dropdown covering the rest of
  `modal_analysis.m`'s toolbox: participation heatmap, single-mode
  participation bar chart, sensitivity heatmap (capped to state counts
  ≤30 — a 87×87 HTML table isn't a readable heatmap, the top-8 list from
  the API is the useful output past that size), mode shape (polar plot,
  phase angle of the top-5 participating states), free motion response
  and step response (both reuse the EMT tab's line-chart renderer, since
  they're the same "value vs time" form).
- **EMT Simulation** plots the perturbed trajectory as up to three
  separate multi-series line charts — states, inputs, outputs — rather
  than one shared axis: they're physically different kinds of quantity
  (per-unit dynamics, fixed exogenous references, named physical outputs
  like P/Q/V), so combining them on one y-axis would risk exactly the
  "two measures of different scale" dual-axis problem the dataviz skill
  rules out. The chart itself is built directly via the SVG DOM API (not
  string templates, since it needs live `mousemove` listeners for a real
  crosshair+tooltip) rather than the `innerHTML` string-building most of
  the rest of this file uses; the three state/input/output pickers are
  genuine `<select multiple>` elements over the full list from `/states`,
  not a text filter, per the explicit ask this round. What can be
  perturbed was later widened from states-only to states *or* inputs — a
  single merged `<select>` (grouped by `<optgroup>`, kind read off the
  chosen `<option>`'s `data-kind`) — since an input step (a permanent
  change to an exogenous reference like `P_ref`, held from t=0) is the
  more standard disturbance test and `timedomain.simulate()`'s `u_exo_fn`
  parameter already existed for exactly this, just unused until now. A
  "Timestep" field controls the output sampling grid (`dt`, not the
  solver's own adaptive internal step size — see {doc}`timedomain
  <modules/timedomain>`), bounded server-side (10–2000 samples) since
  each sample costs a Newton solve when inputs/outputs are requested too.
- **Region of Attraction** renders the grid as an HTML `<table>` (a natural
  fit for a small labeled matrix — row/column headers give the axis
  offsets "for free"), each cell colored by status: good (trending toward
  baseline), critical (trending away), or muted gray (Newton solve failed
  to converge — a distinct "we don't know," never silently folded into
  either verdict, matching {doc}`the ROA module's own design decision
  <modules/stability>`). A native `title` attribute per cell gives the
  early/late distance on hover, the same lightweight-tooltip precedent the
  eigenvalue map already set with per-point `<title>` elements.

### The network editor

An "Edit this network" button next to the preset picker clones the
selected preset's `Network` (via `GET .../network`) into `state.customNetwork`
and shows a slide-down `#network-editor` panel — not a dedicated tab, since
editing has to affect every existing tab identically, which a tab would
only awkwardly achieve. Two dispatcher functions, `runAnalysis()` and
`fetchTopologyActive()`, route every analysis trigger in the file to
`/api/network/*` with `state.customNetwork` as the body whenever it's set,
falling back to the ordinary preset-id endpoints otherwise — every tab's
own button/spinner/error handling is unchanged, only *where the request
goes* differs.

The editor itself is table-per-element-kind (buses/lines/transformers/
loads/DER units) with inline-editable `<input>`/`<select>` cells, reusing
the same `.tablewrap`/`table` CSS every other data table in this file
uses — bulk-editing dozens of fields across 5 element kinds is a CRUD
table's job, not the click-in-diagram detail panel's (built for
inspecting one clicked element, not editing many). The same debounce also
POSTs to `POST /api/network/validate` (see {doc}`network
<modules/network>`'s "Validation" section) and renders every issue it
finds into its own panel — every structural problem at once (a DER not
behind a transformer, a Line/Load on a DER's own bus, an unsupported
slack unit type, ...), not discovered one crash at a time by actually
running an analysis. Separately, a debounced (~400ms) preview re-POSTs
the in-progress network to `/api/network/topology` and re-renders
`renderNetworkDiagram()`; a network that's momentarily invalid mid-edit
(e.g. a line still referencing a just-deleted bus) shows as a plain "not
valid yet" message there, not an error — the validation panel above
already explains what's wrong in detail. A "Recalculate" button
re-clicks whichever analysis tab is currently active, reusing its existing
run button entirely (see `recalculateActiveTab()`). "Reset to preset"
clears `state.customNetwork` and restores the ordinary preset flow —
nothing is ever saved, so there's no discard confirmation, and picking a
different preset from the dropdown while editing exits edit mode the same
way.

FastAPI's automatic 422 for an invalid `Network` body has `detail` as a
*list* of `{loc, msg, type}` objects (pydantic v2's shape), unlike every
hand-written `HTTPException` elsewhere in this app (a plain string) — the
`api()` helper's `formatDetail()` renders that list as one readable
message instead of stringifying it to `"[object Object],..."`.

**Deliberately plain JS, not React** — the migration plan's original
target stack is React/TypeScript via Vite, but this environment has no
Node.js/npm, so there's no way to install or build one. This ships as a
build-free page FastAPI serves directly as a static file — real, working,
and swappable for a proper Vite+React+TS app later without touching the
API contract it calls.

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
