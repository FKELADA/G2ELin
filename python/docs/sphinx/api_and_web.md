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
- **Modal Analysis** — one linearisation per network version (cached
  server-side too, see `analysis._cached`), seven sub-pages: eigenvalue
  map (+ mode table, click to select a mode), participation heatmap,
  single-mode participation, sensitivity heatmap, mode shape, free-motion
  response and step response. The last two hold *channels* — one
  perturbation (free motion) or one input step (step response) each —
  and every channel has any number of *subplots*, each with its own set of
  signals; crosshairs are synced across a channel's subplots.
- **EMT Simulation** — state offset or input step at T0 = 0, integrated
  with the nonlinear model; plots start 0.01 s before T0 (`t_pre`) so the
  initial point x0 is visible. *Overlay the equivalent linearised
  response* adds the linear model's response to the same disturbance
  (`linear_overlay`: `lsim` of the full `(A, B, C, D)` around the same
  operating point, in absolute units) as dotted lines in the same colours.
  The nonlinear model's initial point is not an exact equilibrium (see
  {doc}`modules/timedomain`), so the EMT traces include a small initial
  transient the linear response doesn't have — the difference
  EMT(disturbed) − EMT(undisturbed) is what matches the linear response
  (checked in `tests/test_analysis_options.py`). *Trace live* streams the
  solver steps as before.

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
