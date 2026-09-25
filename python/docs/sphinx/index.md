# G2ELin documentation

**G2ELin** is an open-access tool for the small-signal analysis and
time-domain simulation of power systems with a high share of inverter-based
resources.
It builds a network of synchronous machines, grid-forming and grid-following
converters, infinite buses, lines, transformers and loads; solves its power
flow; linearises the whole thing into one closed-loop state-space model for
modal analysis (eigenvalues, participation factors, root loci); and
integrates the same model nonlinearly in time.

What you can do with it:

- **Build or load a network** — 24 presets (the CIGRE MV benchmark islanded and
  interconnected, WSCC 9-bus, single- and two-machine cases, Kundur's two-area
  system, and IEEE 14, 39 and 118-bus), or draw one from scratch, with
  breakers on every branch, load and unit.
- **Solve power flow** — several solvers, single or batch, on networks that
  may be split into islands.
- **Choose the models** — a machine's exciter, stabiliser and governor, and a
  converter's power-control law (droop, filtered droop, dVOC, VSM, matching),
  per unit; and how much of each element's dynamics to keep, which is what
  makes the same case EMT or RMS.
- **Analyse modes** — eigenvalue maps, participation, sensitivity (including
  which *physical parameter* moves a mode), mode shapes, free and step
  responses, and root loci that sweep any parameter, including each control
  loop's tuning.
- **Simulate in time** — nonlinear runs with network events (breaker
  openings, load steps, phase jumps) and the linearised response overlaid.
- **Choose the model order** — keep every dynamic for an EMT run, or make
  the network quasi-stationary and the units lower-order for an
  electromechanical (RMS) one, with a check that tells you whether the
  reduction is safe for *your* network. See {doc}`modules/reduction`.

It runs as a Python library (`g2elin_core`), as a web interface
(`g2elin_api` + `web/`), or from the notebooks. New here? Start with
{doc}`getting_started`; the {doc}`installation` page has the command-line
version. The rest of this manual is the architecture and API reference,
generated from the code with [Sphinx](https://www.sphinx-doc.org/) and
[autodoc](https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html),
with [Mermaid](https://mermaid.js.org/) diagrams for the parts that are
easier to see than to read as prose.

Source: [github.com/FKELADA/G2ELin](https://github.com/FKELADA/G2ELin) ·
licence: GPL-3.0.

## How the pieces fit together

```{mermaid}
flowchart TD
    subgraph Core["g2elin_core (compute core, no web deps)"]
        NET["network<br/>schema + presets + topology"]
        PU["pu_base"]
        PF["powerflow<br/>pandapower_adapter"]
        OP["operating_point"]
        COMP["components<br/>base + sm/gfm/gfl/line/node/load/ib"]
        IC["interconnect<br/>assemble + network_assembly"]
        PIPE["pipeline"]
        MODAL["modal<br/>analysis"]
        TS["timeseries<br/>scenario"]
        TD["timedomain<br/>emt"]
    end
    subgraph Surface["g2elin_api + web/"]
        API["g2elin_api<br/>FastAPI"]
        WEB["web/index.html<br/>vanilla JS UI"]
        DOCS["this documentation<br/>served at /manual"]
    end

    NET --> PF
    NET --> TS
    PF --> TS
    PU --> OP
    PF --> OP
    OP --> COMP
    NET --> IC
    COMP --> IC
    IC --> PIPE
    IC --> TD
    PIPE --> MODAL

    NET --> API
    PF --> API
    PIPE --> API
    MODAL --> API
    TS --> API
    API --> WEB
    DOCS --> WEB

    click NET "modules/network.html" "network schema + presets + topology"
    click PF "modules/powerflow.html" "static power flow"
    click OP "modules/operating_point.html" "operating point"
    click COMP "modules/components.html" "symbolic component models"
    click IC "modules/interconnect.html" "block interconnection"
    click PIPE "modules/pipeline.html" "linearization pipeline"
    click MODAL "modules/pipeline.html" "modal analysis"
    click TS "modules/timeseries.html" "time-series load flow"
    click TD "modules/timedomain.html" "EMT/nonlinear simulation"
    click API "api_and_web.html" "FastAPI + web UI"
```

Reading it left to right: a {doc}`Network <modules/network>` (typed buses,
lines, transformers, loads, DER units — see `g2elin_core.network.schema`)
feeds {doc}`static power flow <modules/powerflow>` via pandapower, which
produces bus voltages/angles that {doc}`operating_point <modules/operating_point>`
turns into per-component operating points. Those operating points parameterize
the {doc}`symbolic component models <modules/components>` (synchronous
machine, grid-forming/grid-following converters, lines, nodes, loads), which
{doc}`interconnect <modules/interconnect>` wires into one closed-loop system —
consumed either by the linear {doc}`pipeline and modal analysis <modules/pipeline>`
or by the nonlinear {doc}`EMT time-domain solver <modules/timedomain>`.
{doc}`g2elin_api and the
web UI <api_and_web>` sit on top, exposing all of this over HTTP to the
browser tab you're probably reading this from. For a single worked
example that runs every one of these modules end to end on one real
network, with actual computed results and the physical interpretation of
each one, see {doc}`cigre_walkthrough` — or {doc}`_notebooks/cigre_1sm_1gfm_1gfl_walkthrough`
for the same walkthrough as the actual notebook it's written from,
rendered here cell by cell (code, narration, and every plot), not just
the curated highlights. For the reverse direction — building a `Network`
from scratch by hand, respecting the same conventions {doc}`its own
validator <modules/network>` checks for — see
{doc}`_notebooks/random_network`.

## Contents

```{toctree}
:maxdepth: 2
:caption: Getting started

getting_started
installation
```

```{toctree}
:maxdepth: 1
:caption: Validation

validation/cross-validation-2026-09
```

```{toctree}
:maxdepth: 2
:caption: Architecture

modules/network
modules/powerflow
modules/operating_point
modules/components
modules/interconnect
modules/pipeline
modules/reduction
modules/timeseries
modules/timedomain
api_and_web
cigre_walkthrough
_notebooks/cigre_1sm_1gfm_1gfl_walkthrough
_notebooks/random_network
```

## Status at a glance

| Phase | Feature | State |
|---|---|---|
| P0 | Foundations (schema, pu-base, presets) | done |
| P1 | Static power flow (pandapower) | done |
| P2 | Linear small-signal / modal analysis | done — SM (selectable exciter/PSS/governor), GFM (five power-control laws), GFL, IB slack; full toolbox (sensitivity, mode shape, free/step response) |
| P3 | Time-series load flow | done |
| P4 | EMT / nonlinear time-domain simulation | done — analytic Newton Jacobian |
| P7 | Web UI | presets + 5 result tabs (incl. Network) + this doc tab, always-light theme |

This table (and everything linked from it) describes what's implemented and
tested *in this codebase*, not aspirations — see each module page's own
"what this doesn't do" notes where relevant, and `python/README.md` for the
full phase-by-phase account including bugs found and fixed along the way.
```
