# g2elin-core

The Python core of **G2ELin**, a tool for power-system linearisation and
time-domain simulation — see the [repository README](../README.md) for what
the tool is for, and [the manual](https://fkelada-g2elin-docs.static.hf.space)
for the equations and the API reference.

Two packages:

- **`g2elin_core`** — the compute core. Typed network schema, power flow,
  symbolic component models, linear assembly, modal analysis, model-order
  reduction and nonlinear time-domain integration. No web dependencies;
  usable on its own.
- **`g2elin_api`** — an optional FastAPI layer over it (`[api]` extra) that
  also serves the static web interface from `web/`.

## Setup

```powershell
py -m venv .venv
.venv\Scripts\Activate.ps1
pip install -e ".[dev]"           # the core only
pip install -e ".[dev,api]"       # + the FastAPI backend and web interface
pip install -e ".[dev,api,docs]"  # + Sphinx, to build the manual
pip install -e ".[notebook]"      # + Jupyter and matplotlib, for notebooks/tour.ipynb
```

## Run it

```powershell
uvicorn g2elin_api.main:app --reload
```

Then open `http://127.0.0.1:8000/`. One process serves both the JSON API and
the interface. The first modal-analysis request for a given model takes a few
seconds (the symbolic models are derived once, then cached); everything after
that is fast.

From the repository root, `start-windows.bat` does the whole thing —
environment, install, docs, browser — in one step.

## Tests

```powershell
pytest -q
```

## Documentation

```powershell
python tools/build_docs.py              # -> docs/sphinx/_build/html, served by the interface
python tools/build_docs.py --theme rtd  # -> docs/sphinx/_build/rtd, the published copy
python tools/check_docs_facts.py        # the numbers the prose states, against the code
```

## Notebook

`notebooks/tour.ipynb` exercises every module directly, in the order data
flows through the tool, plotting each result. It is committed with its
outputs, so it reads without being run. The EMT section is off by default
(`RUN_EMT_SECTION`) because it is much the slowest cell.

## Layout

```
src/g2elin_core/
  network/            typed schema, preset cases, 2D layout, structural validation
  powerflow/          pandapower adapter
  components/         symbolic linear + nonlinear models (SM, GFM, GFL, IB, line, node, load)
  operating_point.py  solved power flow -> per-component linearisation points
  interconnect/       block/wiring topology + closed-loop assembly
  reduction.py        the model-order catalogue: what each element may give up
  modal/              eigenvalues, participation, classification, sensitivity, adequacy
  timedomain/         nonlinear (EMT) time-domain integration
  timeseries/         time-series load flow
  pipeline.py         end to end: power flow -> linearised closed-loop model
  pu_base.py          AC/DC per-unit base values
src/g2elin_api/
  main.py             preset-keyed endpoints and app wiring
  network_routes.py   /api/network/* -- any analysis on an arbitrary Network JSON
  analysis.py         the analysis logic both of those share
  presets.py          the named preset catalogue
web/                  the static frontend (plain JS)
notebooks/            tour.ipynb
tools/                docs build, preset import, cross-validation, deployment
  ui_tour/            records the captioned video tour of the interface
tests/                the test suite; tests/golden/ holds the MATLAB regression fixtures
docs/
  sphinx/             the manual
  development-log.md  how this was built, entry by entry
```
