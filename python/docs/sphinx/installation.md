# Installation and setup

Everything here runs from a checkout of the repository: there is no package
to install from PyPI, and nothing is written outside the project directory.
The Python side lives in `python/`, and all commands below are run from
there.

## What you need

- **Python 3.10 or newer** ([python.org/downloads](https://www.python.org/downloads/)).
  On Windows, tick *"Add python.exe to PATH"* in the installer.
- About **1.5 GB** of disk space (most of it SciPy, pandapower and numba).
- A browser. Nothing else: no MATLAB, no compiler, no database, no account.

Git is optional — the repository can be downloaded as a ZIP instead.

## The short way

A launcher script at the top of the repository does the whole setup and then
opens the web interface. It is safe to run again: after the first time it
only starts the interface.

- **Windows** — double-click `start-windows.bat`
- **macOS / Linux** — run `./start-macos-linux.sh`

The first run takes a few minutes (it downloads the scientific stack and
builds this manual); afterwards it takes seconds. Both scripts create a
private environment in `python/.venv/` and leave the rest of the machine
untouched. {doc}`getting_started` walks through it in plain language, for
readers who do not work with Python.

## From the command line

```bash
cd python

# 1. A private environment for the project
python -m venv .venv
# Windows:            .venv\Scripts\activate
# macOS / Linux:      source .venv/bin/activate

# 2. The project itself, in editable mode, with what you need
pip install -e ".[api,docs]"
```

`-e` (editable) means the installed package points at the source tree, so
edits take effect without reinstalling.

### Which extras to install

| Extra | Brings in | Needed for |
| --- | --- | --- |
| *(none)* | pydantic, pandapower, numpy/pandas, sympy, scipy, networkx | `g2elin_core` on its own: power flow, modal analysis, EMT, from Python |
| `api` | FastAPI, uvicorn, httpx | the web interface and the HTTP API |
| `docs` | Sphinx, myst-nb, furo, mermaid | building this manual (the interface's Documentation page) |
| `dev` | pytest, pytest-cov | running the test suite |
| `notebook` | JupyterLab, ipykernel, matplotlib | the notebooks in `notebooks/` |

They combine: `pip install -e ".[api,docs,dev,notebook]"` installs
everything.

### Build the manual

```bash
python tools/build_docs.py
```

This writes `docs/sphinx/_build/html/`, which the API serves at `/manual`
and the interface's **Documentation** page shows. Without it that page is
the only thing that will not work; everything else runs regardless.

### Run the web interface

```bash
python -m uvicorn g2elin_api.main:app --port 8000
```

Then open <http://127.0.0.1:8000>. Add `--reload` while editing the code,
and `--host 0.0.0.0` to let other machines on the network reach it.

### Run the tests

```bash
pip install -e ".[dev]"
python -m pytest -q          # the whole suite, a few minutes
python -m pytest -q tests/test_breakers.py     # one file
```

### Use it from Python or a notebook

```python
from g2elin_core.network.presets import wscc9_3sm
from g2elin_core.powerflow import run_power_flow

net = wscc9_3sm()
result = run_power_flow(net)
print(result.bus_table())
```

```bash
pip install -e ".[notebook]"
jupyter lab notebooks/
```

## Keeping it up to date

```bash
git pull
pip install -e ".[api,docs]"   # only if the dependencies changed
python tools/build_docs.py     # only if the documentation changed
```

An editable install needs no reinstall for ordinary code changes.

## If something goes wrong

**`python` is not recognised (Windows).** Python is not on the PATH.
Re-run its installer, choose *Modify*, and tick *"Add python.exe to PATH"* —
or use the `py` launcher (`py -3 -m venv .venv`).

**`running scripts is disabled on this system` (Windows PowerShell).**
Activating an environment from PowerShell is blocked by the execution
policy. Use `cmd.exe`, or call the interpreter directly without activating:
`.venv\Scripts\python -m uvicorn g2elin_api.main:app`.

**`address already in use`.** Something is already on that port — often this
tool, still running in another window. Use another one: `--port 8001`.

**The first power flow or EMT run is slow.** The first call compiles the
symbolic component models (sympy) and pandapower's numba kernels. It is a
one-off per process; later runs are much faster, and the server keeps the
result cached per network.

**The Documentation page shows an error.** The manual has not been built:
run `python tools/build_docs.py` (needs the `docs` extra).

**Deleting everything.** Remove `python/.venv/`. Nothing else is installed
anywhere on the machine.

## Publishing the documentation

The manual is published as a free Hugging Face *Static* Space, in the Read
the Docs theme:

```bash
pip install -e ".[docs]"
hf auth login                                   # once, with a write token
python tools/deploy_docs_space.py --space <username>/g2elin-docs
```

It builds `--theme rtd` into `docs/sphinx/_build/rtd/`, stages it with a
Space card, creates the Space on the first run and updates it afterwards, and
prints the address (`https://<username>-g2elin-docs.static.hf.space`). Two
builds exist because the web interface embeds the default one and reads
Furo's sidebar markup to make its own navigation; the published copy uses
`sphinx_rtd_theme` instead. A static Space is free; the interactive
application needs a backend and so cannot live in one (see below).

## Running it on a server

`python/Dockerfile` builds a container that installs the package, builds the
manual and serves the interface with uvicorn; `render.yaml` at the top of
the repository deploys that image on [Render](https://render.com) as a web
service. `G2ELIN_WARMUP=1` makes the server build a few models at startup so
the first request is not the slow one.
