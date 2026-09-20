"""Builds the Sphinx documentation to ``docs/sphinx/_build/html``.

Run from ``python/`` after ``pip install -e ".[docs]"``:

    python tools/build_docs.py

Thin wrapper around ``sphinx-build`` (rather than requiring ``make``, which
isn't available out of the box on Windows) so this works the same way in
any shell. ``g2elin_api.main`` mounts the output directory at ``/manual``
if it exists — see ``docs/sphinx/api_and_web.md``.

Also stages any notebook this documentation embeds (via ``myst_nb``,
``nb_execution_mode = "off"`` in ``conf.py``) from the project's canonical
``notebooks/`` directory into ``docs/sphinx/_notebooks/`` first.
``notebooks/`` stays the single source of truth (same place
``tour.ipynb``/``random_network.ipynb`` live, git-tracked, opened directly
by a reader with Jupyter) — the staged copy is a build artifact, not a
second copy to keep in sync by hand, the same way ``_build/`` itself is.
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys
from pathlib import Path

from sphinx.cmd.build import main as sphinx_main

ROOT = Path(__file__).resolve().parents[1]
SPHINX_SRC = ROOT / "docs" / "sphinx"
OUT_DIR = SPHINX_SRC / "_build" / "html"
NOTEBOOKS_STAGING_DIR = SPHINX_SRC / "_notebooks"

# Notebooks this documentation embeds directly (rendered in place via
# myst_nb, not just linked) -- add a name here when a new one is wired
# into index.md's toctree. tour.ipynb is the general-purpose
# exploratory notebook and stays link-only (see README.md), not embedded.
EMBEDDED_NOTEBOOKS = [
    "cigre_1sm_1gfm_1gfl_walkthrough.ipynb",
    "random_network.ipynb",
]


def _stage_notebooks() -> None:
    NOTEBOOKS_STAGING_DIR.mkdir(exist_ok=True)
    for name in EMBEDDED_NOTEBOOKS:
        shutil.copy2(ROOT / "notebooks" / name, NOTEBOOKS_STAGING_DIR / name)


THEMES = {
    # name -> (Sphinx theme, output directory). Furo is what the web interface
    # embeds (its Documentation page reads Furo's sidebar); the Read the Docs
    # theme is the published copy (tools/deploy_docs_space.py).
    "furo": ("furo", OUT_DIR),
    "rtd": ("sphinx_rtd_theme", SPHINX_SRC / "_build" / "rtd"),
}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build the G2ELin documentation.")
    parser.add_argument("--theme", choices=sorted(THEMES), default="furo",
                        help="furo (default; embedded in the web interface) or rtd (the published copy)")
    args_ = parser.parse_args(argv)
    theme, out = THEMES[args_.theme]
    os.environ["G2ELIN_DOCS_THEME"] = theme
    _stage_notebooks()
    # Each theme has its own output directory, and Sphinx keeps its build cache
    # inside it, so switching theme never reuses the other one's pages.
    args = ["-b", "html", "-W", "--keep-going", str(SPHINX_SRC), str(out)]
    return sphinx_main(args)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
