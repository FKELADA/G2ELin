"""Sphinx configuration for the g2elin-core / g2elin-api documentation.

Run with ``python tools/build_docs.py`` (from ``python/``) or directly:
``sphinx-build -b html docs/sphinx docs/sphinx/_build/html``.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

project = "G2ELin"
copyright = "2026, Fadi Kelada"
author = "Fadi Kelada"
release = "0.1.0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "myst_nb",  # supersedes myst_parser: same markdown parser, plus .ipynb rendering
    "sphinxcontrib.mermaid",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "myst-nb",
    ".ipynb": "myst-nb",
}

# dollarmath/amsmath: lets every module page write physics as real LaTeX
# (`$...$` inline, `$$...$$` display, `\begin{aligned}...\end{aligned}`
# blocks) instead of prose paraphrase — rendered client-side by MathJax.
# myst_nb reads this same myst_enable_extensions setting (it wraps
# myst_parser rather than replacing its config surface).
myst_enable_extensions = ["colon_fence", "deflist", "dollarmath", "amsmath"]
# Auto-generate heading anchor ids (so "#mount-layout"-style in-page links
# resolve) and don't try to resolve plain markdown [text](path) links as
# internal cross-references — a doc source can link to non-.md source
# files by relative path, which isn't a doc target Sphinx knows about.
myst_heading_anchors = 3
myst_all_links_external = True

# Notebooks embedded in the docs (see cigre_walkthrough_notebook.md) are
# executed and verified separately (jupyter nbconvert --execute, checked
# for zero cell errors) before being staged into this source tree by
# tools/build_docs.py -- "off" renders exactly those already-baked
# outputs rather than re-running a multi-minute EMT simulation on
# every docs build.
nb_execution_mode = "off"

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# Autodoc: show type hints inline, keep source order (matches how each module's
# own docstring narrates it), skip re-exported __init__ noise.
autodoc_typehints = "description"
autodoc_member_order = "bysource"
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
}
napoleon_google_docstring = True
napoleon_numpy_docstring = False

# pandapower/sympy/scipy are heavy optional runtime deps that pull in compiled
# extensions; mock them out so autodoc can import g2elin_core modules for their
# docstrings/signatures without needing a fully built environment.
autodoc_mock_imports = []

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
}

html_theme = "furo"
html_title = "G2ELin documentation"
html_static_path = ["_static"] if Path(__file__).resolve().parent.joinpath("_static").is_dir() else []

# Mermaid diagrams render client-side via mermaid.js (loaded from a CDN by
# sphinxcontrib-mermaid) — no Node.js/mermaid-cli needed at build time, which
# matters in this environment (see README: no Node.js available).
mermaid_version = "10.9.1"
