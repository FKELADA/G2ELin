"""Publishes the documentation as a Hugging Face *Static* Space.

Run from ``python/`` after logging in to Hugging Face once (``hf auth login``,
with a token that has *write* access):

    python tools/deploy_docs_space.py --space <your-hf-username>/g2elin-docs

It builds the Read the Docs-themed copy of the manual
(``build_docs.py --theme rtd``), stages it with a Space card, creates the Space
on the first run and updates it on later ones, and prints the address it is
served at: ``https://<username>-<space-name>.static.hf.space``.

Why a Static Space and not the Docker one ``deploy_hf_space.py`` builds: the
manual is plain HTML, and static Spaces are free, while a Docker Space needs a
paid plan. The interactive application itself needs a backend and cannot be
static -- see ``render.yaml`` for that.

``--stage-only DIR`` assembles the upload in DIR without uploading, e.g. to
open it locally first.
"""

from __future__ import annotations

import argparse
import shutil
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RTD_HTML = ROOT / "docs" / "sphinx" / "_build" / "rtd"
IGNORE = shutil.ignore_patterns(".doctrees", ".buildinfo", ".DS_Store")

CARD = """---
title: G2ELin documentation
emoji: ⚡
colorFrom: blue
colorTo: indigo
sdk: static
pinned: false
license: gpl-3.0
short_description: Small-signal and EMT analysis of power systems
---

Documentation of **G2ELin**: power flow, modal analysis, root loci and EMT
simulation of networks with synchronous machines, grid-forming and
grid-following converters.

Source and installation: <https://github.com/FKELADA/G2ELin>
"""


def build() -> None:
    sys.path.insert(0, str(ROOT / "tools"))
    import build_docs

    code = build_docs.main(["--theme", "rtd"])
    if code:
        sys.exit(f"The documentation build failed (exit {code}).")


def stage(dest: Path) -> None:
    if not (RTD_HTML / "index.html").is_file():
        sys.exit("Built docs not found -- run `python tools/build_docs.py --theme rtd` first.")
    dest.mkdir(parents=True, exist_ok=True)
    shutil.copytree(RTD_HTML, dest, ignore=IGNORE, dirs_exist_ok=True)
    (dest / "README.md").write_text(CARD, encoding="utf-8")
    # Sphinx copies its sources to _sources/ for "View page source"; the
    # theme links to them, so keep them. Underscore-prefixed folders are
    # served as-is by a static Space (no Jekyll to strip them).


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--space", help="Space id, <username>/<space-name> (e.g. FKelada/g2elin-docs)")
    parser.add_argument("--private", action="store_true", help="create the Space as private (first run only)")
    parser.add_argument("--no-build", action="store_true", help="reuse the last rtd build instead of rebuilding")
    parser.add_argument("--stage-only", metavar="DIR", help="only assemble the upload in DIR; don't upload")
    args = parser.parse_args()

    if not args.no_build:
        build()

    if args.stage_only:
        stage(Path(args.stage_only))
        print(f"Staged in {args.stage_only}")
        return 0
    if not args.space or "/" not in args.space:
        parser.error("--space <username>/<space-name> is required")

    from huggingface_hub import HfApi

    api = HfApi()
    who = api.whoami()["name"]
    print(f"Logged in as {who}")
    api.create_repo(
        args.space, repo_type="space", space_sdk="static", private=args.private, exist_ok=True,
    )
    with tempfile.TemporaryDirectory() as tmp:
        stage(Path(tmp))
        api.upload_folder(
            folder_path=tmp, repo_id=args.space, repo_type="space", delete_patterns="*",
            commit_message="Update the documentation",
        )
    owner, name = args.space.split("/", 1)
    print(f"Space:  https://huggingface.co/spaces/{args.space}")
    print(f"Served: https://{owner.lower()}-{name.lower()}.static.hf.space")
    return 0


if __name__ == "__main__":
    sys.exit(main())
