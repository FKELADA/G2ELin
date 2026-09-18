"""Publishes the G2ELin web app as a Hugging Face Docker Space.

Run from ``python/`` after building the docs (``python tools/build_docs.py``)
and logging in to Hugging Face once (``hf auth login``, with a token that
has *write* access):

    python tools/deploy_hf_space.py --space <your-hf-username>/g2elin

The Space is created on the first run and updated on later ones; it then
builds the Docker image itself and serves the app at
``https://<username>-g2elin.hf.space``. ``--stage-only DIR`` assembles the
upload in DIR without uploading, e.g. to inspect it or to ``docker build``
it locally.

What's uploaded is only what the image needs: ``pyproject.toml``,
``src/``, ``web/``, the built Sphinx HTML (without Sphinx's ``.doctrees``
build cache), and ``deploy/hf_space/``'s ``Dockerfile`` + ``README.md`` (the
Space card, which also satisfies ``pyproject.toml``'s ``readme``).
"""

from __future__ import annotations

import argparse
import shutil
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCS_HTML = ROOT / "docs" / "sphinx" / "_build" / "html"
SPACE_FILES = ROOT / "deploy" / "hf_space"
IGNORE = shutil.ignore_patterns("__pycache__", "*.pyc", "*.egg-info", ".doctrees", ".DS_Store")


def stage(dest: Path) -> None:
    if not (DOCS_HTML / "index.html").is_file():
        sys.exit("Built docs not found -- run `python tools/build_docs.py` first.")
    dest.mkdir(parents=True, exist_ok=True)
    shutil.copy2(ROOT / "pyproject.toml", dest / "pyproject.toml")
    shutil.copy2(SPACE_FILES / "Dockerfile", dest / "Dockerfile")
    shutil.copy2(SPACE_FILES / "README.md", dest / "README.md")
    shutil.copytree(ROOT / "src", dest / "src", ignore=IGNORE, dirs_exist_ok=True)
    shutil.copytree(ROOT / "web", dest / "web", ignore=IGNORE, dirs_exist_ok=True)
    shutil.copytree(DOCS_HTML, dest / "docs" / "sphinx" / "_build" / "html", ignore=IGNORE, dirs_exist_ok=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--space", help="Space id, <username>/<space-name> (e.g. fkelada/g2elin)")
    parser.add_argument("--private", action="store_true", help="create the Space as private (first run only)")
    parser.add_argument("--stage-only", metavar="DIR", help="only assemble the upload in DIR; don't upload")
    args = parser.parse_args()

    if args.stage_only:
        stage(Path(args.stage_only))
        print(f"Staged in {args.stage_only}")
        return 0
    if not args.space or "/" not in args.space:
        parser.error("--space <username>/<space-name> is required")

    try:
        from huggingface_hub import HfApi
    except ImportError:
        sys.exit('huggingface_hub is not installed -- `pip install huggingface_hub`')

    api = HfApi()
    try:
        user = api.whoami()["name"]
    except Exception:  # noqa: BLE001 -- any auth failure looks the same to the user
        sys.exit("Not logged in to Hugging Face -- run `hf auth login` (token with write access) first.")
    print(f"Logged in as {user}")

    api.create_repo(args.space, repo_type="space", space_sdk="docker", private=args.private, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        stage(Path(tmp))
        api.upload_folder(
            folder_path=tmp,
            repo_id=args.space,
            repo_type="space",
            commit_message="Deploy G2ELin web app",
            # Mirror the local tree: drop files deleted since the last deploy.
            delete_patterns=["src/*", "web/*", "docs/*"],
        )
    owner, name = args.space.split("/")
    print(f"Uploaded. Build logs: https://huggingface.co/spaces/{args.space}")
    print(f"App (once built):     https://{owner.lower()}-{name.lower().replace('_', '-')}.hf.space")
    return 0


if __name__ == "__main__":
    sys.exit(main())
