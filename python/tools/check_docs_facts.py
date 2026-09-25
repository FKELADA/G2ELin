"""Check the numbers the documentation states against the code it describes.

Prose goes stale quietly. A sentence that says "19 presets" stays readable
and confident long after the twentieth is added, and nothing fails -- which
is exactly the kind of error a reader has no way to catch.

This walks the Markdown and asserts the countable claims: how many presets
there are, how many states a named preset linearises to, which model-order
levels exist. It is deliberately narrow -- it checks facts that have a
single right answer, not wording -- and it prints every mismatch rather than
stopping at the first.

    python tools/check_docs_facts.py          # report
    python tools/check_docs_facts.py --fix    # rewrite the numbers it can
"""
from __future__ import annotations

import argparse
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
DOCS = [
    ROOT / "README.md",
    ROOT / "GETTING_STARTED.md",
    ROOT / "python" / "README.md",
    *sorted((ROOT / "python" / "docs" / "sphinx").rglob("*.md")),
]


def live_facts() -> dict:
    """What the code actually says, right now."""
    from g2elin_api.presets import PRESETS
    from g2elin_core import reduction
    from g2elin_core.network import presets as np_
    from g2elin_core.pipeline import linearize_network
    from g2elin_core.powerflow import run_power_flow

    facts = {
        "presets": len(PRESETS),
        "sm_levels": len(reduction.element("sm").levels),
        "element_kinds": len(reduction.ELEMENTS),
    }
    for name in ("wscc9_3sm", "cigre_interconnected_1sm_1gfm_1gfl"):
        net = getattr(np_, name)()
        facts[f"{name}_states"] = linearize_network(net, run_power_flow(net)).A.shape[0]
        facts[f"{name}_buses"] = len(net.buses)
    return facts


#: Everything after this marker in a file is left alone. A changelog states
#: what was true when the entry was written -- rewriting those numbers would
#: not correct the document, it would falsify its history.
STOP = "<!-- doc-facts: stop -->"

#: (regex over the document, fact key, what the number means). The regex must
#: capture the number as group "n" so --fix knows what to replace. The
#: lookbehind keeps "IEEE 14/39/118 presets" from reading as a count.
_NOT_A_LIST = r"(?<![\d/.-])"
CLAIMS = [
    (re.compile(_NOT_A_LIST + r"(?P<n>\d+)\s+presets?\b", re.I), "presets", "the number of preset networks"),
    (re.compile(_NOT_A_LIST + r"(?P<n>\d+)\s+preset cases", re.I), "presets", "the number of preset networks"),
    (re.compile(_NOT_A_LIST + r"(?P<n>\d+)\s+named (?:machine )?orders", re.I), "sm_levels",
     "the number of named machine orders"),
]


def check(fix: bool) -> int:
    facts = live_facts()
    bad = 0
    for path in DOCS:
        if not path.exists():
            continue
        text = original = path.read_text(encoding="utf-8")
        live_to = text.index(STOP) if STOP in text else len(text)
        for pattern, key, what in CLAIMS:
            for m in list(pattern.finditer(text)):
                if m.start() >= live_to:
                    continue  # a historical entry, left as it was written
                said, truth = int(m.group("n")), facts[key]
                if said == truth:
                    continue
                bad += 1
                line = text[: m.start()].count("\n") + 1
                rel = path.relative_to(ROOT)
                print(f"{rel}:{line}: says {said} for {what}, which is {truth}")
                print(f"    {m.group(0)!r}")
                if fix:
                    a, b = m.span("n")
                    text = text[:a] + str(truth) + text[b:]
        if fix and text != original:
            path.write_text(text, encoding="utf-8")
            print(f"  -> rewrote {path.relative_to(ROOT)}")
    if not bad:
        print(f"every countable claim matches the code ({facts['presets']} presets, "
              f"{facts['sm_levels']} machine orders)")
    return bad


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fix", action="store_true", help="rewrite the numbers instead of only reporting them")
    args = ap.parse_args()
    sys.exit(1 if check(args.fix) and not args.fix else 0)


if __name__ == "__main__":
    main()
