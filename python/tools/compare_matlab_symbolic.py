"""Check every component's linearisation against the MATLAB toolbox's own.

``matlab/Symbolic/A_*.txt`` holds the state matrices the original toolbox
generated, entry by entry, as symbolic expressions in the same parameter and
equilibrium symbols this port uses. That makes them a reference this port can
be checked against without MATLAB, and an *independent* one: they were
produced by different code, in a different language, from equations written
out by hand a second time.

The comparison is numeric but the claim is symbolic. ``A = Fx - Fz Gz^-1 Gx``
is a formula, and both tools compute the same formula; evaluating both at the
same randomly chosen point tests the whole expression rather than one
operating point, and a disagreement anywhere in it shows up. The point does
not have to be a physical equilibrium -- only a place where ``Gz`` inverts.

Run: ``python tools/compare_matlab_symbolic.py``
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from g2elin_core.components.gfl import gfl_dae  # noqa: E402
from g2elin_core.components.gfm import gfm_dae  # noqa: E402
from g2elin_core.components.ib import ib_dae  # noqa: E402
from g2elin_core.components.line import line_dae  # noqa: E402
from g2elin_core.components.load import load_dae  # noqa: E402
from g2elin_core.components.node import node_dae  # noqa: E402
from g2elin_core.components.sm import sm_dae  # noqa: E402

SYMBOLIC = Path(__file__).resolve().parents[2] / "matlab" / "Symbolic"
ENTRY = re.compile(r"A\[(\d+)\]\[(\d+)\]\s*=\s*(.+?)\s*$")
IDENTIFIER = re.compile(r"[A-Za-z_][A-Za-z_0-9]*")

# MATLAB writes `id`; Python can't, so this port calls it `id_`. Same variable.
RENAME = {"id": "id_"}
# Names in these expressions that are functions, not variables -- sympy's own
# must be used for them, or `cos(x)` parses as a symbol called cos applied to
# something.
FUNCTIONS = {"cos", "sin", "tan", "sqrt", "exp", "log", "abs", "atan", "acos", "asin"}


def parse_matlab_a(path: Path) -> dict[tuple[int, int], str]:
    """``{(row, col): expression}``, 1-based as MATLAB wrote it. Entries the
    file omits are zero -- the generator only prints the non-zero ones."""
    out: dict[tuple[int, int], str] = {}
    for line in path.read_text().splitlines():
        m = ENTRY.match(line.strip())
        if m:
            out[(int(m.group(1)), int(m.group(2)))] = m.group(3)
    return out


def matlab_expression(text: str, dae) -> sp.Expr:
    """One MATLAB entry as a sympy expression in *this port's* symbols.

    MATLAB substituted every state, algebraic variable and input by its
    ``X_0`` equilibrium symbol before printing; this port keeps the bare
    symbol and substitutes numbers later, so the two are the same quantity
    under different names.
    """
    # Every identifier is declared a symbol up front. Without that, sympify
    # resolves bare names against Python builtins and sympy's own namespace,
    # where `id` is a function, `E` is Euler's number and `I` is the
    # imaginary unit -- all of which appear here as variable names.
    text = text.replace("^", "**")
    local = {n: sp.Symbol(n) for n in set(IDENTIFIER.findall(text)) - FUNCTIONS}
    expr = sp.sympify(text, locals=local)
    swap = {}
    for sym in (*dae.state_syms, *dae.alg_syms, *dae.input_syms):
        name = sym.name
        matlab_name = next((k for k, v in RENAME.items() if v == name), name)
        swap[sp.Symbol(f"{matlab_name}_0")] = sym
    return expr.subs(swap)


def random_point(symbols, rng) -> dict:
    """A point to evaluate both sides at.

    Positive and O(1): the parameters are per-unit impedances, gains and time
    constants, and a negative inductance or a zero time constant would make
    the expressions singular for reasons that say nothing about either tool.
    """
    return {s: float(rng.uniform(0.3, 2.0)) for s in symbols}


def compare(name: str, dae, filename: str, rng, trials: int = 5) -> dict:
    path = SYMBOLIC / filename
    if not path.exists():
        return {"name": name, "status": "no reference file", "file": filename}

    entries = parse_matlab_a(path)
    n = len(dae.state_syms)
    rows = max(r for r, _ in entries)
    cols = max(c for _, c in entries)
    if rows > n or cols > n:
        return {"name": name, "status": f"size mismatch: reference is {rows}x{cols}, this port {n}x{n}",
                "file": filename}

    reference = sp.zeros(n, n)
    for (i, j), text in entries.items():
        reference[i - 1, j - 1] = matlab_expression(text, dae)

    # Everything linearize() will need a number for: its own parameters, and
    # the states/algebraics/inputs it substitutes an operating point into.
    # Missing even one makes _to_numpy refuse the whole matrix.
    all_syms = sorted(
        set(dae.param_syms()) | set(dae.state_syms) | set(dae.alg_syms)
        | set(dae.input_syms) | reference.free_symbols,
        key=lambda s: s.name,
    )

    worst = 0.0
    worst_cell = None
    checked = 0
    failures: list[str] = []
    for _ in range(trials):
        subs = random_point(all_syms, rng)
        try:
            mine = dae.linearize(subs).A
        except Exception as exc:  # a random point can land on a singular Gz
            failures.append(f"{type(exc).__name__}: {exc}")
            continue
        theirs = np.array(reference.subs(subs).evalf(), dtype=np.float64)
        scale = max(np.abs(theirs).max(), np.abs(mine).max(), 1.0)
        diff = np.abs(mine - theirs) / scale
        if diff.max() > worst:
            worst = float(diff.max())
            worst_cell = np.unravel_index(int(np.argmax(diff)), diff.shape)
        checked += 1

    return {
        "name": name, "file": filename, "n_states": n, "entries": len(entries),
        "points": checked, "worst_relative_error": worst,
        "worst_cell": None if worst_cell is None else (int(worst_cell[0]) + 1, int(worst_cell[1]) + 1),
        "status": "match" if checked and worst < 1e-9 else ("differs" if checked else "no usable point"),
        "why": failures[0] if failures and not checked else "",
    }


def best_variant(name: str, variants: dict, filename: str, rng) -> dict:
    """Compare every variant of a component and keep the one that matches.

    The toolbox generated one file per component, for whichever slack setting
    that component happened to have in the case it was generated from -- the
    machine as the slack unit, the converters as ordinary ones. Nothing in
    the file records which, so the match itself identifies it: a wrong
    variant differs in whole columns, not in the last digit.
    """
    tried = {label: compare(f"{name} [{label}]", dae, filename, rng) for label, dae in variants.items()}
    matched = [(label, r) for label, r in tried.items() if r["status"] == "match"]
    if len(matched) == 1:
        label, r = matched[0]
        r["name"] = name
        r["variant"] = label
        r["other_variants"] = {k: v["worst_relative_error"] for k, v in tried.items() if k != label}
        return r
    # Either none matched (a real disagreement) or both did (the variants are
    # identical here, which is fine); report the first either way.
    label, r = next(iter(tried.items()))
    r["name"] = name
    r["variant"] = label if len(matched) != len(tried) else "either"
    return r


def main() -> int:
    rng = np.random.default_rng(20260923)
    results = [
        best_variant("Synchronous machine",
                     {"slack": sm_dae(True), "non-slack": sm_dae(False)}, "A_SG0.txt", rng),
        best_variant("Grid-forming converter (droop)", {"non-slack": gfm_dae()}, "A_GFM0_Droop.txt", rng),
        best_variant("Grid-following converter", {"non-slack": gfl_dae()}, "A_GFL0.txt", rng),
        best_variant("Infinite bus", {"slack": ib_dae(True)}, "A_IB0.txt", rng),
        best_variant("Line", {"-": line_dae()}, "A_Line0.txt", rng),
        best_variant("Load", {"-": load_dae()}, "A_Load0.txt", rng),
        best_variant("Node", {"-": node_dae()}, "A_Node0.txt", rng),
    ]

    width = max(len(r["name"]) for r in results)
    print(f"{'component'.ljust(width)}  variant    states  entries  points  worst rel. err  verdict")
    for r in results:
        print(f"{r['name'].ljust(width)}  {r.get('variant', '-'):<9}  {r.get('n_states', '-'):>6}  "
              f"{r.get('entries', '-'):>7}  {r.get('points', '-'):>6}  "
              f"{r.get('worst_relative_error', float('nan')):>14.3e}  {r['status']}"
              + (f"  ({r['why'][:60]})" if r.get("why") else ""))
    print()
    print(f"{sum(r['status'] == 'match' for r in results)} of {len(results)} components match the "
          "MATLAB toolbox's own state matrices.")
    return 0 if all(r["status"] == "match" for r in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
