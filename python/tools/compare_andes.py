"""Check the electromechanical modes against ANDES, an independent tool.

ANDES (https://github.com/CURENT/andes) is an open-source power system
simulator with its own eigenvalue analysis, its own machine models in
standard time-constant form (GENROU), and its own Kundur two-area case. None
of its code shares an ancestor with this one, which is what makes it a real
second opinion rather than a second run.

**Most of the work in a comparison like this is making the two cases the
same case.** Two tools disagreeing usually means they were asked different
questions, and every difference below was found by looking rather than
assumed:

- ANDES's ``kundur_full`` carries the book's inertias and reactances but a
  *different network*: three circuits on the 7-8 tie where the book has two,
  and no shunt capacitor banks at the load buses.
- Its machines are driven by EXDC2 and TGOV1. This tool's own Kundur preset
  uses a fast exciter with no stabiliser, which is what Example 12.6 is
  about -- a different machine to drive, not a different machine.

So this script lines the two up step by step and prints what each step is
worth. The point is not a single number but which differences matter.

Needs ``pip install andes``. Run: ``python tools/compare_andes.py``
"""

from __future__ import annotations

import math
import os
import sys
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from g2elin_core import reduction as R  # noqa: E402
from g2elin_core.modal import analyze, classify_modes  # noqa: E402
from g2elin_core.network.presets import kundur_two_area_classic  # noqa: E402
from g2elin_core.network.schema import ModelOptions  # noqa: E402
from g2elin_core.pipeline import linearize_network  # noqa: E402
from g2elin_core.powerflow import run_power_flow  # noqa: E402

# ANDES's EXDC2 as its Kundur case sets it. This tool's AVR is the same IEEE
# DC1A structure -- transducer, amplifier, exciter, rate feedback -- so the
# parameters carry across one for one.
EXDC2 = {"Tr": 0.02, "Ka": 20.0, "Ta": 0.02, "Ke": 1.0, "Te": 0.83, "Kfd": 0.0, "Tfd": 1.246}
# TGOV1's droop and first time constant. Its lead-lag (T2 = 2.1, T3 = 7.0) has
# no counterpart in this tool's first-order governor, which is the one piece
# of the comparison that stays approximate.
TGOV1 = {"TG": 0.49, "mp": 0.05}

BAND = (0.3, 2.0)  # the electromechanical band for this network


def electromechanical(system, modal) -> list[tuple[float, float]]:
    kinds = classify_modes(modal)
    out = []
    for j, lam in enumerate(modal.eigenvalues):
        if lam.imag <= 1e-8:
            continue
        hz = abs(lam) / (2 * math.pi)
        if not (BAND[0] <= hz <= BAND[1]) or kinds[j].category != R.SYNCHRONISATION:
            continue
        out.append((hz, -lam.real / abs(lam) * 100))
    return sorted(out)


def g2elin_modes(*, exciter: bool, extra_tie: int, drop_shunts: bool) -> list[tuple[float, float]]:
    net = kundur_two_area_classic()
    if exciter:
        net = net.model_copy(update={"der_units": [
            d.model_copy(update={"params": {**d.params, **EXDC2, **TGOV1}}) for d in net.der_units
        ]})
    update = {"models": ModelOptions(network_level="quasi_stationary", sm_level="order6")}
    if extra_tie:
        tie = next(ln for ln in net.lines if (ln.from_bus, ln.to_bus) == (7, 8))
        update["lines"] = list(net.lines) + [tie.model_copy() for _ in range(extra_tie)]
    if drop_shunts:
        update["shunts"] = []
    net = net.model_copy(update=update)
    system = linearize_network(net, run_power_flow(net))
    return electromechanical(system, analyze(system.A, system.state_names))


def andes_modes() -> list[tuple[float, float]]:
    import andes

    andes.config_logger(stream_level=40)
    case = Path(andes.__file__).parent / "cases" / "kundur" / "kundur_full.xlsx"
    ss = andes.run(str(case), routine="eig", no_output=True, default_config=True)
    return sorted(
        (abs(l) / (2 * np.pi), -l.real / abs(l) * 100)
        for l in np.asarray(ss.EIG.mu)
        if l.imag > 1e-8 and BAND[0] <= abs(l) / (2 * np.pi) <= BAND[1]
    )


def row(label: str, modes: list[tuple[float, float]], reference=None) -> str:
    cells = "  ".join(f"{hz:7.4f} Hz {d:6.2f}%" for hz, d in modes[:3])
    line = f"{label:44s}  {cells}"
    if reference and len(modes) >= 3:
        diff = "  ".join(
            f"{(m[0] - r[0]) / r[0] * 100:+6.2f}% {m[1] - r[1]:+6.2f}pp"
            for m, r in zip(modes[:3], reference[:3])
        )
        line += f"\n{'':44s}  {diff}"
    return line


def main() -> int:
    try:
        reference = andes_modes()
    except ImportError:
        print("ANDES is not installed -- pip install andes")
        return 2

    print("Kundur two-area, electromechanical modes (one inter-area, two local)\n")
    print(row("ANDES kundur_full (GENROU + EXDC2 + TGOV1)", reference))
    print()
    for label, kwargs in [
        ("this tool, as the preset ships", dict(exciter=False, extra_tie=0, drop_shunts=False)),
        ("  + ANDES's exciter and governor", dict(exciter=True, extra_tie=0, drop_shunts=False)),
        ("  + ANDES's third tie circuit", dict(exciter=True, extra_tie=1, drop_shunts=False)),
        ("  + ANDES's missing shunt banks", dict(exciter=True, extra_tie=1, drop_shunts=True)),
    ]:
        print(row(label, g2elin_modes(**kwargs), reference))
    print(
        "\nThe exciter governs the damping, the tie circuit governs the inter-area\n"
        "frequency, and the shunt banks pull it back the other way -- which is\n"
        "why they have to be matched before any of this means anything."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
