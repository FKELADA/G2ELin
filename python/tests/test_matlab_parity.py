"""The component linearisations, against the MATLAB toolbox's own.

``matlab/Symbolic/A_*.txt`` are the state matrices the original toolbox
generated. They are an independent reference -- different code, different
language, the equations written out a second time -- and they were sitting
in the repository unused.

This is the check that the port is faithful. It says nothing about whether
the equations are right (both tools implement the same physics from the same
source); for that see docs/validation/cross-validation-2026-09.md.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

TOOLS = Path(__file__).resolve().parents[1] / "tools"
sys.path.insert(0, str(TOOLS))

from compare_matlab_symbolic import SYMBOLIC, best_variant  # noqa: E402

import numpy as np  # noqa: E402

from g2elin_core.components.gfl import gfl_dae  # noqa: E402
from g2elin_core.components.gfm import gfm_dae  # noqa: E402
from g2elin_core.components.ib import ib_dae  # noqa: E402
from g2elin_core.components.line import line_dae  # noqa: E402
from g2elin_core.components.load import load_dae  # noqa: E402
from g2elin_core.components.node import node_dae  # noqa: E402
from g2elin_core.components.sm import sm_dae  # noqa: E402

pytestmark = pytest.mark.skipif(
    not SYMBOLIC.exists(),
    reason="the MATLAB toolbox's symbolic output isn't present in this checkout",
)

CASES = [
    ("synchronous machine", {"slack": sm_dae(True), "non-slack": sm_dae(False)}, "A_SG0.txt"),
    ("grid-forming converter", {"non-slack": gfm_dae()}, "A_GFM0_Droop.txt"),
    ("grid-following converter", {"non-slack": gfl_dae()}, "A_GFL0.txt"),
    ("infinite bus", {"slack": ib_dae(True)}, "A_IB0.txt"),
    ("line", {"-": line_dae()}, "A_Line0.txt"),
    ("load", {"-": load_dae()}, "A_Load0.txt"),
    ("node", {"-": node_dae()}, "A_Node0.txt"),
]


@pytest.mark.parametrize("name,variants,filename", CASES, ids=[c[0] for c in CASES])
def test_the_state_matrix_matches_the_matlab_toolbox(name, variants, filename):
    """Both tools compute A = Fx - Fz Gz^-1 Gx, which is a formula; evaluating
    both at the same random point tests the whole expression rather than one
    operating point."""
    rng = np.random.default_rng(20260923)
    result = best_variant(name, variants, filename, rng)
    assert result["status"] == "match", (
        f"{name}: worst relative error {result.get('worst_relative_error')} "
        f"at entry {result.get('worst_cell')} of {filename}"
    )
    assert result["points"] > 0, "no point was usable -- the comparison proved nothing"
    assert result["worst_relative_error"] < 1e-12
