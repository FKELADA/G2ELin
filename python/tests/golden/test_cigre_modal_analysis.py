"""Linear small-signal modal analysis regression test for the CIGRE islanded
1SM+2GFM+1GFL preset. Same status as the other golden tests: no MATLAB
reference yet, so this checks structural correctness. This case is the
first to exercise GFM and GFL in the interconnection (WSCC-9/3SM is SM-only).
"""

from __future__ import annotations

import numpy as np
import pytest

from g2elin_core.network.presets import cigre_islanded_1sm_2gfm_1gfl
from g2elin_core.pipeline import linearize_network
from g2elin_core.powerflow import run_power_flow
from g2elin_core.modal import analyze


@pytest.fixture(scope="module")
def system():
    network = cigre_islanded_1sm_2gfm_1gfl()
    result = run_power_flow(network)
    assert result.converged
    return linearize_network(network, result)


@pytest.fixture(scope="module")
def modal(system):
    return analyze(system.A, system.state_names)


def test_state_count(system):
    # 1 SM x 19 (PSS on) + 2 GFM x 15 (Droop) + 1 GFL x 14
    # + 14 nodes x 2 + 13 lines x 2 + 13 loads x 2
    # + 1 for the reference frame's own angle (components/frame.py)
    expected = 19 + 2 * 15 + 14 + 14 * 2 + 13 * 2 + 13 * 2 + 1
    assert system.A.shape == (expected, expected)


def test_no_nan_or_inf(system):
    import numpy as np

    assert np.isfinite(system.A).all()
    assert np.isfinite(system.B).all()
    assert np.isfinite(system.C).all()
    assert np.isfinite(system.D).all()


def test_system_is_stable(modal):
    # Excluding the reference-angle modes, which are the model's own free
    # coordinates and sit on the imaginary axis by construction (they land a
    # hair either side of it numerically) -- see modal.reference_angle_modes.
    from g2elin_core.modal import reference_angle_modes

    ref = set(reference_angle_modes(modal))
    physical = [z for j, z in enumerate(modal.eigenvalues) if j not in ref]
    assert len(ref) == 2 and (np.array(physical).real < 1e-6).all()


def test_participation_columns_sum_to_one(modal):
    import numpy as np

    sums = modal.participation.sum(axis=0)
    nonzero = sums[sums > 0]
    assert np.allclose(nonzero, 1.0, atol=1e-6)
