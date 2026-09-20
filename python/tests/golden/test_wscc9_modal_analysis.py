"""Linear small-signal modal analysis regression test for WSCC 9-bus / 3-SM.

Same status as the power-flow golden tests: no MATLAB-exported reference
yet, so this checks structural correctness (state count, stability, the
presence of an electromechanical mode in a physically sane frequency range)
rather than exact eigenvalues.
"""

from __future__ import annotations

import pytest

from g2elin_core.network.presets import wscc9_3sm
from g2elin_core.pipeline import linearize_network
from g2elin_core.powerflow import run_power_flow
from g2elin_core.modal import analyze


@pytest.fixture(scope="module")
def system():
    network = wscc9_3sm()
    result = run_power_flow(network)
    assert result.converged
    return linearize_network(network, result)


@pytest.fixture(scope="module")
def modal(system):
    return analyze(system.A, system.state_names)


def test_state_count(system):
    # 3 SM x 19 states (PSS on) + 6 nodes x 2 + 6 lines x 2 + 3 loads x 2,
    # + 1 for the reference frame's own angle (components/frame.py)
    expected = 3 * 19 + 6 * 2 + 6 * 2 + 3 * 2 + 1
    assert system.A.shape == (expected, expected)


def test_no_nan_or_inf(system):
    import numpy as np

    assert np.isfinite(system.A).all()
    assert np.isfinite(system.B).all()
    assert np.isfinite(system.C).all()
    assert np.isfinite(system.D).all()


def test_system_is_stable(modal):
    # A well-damped power system operating near nominal should have no
    # right-half-plane eigenvalues (allow a small numerical tolerance).
    assert (modal.eigenvalues.real < 1e-6).all()


def test_has_electromechanical_modes(modal):
    # 2 synchronous machines swinging against each other/the slack should
    # produce complex mode pairs in the classic 0.1-3 Hz electromechanical
    # range (Kundur), clearly separated from the ~50-60 Hz electrical modes.
    complex_modes = modal.eigenvalues[modal.eigenvalues.imag > 1e-6]
    freqs_hz = complex_modes.imag / (2 * 3.141592653589793)
    assert ((freqs_hz > 0.1) & (freqs_hz < 3.0)).any()


def test_participation_columns_sum_to_one(modal):
    import numpy as np

    sums = modal.participation.sum(axis=0)
    nonzero = sums[sums > 0]
    assert np.allclose(nonzero, 1.0, atol=1e-6)
