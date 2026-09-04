"""Tests for modal/toolbox.py -- eigenvalue sensitivity, mode shape, free
motion response, and step response, ported from the parts of
Functions/modal_analysis.m that modal/analysis.py's own docstring flagged
as "not ported yet". All pure linear algebra on an already-validated
eigendecomposition (see test_wscc9_modal_analysis.py for that), so these
tests check the *new* math's own properties rather than re-deriving what
modal analysis already validates.
"""

from __future__ import annotations

import numpy as np
import pytest

from g2elin_core.modal import analyze, eigenvalue_sensitivity, free_response, mode_shape, step_response
from g2elin_core.network.presets import wscc9_3sm
from g2elin_core.pipeline import linearize_network
from g2elin_core.powerflow import run_power_flow


@pytest.fixture(scope="module")
def system_and_modal():
    network = wscc9_3sm()
    result = run_power_flow(network)
    system = linearize_network(network, result)
    modal = analyze(system.A, system.state_names)
    return system, modal


def test_sensitivity_masked_to_nonzero_A_entries(system_and_modal):
    system, modal = system_and_modal
    mode = len(modal.eigenvalues) // 2
    sens = eigenvalue_sensitivity(modal, system.A, mode)
    zero_mask = system.A == 0
    # Every A-zero entry must be masked out (this is the point of the mask,
    # not incidental -- see the module docstring on why).
    assert np.all(sens.matrix[zero_mask] == 0)
    assert len(sens.top) <= 8
    # The top entries must actually be the largest values in the matrix.
    flat_sorted = np.sort(sens.matrix, axis=None)[::-1]
    assert sens.top[0].value == pytest.approx(flat_sorted[0], rel=1e-9)


def test_mode_shape_top5_have_unit_magnitude_by_construction(system_and_modal):
    # The original tool's mode shape plot fixes magnitude at 1 and only
    # varies phase -- this just confirms the states/angles line up with the
    # actual top-5 participation-factor states for that mode (not that the
    # magnitude is 1, which is true by construction, not something to test).
    _, modal = system_and_modal
    mode = 10
    ms = mode_shape(modal, mode, n_top=5)
    assert len(ms.states) == 5
    part_col = np.abs(modal.participation[:, mode])
    expected_top5 = {modal.state_names[i] for i in np.argsort(part_col)[::-1][:5]}
    assert set(ms.states) == expected_top5
    assert all(-180.0 <= a <= 180.0 for a in ms.angles_deg)


def test_free_response_zero_offset_is_flat(system_and_modal):
    # A zero perturbation should give an exactly-zero response at every
    # state and time -- the trivial sanity check that modal expansion and
    # the equilibrium (deviation-form) convention agree.
    _, modal = system_and_modal
    t = np.linspace(0.0, 1.0, 20)
    x_t = free_response(modal, 0, 0.0, t)
    assert np.allclose(x_t, 0.0, atol=1e-9)


def test_free_response_matches_direct_state_space_integration(system_and_modal):
    # The whole point of modal expansion is that it's a closed-form
    # shortcut for x(t) = expm(A*t) @ x0 -- cross-check against that directly
    # (scipy.linalg.expm), not just "it runs and is finite".
    system, modal = system_and_modal
    idx = 5
    offset = 0.02
    t_final = 0.3
    x0 = np.zeros(system.A.shape[0])
    x0[idx] = offset

    import scipy.linalg

    x_direct = scipy.linalg.expm(system.A * t_final) @ x0
    t = np.array([0.0, t_final])
    x_modal = free_response(modal, idx, offset, t)
    assert np.allclose(x_modal[:, -1], x_direct, atol=1e-6, rtol=1e-4)


def test_step_response_matches_direct_state_space_step(system_and_modal):
    system, _ = system_and_modal
    import scipy.signal

    i = 0
    j = 0
    t = np.linspace(0.0, 0.5, 30)
    y = step_response(system, system.input_names[i], system.output_names[j], amplitude=1.0, t=t)

    sys = scipy.signal.StateSpace(
        system.A, system.B[:, i : i + 1], system.C[j : j + 1, :], system.D[j : j + 1, i : i + 1]
    )
    _, y_ref = scipy.signal.step(sys, T=t)
    assert np.allclose(y, y_ref, atol=1e-8)


def test_step_response_amplitude_scales_linearly(system_and_modal):
    system, _ = system_and_modal
    t = np.linspace(0.0, 0.3, 20)
    y1 = step_response(system, system.input_names[0], system.output_names[0], amplitude=1.0, t=t)
    y2 = step_response(system, system.input_names[0], system.output_names[0], amplitude=2.5, t=t)
    assert np.allclose(y2, 2.5 * y1, atol=1e-10)
