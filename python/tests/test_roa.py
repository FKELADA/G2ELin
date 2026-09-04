"""Lyapunov region-of-attraction tracing tests — feature 3.2.

Deliberately small (3x3 grid, ~100s): this is new infrastructure (see
stability/roa.py's module docstring) whose per-point cost is a multi-second
nonlinear trajectory, so the bar is "the mechanism produces physically
sensible classifications on a case we can reason about by hand" rather
than a finely resolved boundary. WSCC-9 is well-damped, so a moderate
(theta, omega) perturbation on a non-slack machine should mostly recover.
"""

from __future__ import annotations

import numpy as np
import pytest

from g2elin_core.network.presets import wscc9_3sm
from g2elin_core.powerflow import run_power_flow
from g2elin_core.stability import RoaAxis, find_state_index, trace_roa_grid
from g2elin_core.timedomain import build_nonlinear_network


@pytest.fixture(scope="module")
def model():
    network = wscc9_3sm()
    result = run_power_flow(network)
    return build_nonlinear_network(network, result)


@pytest.fixture(scope="module")
def grid_result(model):
    idx_w = find_state_index(model, "dw_r_{SM_2}")
    idx_th = find_state_index(model, "theta_{SM_2}")
    axis_x = RoaAxis(label="d(theta) [rad]", state_index=idx_th, offsets=np.linspace(-0.4, 0.4, 3))
    axis_y = RoaAxis(label="d(omega) [pu]", state_index=idx_w, offsets=np.linspace(-0.02, 0.02, 3))
    return trace_roa_grid(
        model, axis_x=axis_x, axis_y=axis_y, t_final=1.5, t_early=0.1,
        simulate_kwargs=dict(rtol=1e-3, atol=1e-5, first_step=1e-8),
    )


def test_find_state_index_exact_and_ambiguous(model):
    idx = find_state_index(model, "dw_r_{SM_2}")
    assert model.state_names[idx] == "dw_r_{SM_2}"
    with pytest.raises(ValueError, match="more than one"):
        find_state_index(model, "dw_r_")  # matches all 3 SMs
    with pytest.raises(ValueError, match="no state name"):
        find_state_index(model, "not_a_real_state")


def test_zero_perturbation_matches_baseline(grid_result):
    # The center grid point (offset 0, 0) is exactly the baseline trajectory
    # against itself -- distance should be ~0 at both checkpoints, trivially
    # "in the ROA". This is the sanity check that the distance metric and
    # baseline-vs-perturbed comparison are wired up correctly, independent
    # of any judgment call about what counts as "recovering".
    cy = cx = 1  # middle of a 3-wide grid
    assert not grid_result.failed[cy, cx]
    assert grid_result.early_distance[cy, cx] == pytest.approx(0.0, abs=1e-9)
    assert grid_result.late_distance[cy, cx] == pytest.approx(0.0, abs=1e-9)
    assert grid_result.in_roa[cy, cx]


def test_most_moderate_perturbations_recover(grid_result):
    # WSCC-9 is well-damped (see the golden modal-analysis tests: 12-49%
    # damping on the electromechanical modes). A a +-0.4 rad / +-0.02 pu
    # perturbation on one non-slack machine is a real but moderate swing,
    # not a bolted fault -- most of the 8 non-center grid points should
    # show the distance-to-baseline shrinking from t_early to t_final.
    valid = ~grid_result.failed
    assert valid.sum() >= 7  # at most the single largest corner may fail to solve (see below)
    assert grid_result.in_roa[valid].sum() >= valid.sum() - 1


def test_solver_failure_is_tracked_not_misclassified(grid_result):
    """A Newton solve failing (seen during development at the single most
    extreme corner: theta and omega both perturbed to their most negative
    offset together) is a distinct outcome from "diverged", and must not
    be silently folded into `in_roa=False` as if it were a real physical
    classification -- that would conflate "we don't know" with "no". Every
    failed point's distances are NaN and it's excluded from the recovery
    check above; this just confirms that bookkeeping.
    """
    for iy in range(grid_result.failed.shape[0]):
        for ix in range(grid_result.failed.shape[1]):
            if grid_result.failed[iy, ix]:
                assert np.isnan(grid_result.early_distance[iy, ix])
                assert np.isnan(grid_result.late_distance[iy, ix])
