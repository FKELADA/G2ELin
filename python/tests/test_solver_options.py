"""The time-domain solver: the analytic Jacobian, the method choice, and
fixed stepping.

The Jacobian tests matter most. An implicit solver given no ``jac``
finite-differences one, at a full coupled Newton solve per state -- a cost
SciPy doesn't even report in ``nfev``. Getting that matrix wrong would not
raise; it would quietly make the solver take smaller steps, or converge to
something else.
"""

from __future__ import annotations

import time

import numpy as np
import pytest

from g2elin_core.network.presets import wscc9_3sm
from g2elin_core.network.schema import ModelOptions, Network
from g2elin_core.powerflow import run_power_flow
from g2elin_core.timedomain import (
    DEFAULT_SOLVER, SOLVERS, build_nonlinear_network, simulate, simulate_fixed_step,
)


def model_for(**models) -> tuple:
    net = wscc9_3sm().model_copy(update={
        "nodes_share_first_line_b": False, "min_node_b_pu": 1e-4,
        "models": ModelOptions(**models),
    })
    model = build_nonlinear_network(net, run_power_flow(net))
    x0 = model.initial_state().copy()
    x0[model.state_names.index("dw_r_{SM_2}")] += 1e-3
    return model, x0


@pytest.fixture(scope="module")
def emt_model():
    return model_for()


@pytest.fixture(scope="module")
def rms_model():
    return model_for(network_level="quasi_stationary", sm_level="order6")


# --- the analytic Jacobian ----------------------------------------------------


@pytest.mark.parametrize("reduced", [False, True])
def test_the_ode_jacobian_matches_a_finite_difference(reduced):
    """d(xdot)/dx is assembled from the constraint derivatives rather than
    differenced. It has to agree with the difference it replaces."""
    model, x0 = model_for(network_level="quasi_stationary", sm_level="order6") if reduced else model_for()
    u = model.default_u_exo()
    x = x0 + 1e-4
    z, uu = model.solve_algebraic(x, u)
    analytic = model.ode_jacobian(x, u, z, uu)

    def f(xv):
        xdot, _, _ = model.rhs(xv, u, z, uu)
        return xdot

    n = len(x)
    difference = np.zeros((n, n))
    h = 1e-7
    for j in range(n):
        e = np.zeros(n)
        e[j] = h
        difference[:, j] = (f(x + e) - f(x - e)) / (2 * h)

    scale = np.abs(difference).max()
    assert np.abs(analytic - difference).max() / scale < 1e-6


def test_an_implicit_solver_is_actually_given_the_jacobian(emt_model):
    """SciPy's ``nfev`` counts only the integrator's own evaluations -- the
    ones its finite-difference Jacobian makes are invisible in it. So the
    two numbers being *equal* is exactly the property to assert: every RHS
    call is one the integrator asked for, none went into building a
    Jacobian. Measured before this existed: 257 real calls against an
    ``nfev`` of 81."""
    model, x0 = emt_model
    calls = {"n": 0}
    real_rhs = model.rhs

    def counting_rhs(x, u_exo, z=None, uu=None):
        calls["n"] += 1
        return real_rhs(x, u_exo, z, uu)

    model.rhs = counting_rhs
    try:
        result = simulate(model, (0.0, 0.2), x0=x0, method="Radau", rtol=1e-4, atol=1e-6)
    finally:
        model.rhs = real_rhs
    assert result.scipy_result.njev > 0, "the solver never asked for a Jacobian -- test proves nothing"
    assert calls["n"] == result.scipy_result.nfev
    assert result.t[-1] == pytest.approx(0.2)


# --- the solver catalogue -----------------------------------------------------


def test_the_catalogue_is_consistent():
    assert DEFAULT_SOLVER in SOLVERS
    for name, info in SOLVERS.items():
        assert {"implicit", "label", "note"} <= set(info)
        assert info["note"], f"{name} has no guidance"


@pytest.mark.parametrize("method", ["Radau", "BDF", "LSODA"])
def test_every_implicit_method_agrees_on_the_answer(method, rms_model):
    """Different integrators, same trajectory -- otherwise one of them is
    wrong, and offering the choice would be offering a coin flip."""
    model, x0 = rms_model
    reference = simulate(model, (0.0, 0.3), x0=x0, method="BDF", rtol=1e-8, atol=1e-10)
    got = simulate(model, (0.0, 0.3), x0=x0, method=method, rtol=1e-6, atol=1e-8)
    j = model.state_names.index("dw_r_{SM_1}")
    on_reference_grid = np.interp(reference.t, got.t, got.x[j])
    # Relative to the signal's own size: an absolute bound would just be a
    # statement about how big this particular disturbance happens to be.
    scale = np.abs(reference.x[j]).max()
    assert np.abs(on_reference_grid - reference.x[j]).max() / scale < 1e-2


def test_a_looser_tolerance_costs_fewer_steps(rms_model):
    model, x0 = rms_model
    tight = simulate(model, (0.0, 0.3), x0=x0, method="BDF", rtol=1e-8, atol=1e-10)
    loose = simulate(model, (0.0, 0.3), x0=x0, method="BDF", rtol=1e-3, atol=1e-5)
    assert len(loose.t) < len(tight.t)


# --- fixed stepping -----------------------------------------------------------


def test_fixed_stepping_lands_exactly_on_its_own_grid(rms_model):
    model, x0 = rms_model
    result = simulate_fixed_step(model, (0.0, 0.1), 1e-3, x0=x0)
    assert len(result.t) == 101
    np.testing.assert_allclose(np.diff(result.t), 1e-3, rtol=1e-9)


def test_fixed_stepping_converges_at_second_order(rms_model):
    """Trapezoidal is 2nd order, so halving the step should cut the error by
    about four. A first-order slip (an implicit-Euler bug) would show as two."""
    model, x0 = rms_model
    reference = simulate(model, (0.0, 0.2), x0=x0, method="BDF", rtol=1e-10, atol=1e-12)
    j = model.state_names.index("dw_r_{SM_1}")

    def error(step: float) -> float:
        got = simulate_fixed_step(model, (0.0, 0.2), step, x0=x0)
        return float(np.abs(np.interp(reference.t, got.t, got.x[j]) - reference.x[j]).max())

    coarse, fine = error(2e-3), error(1e-3)
    assert fine < coarse
    assert 2.5 < coarse / fine < 6.0, f"order looks like {np.log2(coarse / fine):.2f}, expected ~2"


def test_fixed_and_variable_stepping_agree(rms_model):
    model, x0 = rms_model
    variable = simulate(model, (0.0, 0.2), x0=x0, method="BDF", rtol=1e-8, atol=1e-10)
    fixed = simulate_fixed_step(model, (0.0, 0.2), 2e-4, x0=x0)
    j = model.state_names.index("dw_r_{SM_1}")
    on_grid = np.interp(variable.t, fixed.t, fixed.x[j])
    assert np.abs(on_grid - variable.x[j]).max() < 1e-7


def test_a_non_positive_fixed_step_is_refused(rms_model):
    model, x0 = rms_model
    with pytest.raises(ValueError, match="positive"):
        simulate_fixed_step(model, (0.0, 0.1), 0.0, x0=x0)
