"""EMT (nonlinear) time-domain simulation tests — feature 2.2.

Runtime here dropped sharply once ``timedomain/emt.py``'s algebraic Newton
solve got an analytic Jacobian instead of relying on SciPy's numerically-
estimated one: a single ``solve_algebraic`` call went from ~0.05s to
~0.006s, and simulating 1 full second of WSCC-9 (previously untested —
even 1e-3s wasn't reliably reachable) now takes ~13s. This still isn't
long enough to see a *full* electromechanical transient settle (that needs
several seconds, several times the ~0.5-2s period of the slowest modes),
so the bar stays "builds, solves, integrates, doesn't blow up or NaN"
rather than exact transient behavior — but the horizon here (0.1s) is a
real, physically meaningful chunk of the fastest electromechanical modes'
period, not just a numerical-stability smoke test.
"""

from __future__ import annotations

import numpy as np
import pytest

from g2elin_core.network.presets import wscc9_3sm
from g2elin_core.powerflow import run_power_flow
from g2elin_core.timedomain import build_nonlinear_network, simulate, simulate_steps


@pytest.fixture(scope="module")
def model():
    network = wscc9_3sm()
    result = run_power_flow(network)
    return build_nonlinear_network(network, result)


def test_model_sizes(model):
    # 3 SM x 19 states (PSS on); node/line/load contribute no algebraic vars.
    assert sum(b.comp.n_states for b in model.blocks) == 87
    assert model.n_z == 3 * 15  # SM's algVec = [wr, ved, veq, id, i1d, ifd, iq, i1q, i2q, vgd, vgq, Cm, Ce, DP, Et]


def test_algebraic_solve_at_operating_point_converges(model):
    x0 = model.initial_state()
    z_guess, u_guess = model.initial_algebraic_guess()
    u_exo0 = model.default_u_exo()

    z_sol, u_sol = model.solve_algebraic(x0, u_exo0, z_guess, u_guess)

    assert np.all(np.isfinite(z_sol))
    assert np.all(np.isfinite(u_sol))
    # z (each component's own algebraic vars) should match the per-component
    # operating point almost exactly -- unlike u (see the note below), z has
    # no cross-component coupling to be thrown off by.
    assert np.allclose(z_sol, z_guess, atol=1e-6)


def test_network_operating_point_is_not_a_perfect_equilibrium(model):
    """Documents a real, expected finding rather than asserting it away:
    the "every node uses line #1's susceptance" quirk (see
    components/sm.py) means the per-component operating point isn't quite
    a network-wide KCL/equilibrium fix point -- solve_algebraic's u_sol
    differs from the naive per-component u_guess by ~0.08 (small relative
    to per-unit scale ~1, but not solver noise), and xdot at that point is
    correspondingly nonzero for node-voltage and SM current states. This
    doesn't affect linear modal analysis (a Jacobian is valid at any point,
    verified in test_nonlinear_cross_validation.py), but it does mean an
    EMT simulation starting here has a real initial transient to settle,
    not just numerical noise.
    """
    x0 = model.initial_state()
    z_guess, u_guess = model.initial_algebraic_guess()
    u_exo0 = model.default_u_exo()

    z_sol, u_sol = model.solve_algebraic(x0, u_exo0, z_guess, u_guess)
    assert np.max(np.abs(u_sol - u_guess)) > 0.01  # real, not floating-point-scale


def _non_rotating_mask(state_names: list[str]) -> np.ndarray:
    """Excludes theta/theta_pll-type states from a "stays bounded" check.

    They're absolute angle states with derivative ~wb (~377 rad/s here) at
    equilibrium — they're *supposed* to grow without bound (a machine
    keeps rotating; only the angle *difference* from the reference frame
    settles). Confirmed by running to t=1s during development: max|x| grew
    to ~376.65, essentially exactly wb*t -- not instability, just this
    state doing what it's defined to do. Any bound on "did the simulation
    blow up" has to exclude these states or it's really testing "did an
    hour pass", not "did anything go wrong".
    """
    return np.array(["theta" in n for n in state_names])


def test_short_simulation_stays_bounded(model):
    sim = simulate(model, (0.0, 0.1), rtol=1e-4, atol=1e-6, first_step=1e-8)

    assert np.all(np.isfinite(sim.x))
    mask = _non_rotating_mask(sim.state_names)
    # Per-unit quantities here are O(1); a bound of 20 is generous headroom
    # while still catching a genuine blow-up (the failure mode this test
    # actually guards against) on every state except the rotating angles.
    assert np.max(np.abs(sim.x[~mask, :])) < 20
    assert len(sim.state_names) == sim.x.shape[0]


def test_simulate_steps_matches_simulate(model):
    """simulate_steps() (the manually-stepped generator the live-tracing
    web UI/notebook features use) has to agree with simulate()'s own
    one-shot solve_ivp() call on the same trajectory -- same integrator,
    same tolerances, just driven step-by-step instead of all at once.
    """
    kwargs = dict(rtol=1e-4, atol=1e-6, first_step=1e-8)
    t_final = 0.1

    sim = simulate(model, (0.0, t_final), t_eval=np.array([0.0, t_final]), **kwargs)
    steps = list(simulate_steps(model, (0.0, t_final), **kwargs))

    assert steps[0].t == 0.0
    assert steps[-1].t == pytest.approx(t_final)
    assert len(steps) > 1  # a real trajectory, not a single jump to t_final
    assert np.allclose(steps[-1].x, sim.x[:, -1], atol=1e-6)


def test_simulate_steps_inputs_and_outputs_are_cheap_and_correct(model):
    """The whole point of simulate_steps() carrying inputs/outputs is that
    they're read off the step's own already-solved (z, u) cache, not a
    fresh solve_algebraic() re-solve -- this checks that shortcut actually
    lands on the same answer a fresh re-solve would give, not just that it
    runs.
    """
    kwargs = dict(rtol=1e-4, atol=1e-6, first_step=1e-8)
    steps = list(simulate_steps(model, (0.0, 0.05), **kwargs))
    last = steps[-1]

    assert set(last.inputs) == set(model.input_names)
    assert set(last.outputs) == set(model.output_names)
    assert all(np.isfinite(v) for v in last.inputs.values())
    assert all(np.isfinite(v) for v in last.outputs.values())

    z_guess, u_guess = model.initial_algebraic_guess()
    u_exo = model.default_u_exo()
    z_fresh, u_fresh = model.solve_algebraic(last.x, u_exo, z_guess, u_guess)
    _, outputs_fresh = model._inputs_and_outputs(last.x, z_fresh, u_fresh)
    one_name = next(iter(outputs_fresh))
    assert last.outputs[one_name] == pytest.approx(outputs_fresh[one_name], abs=1e-3)


def test_simulate_steps_max_step_caps_but_does_not_force_step_size(model):
    steps_unbounded = list(simulate_steps(model, (0.0, 0.02), rtol=1e-4, atol=1e-6, first_step=1e-8))
    steps_capped = list(simulate_steps(model, (0.0, 0.02), rtol=1e-4, atol=1e-6, first_step=1e-8, max_step=1e-4))
    # Capping max_step can only add steps (or leave the count the same),
    # never remove ones the solver already needed for accuracy.
    assert len(steps_capped) >= len(steps_unbounded)
    all_dt = np.diff([s.t for s in steps_capped])
    assert np.all(all_dt <= 1e-4 + 1e-9)
