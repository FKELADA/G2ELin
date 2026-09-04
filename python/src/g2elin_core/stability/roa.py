"""Simulation-based Lyapunov region-of-attraction (ROA) tracing — feature 3.2.

Perturbs the operating point along two chosen state directions, integrates
the *nonlinear* EMT model (:mod:`g2elin_core.timedomain`) briefly from each
perturbed point, and classifies whether the trajectory is trending back
toward the (unperturbed) equilibrium's own trajectory or away from it. This
is the "cheap simulation-based/trajectory-sampling ROA estimator" the
migration plan proposed as a fast, dependency-light method — as opposed to
a formal SOS/Lyapunov-function certification (e.g. via ``pydrake``), which
isn't implemented here.

This module exists *because* ``timedomain/emt.py`` got an analytic Newton
Jacobian: simulating enough multi-second trajectories for a grid sweep
wasn't practical on the numerically-estimated-Jacobian solver (~0.05s per
algebraic solve; a 1s trajectory wasn't validated). With the analytic one
(~0.006s per solve), a 1s WSCC-9 trajectory takes ~13s, which makes a
modest grid (order 25-50 points) a few-minutes job instead of hours.

**Why "trending" rather than "converged by a fixed time".** WSCC-9's
electromechanical modes have ~0.5-1.4s periods (0.7-2.2 Hz) with moderate
damping — a single fixed short horizon catches an early, still-oscillating
part of the transient for perfectly stable cases too, so "is the final
distance below some threshold" misclassifies stable points as unstable
just because they haven't finished settling yet (confirmed during
development: even a mild +-0.02 pu speed perturbation "failed" a
fixed-endpoint check at t=0.5s). Comparing the distance-to-baseline at two
points in time — early (after the fast sub-millisecond transients have
died out) vs. late — and checking whether it *shrank* is a much better
proxy for "is this converging" from a short, affordable trajectory, at the
cost of being qualitative (a decaying-but-still-large-amplitude oscillation
reads as "in the ROA" here even though it may take much longer than
``t_final`` to actually settle).

**The "theta problem".** A rotating machine's absolute angle state
(``theta``/``theta_pll``) grows without bound even at a stable operating
point — it's an absolute angle, not a deviation (see
``test_emt_simulation.py``'s ``_non_rotating_mask``). Any distance metric
used here has to exclude those states, or every classification would read
"diverging" including the truly stable ones. This module excludes them
(:func:`_non_rotating_mask`) and relies on the *other* physically coupled
states (fluxes, currents, mechanical power, AVR states, node voltages) to
reveal a loss of synchronism indirectly — a real simplification, not a
rigorous treatment of angle differences between machines, documented here
rather than silently assumed.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from g2elin_core.timedomain.emt import NonlinearNetworkModel, Vec, simulate

DEFAULT_SIMULATE_KWARGS = dict(rtol=1e-4, atol=1e-6, first_step=1e-8)


def _non_rotating_mask(state_names: list[str]) -> np.ndarray:
    return np.array(["theta" not in n for n in state_names])


def find_state_index(model: NonlinearNetworkModel, name_contains: str) -> int:
    """Locates a state by a (unique) substring of its name, e.g.
    ``"dw_r_{SM_2}"`` or just ``"SM_2}"`` if that's unambiguous. Raises if
    zero or more than one state matches, rather than silently guessing.
    """
    matches = [i for i, n in enumerate(model.state_names) if name_contains in n]
    if len(matches) == 0:
        raise ValueError(f"no state name contains {name_contains!r}")
    if len(matches) > 1:
        names = [model.state_names[i] for i in matches]
        raise ValueError(f"{name_contains!r} matches more than one state: {names}")
    return matches[0]


@dataclass
class RoaAxis:
    label: str
    state_index: int
    offsets: Vec  # perturbation values added to x0[state_index], e.g. np.linspace(-0.3, 0.3, 7)


@dataclass
class RoaGridResult:
    axis_x: RoaAxis
    axis_y: RoaAxis
    in_roa: np.ndarray  # bool, shape (len(y.offsets), len(x.offsets)) -- trending toward baseline
    early_distance: np.ndarray  # float, same shape
    late_distance: np.ndarray  # float, same shape (NaN where the solve failed)
    failed: np.ndarray  # bool, same shape
    baseline_trajectory: object  # EmtSimulationResult, the unperturbed reference trajectory
    t_eval: Vec


def trace_roa_grid(
    model: NonlinearNetworkModel,
    *,
    axis_x: RoaAxis,
    axis_y: RoaAxis,
    t_final: float,
    t_early: float,
    simulate_kwargs: dict | None = None,
) -> RoaGridResult:
    """Samples a grid of initial-condition perturbations along ``axis_x``/
    ``axis_y``; for each, integrates to ``t_final`` and classifies it as
    "in the region of attraction" if its distance to the unperturbed
    trajectory's own path (excluding rotating-angle states — see the
    module docstring) is *smaller* at ``t_final`` than at ``t_early``.

    ``t_early`` should be well past the sub-millisecond fast transients
    (order 1e-2 to 1e-1 s worked for WSCC-9's electromechanical modes) but
    well before ``t_final``, so the comparison actually reflects the slow
    dynamics this is meant to characterize rather than the initial
    algebraic-consistency transient every trajectory here starts with
    (see ``test_emt_simulation.py``'s equilibrium-isn't-perfect note).
    """
    if not (0 < t_early < t_final):
        raise ValueError(f"need 0 < t_early < t_final, got t_early={t_early}, t_final={t_final}")

    kwargs = {**DEFAULT_SIMULATE_KWARGS, **(simulate_kwargs or {})}
    x0 = model.initial_state()
    mask = _non_rotating_mask(model.state_names)
    t_eval = np.array([t_early, t_final])

    baseline = simulate(model, (0.0, t_final), x0=x0, t_eval=t_eval, **kwargs)

    ny, nx = len(axis_y.offsets), len(axis_x.offsets)
    in_roa = np.zeros((ny, nx), dtype=bool)
    early_distance = np.full((ny, nx), np.nan)
    late_distance = np.full((ny, nx), np.nan)
    failed = np.zeros((ny, nx), dtype=bool)

    for iy, dy in enumerate(axis_y.offsets):
        for ix, dx in enumerate(axis_x.offsets):
            xp = x0.copy()
            xp[axis_x.state_index] += dx
            xp[axis_y.state_index] += dy
            try:
                sim = simulate(model, (0.0, t_final), x0=xp, t_eval=t_eval, **kwargs)
                if not np.all(np.isfinite(sim.x)):
                    raise RuntimeError("non-finite state along trajectory")
                diff = (sim.x - baseline.x)[mask, :]
                d_early = float(np.linalg.norm(diff[:, 0]))
                d_late = float(np.linalg.norm(diff[:, 1]))
                early_distance[iy, ix] = d_early
                late_distance[iy, ix] = d_late
                # <=, not <: the zero-perturbation grid point has d_early ==
                # d_late == 0 exactly (it's the baseline compared to itself),
                # and a strict "did it shrink" reads that as failing to
                # recover -- "did not grow" is what's actually being asked.
                in_roa[iy, ix] = d_late <= d_early
            except Exception:
                failed[iy, ix] = True

    return RoaGridResult(
        axis_x=axis_x, axis_y=axis_y, in_roa=in_roa, early_distance=early_distance,
        late_distance=late_distance, failed=failed, baseline_trajectory=baseline, t_eval=t_eval,
    )
