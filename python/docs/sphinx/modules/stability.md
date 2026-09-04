# `stability` — Lyapunov region-of-attraction tracing (simulation-based)

Feature 3.2 from the migration plan. Perturbs the operating point along two
chosen state directions, integrates the *nonlinear* EMT model
({doc}`timedomain <timedomain>`) briefly from each perturbed point, and
classifies whether the trajectory is trending back toward the (unperturbed)
equilibrium's own trajectory or away from it — the "cheap
simulation-based/trajectory-sampling ROA estimator" the migration plan
proposed as a fast, dependency-light method, as opposed to a formal
SOS/Lyapunov-function certification (e.g. via `pydrake`), which isn't
implemented.

This module exists *because* {doc}`timedomain <timedomain>` got an
analytic Newton Jacobian: simulating enough multi-second trajectories for a
grid sweep wasn't practical on the numerically-estimated-Jacobian solver.
With the analytic one, a 1s WSCC-9 trajectory takes ~13s, making a modest
grid (order 25–50 points) a few-minutes job instead of hours.

## The classification loop

```{mermaid}
flowchart TD
    X0["x0 = model.initial_state()"] --> BASE["baseline = simulate(x0, [t_early, t_final])"]
    GRID["for each (dx, dy) in axis_x.offsets x axis_y.offsets"] --> PERT["xp = x0; xp[idx_x]+=dx; xp[idx_y]+=dy"]
    PERT --> SIM["sim = simulate(xp, [t_early, t_final])"]
    SIM -->|non-finite or solver error| FAIL["failed[iy,ix] = True<br/>distances = NaN"]
    SIM -->|ok| DIFF["diff = (sim.x - baseline.x)[non_rotating_mask]"]
    DIFF --> DE["d_early = norm(diff at t_early)"]
    DIFF --> DL["d_late = norm(diff at t_final)"]
    DE --> CLASS["in_roa[iy,ix] = (d_late <= d_early)"]
    DL --> CLASS
    CLASS --> RESULT["RoaGridResult<br/>in_roa, early_distance, late_distance, failed"]
    FAIL --> RESULT
    BASE --> DIFF
```

## Two real design decisions, driven by wrong results, not guessed upfront

- **Trend, not fixed-time convergence.** An early version compared the
  distance at a single fixed end time against a threshold. It
  misclassified clearly-stable, only-moderately-perturbed points as
  "unstable" simply because WSCC-9's ~0.5–1.4s-period electromechanical
  modes hadn't finished a full oscillation yet by the time checked (a mild
  ±0.02 pu speed perturbation "failed" a fixed-endpoint check at t=0.5s).
  Comparing distance at two checkpoints — early, after the fast
  sub-millisecond transients settle; late, near the end of the window —
  and checking whether it *shrank* is a much better proxy for "is this
  converging" from an affordable short trajectory.
- **A solver failure is tracked, not misclassified as instability.** A
  Newton solve failing (seen at the most extreme grid corner) is a
  distinct outcome from "diverged" and must not be silently folded into
  `in_roa=False` — that would conflate "we don't know" with "no". Failed
  points get `NaN` distances and are excluded from any recovery-rate
  statistic, not counted as evidence either way.

## The "theta problem"

A rotating machine's absolute angle state (`theta`/`theta_pll`) grows
without bound even at a stable operating point (`dtheta/dt = wb`, ~377
rad/s) — it's an absolute angle, not a deviation. Any distance metric here
excludes those states (`_non_rotating_mask`) and relies on the *other*
physically coupled states (fluxes, currents, mechanical power, AVR states,
node voltages) to reveal a loss of synchronism indirectly — a real
simplification, not a rigorous treatment of angle differences between
machines, and documented as such rather than silently assumed.

Validated on a small (3×3), affordable (~100s) grid on WSCC-9, perturbing a
non-slack machine's rotor angle and speed deviation: the center
(zero-perturbation) point trivially matches its own baseline, and 7 of 8
surrounding points show real, substantial distance shrinkage consistent
with the system's known 12–49% modal damping.

## Physics: a trajectory-sampling ROA estimate

For each grid perturbation $(\delta_x, \delta_y)$ along two chosen state
axes, this integrates the nonlinear model
({doc}`timedomain <timedomain>`) from the perturbed initial condition
$x_p(0) = x_0 + \delta_x \hat e_x + \delta_y \hat e_y$ and compares it —
at two times, not one — against the unperturbed baseline trajectory
$x_b(t)$ (itself started from the exact operating point $x_0$, so
$x_b(t)$ isn't constant; it's the network's own natural response, which
a real equilibrium should stay bounded near). Excluding the unbounded
rotating-angle states via a mask $\mu$ (see the "theta problem" below),
the distance metric is

$$
d(t) = \big\lVert \big(x_p(t) - x_b(t)\big) \odot \mu \big\rVert_2
$$

evaluated at $t_{early}$ (past the sub-millisecond fast transients) and
$t_{final}$. The classification is a **trend**, not a threshold:

$$
\text{in\_ROA}(\delta_x,\delta_y) = \big[\,d(t_{final}) \le d(t_{early})\,\big]
$$

— the perturbed trajectory is judged to be *returning toward* the
equilibrium's own path if its distance from it shrank between the two
checkpoints, regardless of whether it has fully settled by
$t_{final}$ (see the module's own account of why a fixed-time threshold
misclassified genuinely stable, still-oscillating points as unstable).

**The angle mask**: $\mu_k = 0$ for any state $k$ whose name contains
`theta` ($\theta$, $\theta_{pll}$, $\theta_{up}$ — every rotating angle
in every component above), $\mu_k=1$ otherwise. These states grow
without bound even at a perfectly stable operating point
($\dot\theta=\omega_b\,\omega \approx \omega_b$, not a deviation), so
including them in $d(t)$ would make *every* trajectory read as
"diverging" regardless of true stability — a real modeling
simplification (a rigorous rotor-angle-difference treatment between
machines is not attempted here), stated as such rather than silently
assumed.

## Reference

```{eval-rst}
.. automodule:: g2elin_core.stability.roa
```
