# `timedomain` — EMT (nonlinear) time-domain simulation

Feature 2.2 from the migration plan. Reuses exactly the nonlinear
`f`/`g`/`h` callables built in {doc}`components <components>` and the same
interconnection topology used for linear modal analysis
({doc}`interconnect <interconnect>`) — the wiring between components
(KVL/KCL in the dq frame) is exact regardless of whether the components
themselves are linearized, so it doesn't need to be re-derived. Only the
per-component math changes, from "substitute numbers into a Jacobian" to
"call a nonlinear function and Newton-solve the coupled algebraic system."

**On the name.** This integrates the network's actual nonlinear
differential-algebraic equations forward in time — not a linearized
small-signal approximation — which is the property that makes this "EMT
simulation" for this project's purposes. It's a dq-frame, averaged-
converter, single-rotating-reference-frame formulation, so it doesn't
capture switching-level converter physics, abc-frame/unbalanced faults, or
distributed/traveling-wave line effects.

## Why a single *global* Newton solve, not one per component

Each component's own algebraic variables `z_i` are governed by its own
`g_i(x_i, z_i, u_i) = 0`. But a component's inputs `u_i` are wired to
*other* components' outputs `y_j = h_j(x_j, z_j, u_j)` — and since `h_j`
can itself depend on `u_j`, the full input vector `u` and every
component's `z` have to be solved **together**:

```
0 = g(x, z, u)                     — every component's own algebraic eqs
u = F @ u_exo + G @ h(x, z, u)     — the interconnection, exact and linear
```

at every state `x` the ODE integrator asks for. There's no way to solve
one component's `z_i` in isolation and plug it into the next.

```{mermaid}
flowchart TD
    XU["x (current ODE state), u_exo"] --> RES["_residual(z, u)<br/>= [ g(x,z,u) ; u - (F@u_exo + G@h(x,z,u)) ]"]
    RES --> NEWTON["scipy.optimize.root(method='hybr')<br/>with analytic jac=_residual_jacobian"]
    JAC["_residual_jacobian(z, u)<br/>assembled per-block from Gz/Gu/Hz/Hu<br/>+ interconnection's G matrix"] --> NEWTON
    NEWTON --> SOL["(z, u) satisfying the coupled system"]
    SOL --> RHS["rhs(): xdot = concat(f_i(x_i, z_i, u_i))"]
    RHS --> IVP["scipy.integrate.solve_ivp<br/>(Radau — stiff, L-stable)"]
    IVP -->|next step, warm-started z/u guess| XU
    IVP --> RESULT["EmtSimulationResult<br/>t, x, state_names"]
```

**The analytic Jacobian is the load-bearing optimization here.** An
earlier version used SciPy's numerically-estimated Jacobian:
`method="hybr"` converged fine exactly at the operating point but reliably
failed ("not making good progress") for the nearby points an ODE
integrator's own adaptive stepping evaluates; `method="lm"` was robust but
slow (~0.05s/solve, re-estimating its own Jacobian every iteration). The
analytic `_residual_jacobian` — assembled the same way the linear A/B/C/D
elimination combines per-component Jacobians with the interconnection's
`G` matrix, just kept as an explicit matrix instead of solved in closed
form — let `"hybr"` come back both faster (~0.006s/solve, ~9x) and at
least as robust.

## `NonlinearBlockComp` and `NonlinearNetworkModel`

```{mermaid}
classDiagram
    class NonlinearBlockComp {
        n_states n_z n_us n_ug n_out_s n_out_g : int
        f(x,z,u) g(x,z,u) h(x,z,u)
        Gz Gx Gu Hz Hx Hu : Callable
        x0 z0 u0 : ndarray
    }
    class NonlinearNetworkModel {
        blocks : list~Block~
        topology : Topology
        initial_state() ndarray
        solve_algebraic(x, u_exo) tuple
        rhs(x, u_exo) tuple
        state_names : list~str~
    }
    NonlinearNetworkModel "1" *-- "many" NonlinearBlockComp : via Block.comp
```

`NonlinearBlockComp` satisfies the same `PortSpec` protocol
`LinearComponent` does (see {doc}`interconnect <interconnect>`), so it
slots into the exact same `Block`/`Wiring`/`compute_topology()` machinery —
`NonlinearNetworkModel` is the nonlinear sibling of `AssembledSystem`,
built by `build_nonlinear_network()` the same way `pipeline.linearize_network()`
builds the linear one, down to sharing the same "unsupported DER unit type"
and slack-must-be-SM-or-IB restrictions.

`simulate()` integrates the DAE with `scipy.integrate.solve_ivp`, defaulting
to `Radau` (implicit, L-stable) because the system's stiffest modes reach
~1e6–1e7 rad/s. Each RHS evaluation warm-starts its Newton solve from the
previous step's `(z, u)` solution (`cache` in `simulate()`), which is both
a speed optimization and — since the initial per-component operating-point
guess and the true coupled-network equilibrium differ by a small but real
amount (~0.08 for WSCC-9/CIGRE, from the node-susceptance MATLAB quirk) —
closer to the actual solution than the naive per-component guess would be.

**`solve_algebraic` now retries with `"lm"` on an `"hybr"` failure**,
found necessary (not preemptive) when wiring the SMIB presets in: SMIB's
tighter 2-DER coupling produces an initial-guess offset of ~0.23 pu in the
reactive-current direction — large enough that `"hybr"` (even with the
analytic Jacobian) can't reach the true equilibrium from that starting
point, though `"lm"` finds it to machine precision from the identical
guess, confirming a real equilibrium exists nearby and this is a basin-of-
convergence issue, not a modeling bug.

**`u_exo_fn` is now actually used.** The module docstring for `simulate()`
always described it — "lets a caller drive a disturbance (e.g. a `P_ref`
step)" — but nothing called it with a non-default function until the EMT
tab gained the ability to perturb an *input* (a permanent step in one
exogenous reference, held from t=0) as an alternative to perturbing a
*state* (an initial-condition offset). `g2elin_api.main.run_emt` builds a
closure returning the stepped `default_u_exo()` vector for every `t`, and
threads the *same* closure into `recover_inputs_and_outputs()` too — using
the unperturbed default there instead would silently contradict the
trajectory that was actually integrated.

## Recovering inputs and outputs after the fact

`solve_ivp` only returns state trajectories — the algebraic variables
`z`/full inputs `u` it solves for internally live at the integrator's own
adaptive step points, not necessarily the (possibly interpolated) sample
times it returns. `NonlinearNetworkModel.recover_inputs_and_outputs(t, x)`
re-solves the coupled algebraic system once per already-integrated sample
(warm-started from the previous one, same as during integration) to
recover a `(z, u)` consistent with that exact `(t[i], x[:, i])`, then
evaluates each block's own named outputs `h(x, z, u)[:n_out_s]` there —
the same `output_names` `AssembledSystem` exposes for the linear model,
since both come from the same `ComponentDAE.output_names`. Costs about as
much again as the integration itself, so the EMT API endpoint only calls
it when the caller actually asked for input/output trajectories (see
{doc}`../api_and_web`).

## Physics: the coupled nonlinear DAE

Every component contributes its own **nonlinear** differential-algebraic
system — literally the same `diffeqVec`/`algeqVec`/`outputeqVec`
{doc}`components <components>` documents, evaluated as functions instead
of linearized:

$$
\dot x_i = f_i(x_i, z_i, u_i)
\qquad
0 = g_i(x_i, z_i, u_i)
\qquad
y_i = h_i(x_i, z_i, u_i)
$$

Stacked over every block and coupled through the exact same
interconnection selector matrices {doc}`interconnect <interconnect>`
builds ($u = F\,u_{exo} + G\,y$, $y=h(x,z,u)$ now nonlinear instead of
$Cx+Du$), the full network's algebraic state at a fixed $x$ solves the
**residual**

$$
R(z,u) = \begin{bmatrix} g(x,z,u) \\[2pt] u - \big(F\,u_{exo} + G\,h(x,z,u)\big) \end{bmatrix} = 0
$$

via Newton's method, using the **analytic** residual Jacobian (not a
numerically-estimated one — see below for why that distinction matters)

$$
\frac{\partial R}{\partial(z,u)} =
\begin{bmatrix} G_z & G_u \\[2pt] -G\,H_z & I - G\,H_u \end{bmatrix}
$$

assembled per-block from each component's own $G_z, G_u, H_z, H_u$
Jacobians (evaluated at the *current* Newton iterate, not once at a fixed
operating point — the crucial difference from {doc}`components
<components>`'s linearization, which only ever needs one such evaluation).
Once $(z,u)$ solves $R=0$ at the current $x$, the state derivative
$\dot x = f(x,z,u)$ is the right-hand side `scipy.integrate.solve_ivp`
(implicit `Radau`) integrates forward — one full Newton solve per RHS
evaluation the ODE integrator requests, warm-started from the previous
step's $(z,u)$.

## Choosing a solver

The integration settings live on the request (`EmtRequest`), not the
network: they are about how a trajectory is computed, not about what is
being modelled. `GET /api/solvers` serves the catalogue the web UI's picker
is built from.

### The Jacobian, first

Radau, BDF and LSODA are implicit, so each needs `d(xdot)/dx`. Given none,
SciPy finite-differences it — **one full coupled Newton solve per state**,
and SciPy does not count those in `sol.nfev`, so the cost is invisible in
the diagnostics. Measured on WSCC-9 at full order: 257 real RHS evaluations
against a reported `nfev` of 81, so 68% of the work was building a matrix.

{meth}`~g2elin_core.timedomain.emt.NonlinearNetworkModel.ode_jacobian`
assembles it instead. Differentiating the constraints along the trajectory
gives a linear system whose matrix is the *same* Newton Jacobian the
algebraic solve already factorises, so it costs one sparse factorisation
plus `n_x` back-substitutions. It is passed automatically to every implicit
method; there is no setting for it.

What that changed, measured:

| Case | Before | After |
|---|---|---|
| WSCC-9, full order, Radau | 650 s (320,714 RHS calls) | 0.84 s (141 calls) |
| IEEE-39, full order | 456 s | 267 s |
| IEEE-118, full order | did not finish in 30 min | 189 s |

The first row is the pathology in miniature: without a usable Jacobian,
Radau's inner Newton fails, it rejects the step, shrinks it, rebuilds the
Jacobian, and spirals.

### Method

| Method | When |
|---|---|
| `Radau` | the default. L-stable, copes with the full model's 1e7 rad/s modes. |
| `BDF` | usually the fastest. Measured 1.4–5× quicker than Radau on the same run. Less robust across a sharp event. |
| `LSODA` | detects stiffness and switches. A reasonable "don't know" choice. |
| `RK45`, `DOP853` | explicit, no Jacobian at all — but they must resolve every fast mode. |

```{warning}
The explicit methods are offered for completeness and were slower in every
case measured. On the full EMT model they diverge outright; even on a
quasi-stationary model whose fastest mode is ~6 krad/s, RK45 took 952 steps
and 61 s where BDF took 39 steps and 1.5 s. Model-order reduction removes
several orders of stiffness but not enough to make an explicit method
competitive.
```

### Tolerance

`rtol` is the largest single lever after the method. Loosening it from 1e-4
to 1e-3 roughly halved one run's time — and moved the peak of the signal by
7%, which is the trade in one sentence. 1e-3 for screening, 1e-6 for a
number that goes in a paper.

`max_step` is best left unset. It bounds the solver's step, which only
forces steps it did not need; the *output* spacing is a separate setting.

### Fixed stepping

{func}`~g2elin_core.timedomain.emt.simulate_fixed_step` integrates on a step
you choose, with the trapezoidal rule — how an EMT program works. It gives
up error control and gets back a cost you can predict before starting, output
exactly on the grid you asked for with no interpolation, and a result that
does not move when a tolerance is nudged.

It is not usually *faster*: on a smooth trajectory a variable-step solver
takes big steps where nothing happens and wins easily (measured: 5.2 s at a
2 ms fixed step against 2.8 s for BDF over the same run). Reach for it when
you want reproducibility, a direct comparison against PSCAD or EMTP, or a
run whose duration you can promise.

Trapezoidal is A-stable, so too large a step cannot make it blow up — but it
can make it *ring*, the numerical oscillation EMT programs damp
deliberately. A trace oscillating at exactly two samples per cycle is that,
and the fix is a smaller step, not a smaller tolerance.

## Reference

```{eval-rst}
.. automodule:: g2elin_core.timedomain.emt
```
