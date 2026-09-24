# `components` — symbolic nonlinear-then-linearized models

Ports the pattern used throughout `Functions/sym*.m`: build a full
**nonlinear** DAE symbolically with [sympy](https://www.sympy.org/)
(differential equations, algebraic constraints, outputs), then either
linearize it around an operating point (what every MATLAB `sym*.m` file
did, and all it did) or keep it nonlinear and lambdify it to a fast numeric
callable (new — this is what makes {doc}`EMT simulation <timedomain>`
possible without re-deriving any physics).

## The shared machinery (`components.base`)

```{mermaid}
classDiagram
    class ComponentDAE {
        Fx Fz Fu Gx Gz Gu Hx Hz Hu : sp.Matrix
        state_syms alg_syms input_syms : list~Symbol~
        diffeq_exprs algeq_exprs output_exprs : list~Expr~
        param_syms() list~Symbol~
        point_from_subs(subs) tuple
        linearize(subs) LinearComponent
        nonlinear_funcs() NonlinearFuncs
        nonlinear_jacobians() NonlinearJacobians
    }
    class LinearComponent {
        A B C D : ndarray
        n_us n_ug n_out_s n_out_g : int
        state_names input_names output_names
    }
    class NonlinearFuncs {
        f(x,z,u,p) ndarray
        g(x,z,u,p) ndarray
        h(x,z,u,p) ndarray
    }
    class NonlinearJacobians {
        Gz(x,z,u,p) Gx(x,z,u,p) Gu(x,z,u,p)
        Hz(x,z,u,p) Hx(x,z,u,p) Hu(x,z,u,p)
    }
    ComponentDAE --> LinearComponent : linearize(subs)
    ComponentDAE --> NonlinearFuncs : nonlinear_funcs()
    ComponentDAE --> NonlinearJacobians : nonlinear_jacobians()
```

`build_dae()` is the one function every `sym*`-ported module calls: given
symbolic state/algebraic/input vectors and the raw nonlinear
`diffeq`/`algeq`/`output` expressions, it takes every Jacobian
(`Fx = d(diffeq)/d(state)`, etc. — mirroring the `jacobian(diffeqVec,
stateVec)` calls repeated in every MATLAB `sym*.m` file) and returns one
`ComponentDAE`. Three consumers read from that same object:

- **`linearize(subs)`** substitutes numeric parameter + operating-point
  values and eliminates the algebraic variables (`A = Fx - Fz@inv(Gz)@Gx`,
  etc.) — this is *all* the MATLAB source ever did.
- **`nonlinear_funcs()`** lambdifies the raw `diffeq`/`algeq`/`output`
  expressions (unsubstituted) into `(x, z, u, p) -> ndarray` callables —
  the same equations, kept instead of discarded after differentiation.
- **`nonlinear_jacobians()`** lambdifies `Gz`/`Gx`/`Gu`/`Hz`/`Hx`/`Hu`
  themselves (not substituted at one fixed point) into `(x,z,u,p) ->
  ndarray` callables — what a Newton solver needs to evaluate the
  algebraic Jacobian at whatever point it's currently iterating on, not
  just once at a fixed operating point. This is the analytic Jacobian
  {doc}`timedomain`'s coupled Newton solve uses.

One deliberate deviation from the MATLAB source is worth calling out: the
original toolbox substitutes numeric values *before* forming
`Ai = Fx - Fz*inv(Gz)*Gx` (symbolic matrix inversion). Substitution and
matrix algebra commute, so this port keeps `Fx..Hu` symbolic (cheap — just
differentiation) and defers substitution to `linearize()`, where the `Gz`
inverse runs numerically in numpy instead of symbolically in sympy —
avoiding symbolic inversion of a 15×15+ matrix. Same result, much faster.

## Per-component modules

Each ported component follows the same shape: a module-level `*_dae()`
function (cached with `@lru_cache`, since building the symbolic DAE is the
expensive part and every instance of a given component type shares the same
equations), an `*OperatingPoint` class that turns power-flow results into
substitution values, and `linearize_*()` / `*_nonlinear_point()` /
`*_nonlinear_funcs()` / `*_nonlinear_jacobians()` wrapper functions.

```{mermaid}
flowchart LR
    subgraph Sources["Functions/sym*.m"]
        symSG["symSG.m + SG_subs.m"]
        symGFM["symGFM_types.m + GFM_subs.m"]
        symGFL["symGFL.m + GFL_subs.m"]
        symIB["symIB.m + IB_subs.m"]
        symLine["symLine.m"]
        symNode["symNode.m"]
        symLoad["symLoad.m"]
    end
    subgraph Python["components/*.py"]
        sm["sm.py<br/>synchronous machine"]
        gfm["gfm.py<br/>grid-forming"]
        gfl["gfl.py<br/>grid-following"]
        ib["ib.py<br/>infinite bus"]
        frame["frame.py<br/>reference frame"]
        line["line.py"]
        node["node.py"]
        load["load.py"]
    end
    symSG --> sm --> DAE["build_dae() -> ComponentDAE"]
    symGFM --> gfm --> DAE
    symGFL --> gfl --> DAE
    symIB --> ib --> DAE
    symLine --> line --> DAE
    symNode --> node --> DAE
    symLoad --> load --> DAE
```

**Coverage note**: `gfm.py` implements only the `'Droop'` branch of
`symGFM_types.m` (the CIGRE preset's GFM controller); `'Droop+filter'`,
`'Droop+PLL'`, `'dVOC'`, `'VSM'` and `'Matching'` aren't ported.
`ib.py` (infinite bus) was built and tested in isolation for a while
before being wired into {doc}`interconnect <interconnect>` — no MATLAB
script in the source repo defines a *complete* infinite-bus DER case
(`preset_networks.m` names `SM_SMIB`/`GFM_SMIB`/`GFL_SMIB` but their
`switch` bodies are empty), so there was no real preset to validate it
against until the SMIB presets ({doc}`network <network>`) reconstructed
one from the real topology data that *does* exist (`SMIB_raw.m`).

A major bug was found and fixed while building `sm.py`'s
`SmOperatingPoint`: a double-rotation error in the terminal-voltage
phasor (`v0 * exp(j*(angle_terminal_rad - theta0))` instead of
`v0 * exp(-j*theta0)`) silently corrupted every non-slack synchronous
machine's operating point, invisible on the slack unit only because
pandapower pins its angle to exactly 0. See the `SmOperatingPoint`
docstring below and its internal consistency assertion (comparing
power-flow-derived `id0`/`iq0` against independently-derived values from
the flux-current relations) — that assertion is what caught it.

## Physics: the full nonlinear equations

Everything below is the literal content of each `*_dae()` function — the
same `diffeqVec`/`algeqVec`/`outputeqVec` sympy builds, transcribed to
LaTeX rather than paraphrased, so this can be read as the model
specification and cross-checked line-by-line against the source. Two
notational conventions hold throughout, matching {doc}`pu_base
<operating_point>`'s amplitude-invariant convention and the MATLAB
source's own dq-frame formulation:

- Everything is per-unit on the network's own `sn_mva`/bus `vn_kv` base
  ({doc}`pu_base <operating_point>`), and every state/algebraic/input
  symbol below is exactly the sympy symbol name used in the code, just
  set in math instead of `snake_case` (e.g. `Dwr` → $\Delta\omega_r$,
  `phi_d` → $\psi_d$, `theta_g` → $\theta_g$).
- Two reference frames are in play. A **global** dq frame carries every
  node, line and load (their states already have a `_g` suffix), turning at
  $\omega_g$ with angle $\theta_g$. Each DER (SM/GFM/GFL/IB) instead has its
  **own** rotating frame, at its own angle $\theta$ (SM, GFM), $\theta_{pll}$
  (GFL, PLL-locked) or $\theta_{up}$ (IB) — the algebraic equations below
  rotate a DER's own terminal quantities into and out of the global frame by
  the angle difference $\theta_g - \theta$.
- The global frame is a **block of its own** ({doc}`frame.py <components>`):
  one state, the frame angle, following the speed of the unit its island is
  referenced to, and nothing else. It replaces the toolbox's convention of
  reading $\theta_g$/$\omega_g$ off the slack unit, which made that unit
  impossible to disconnect: with the frame standing apart, any unit can be
  tripped and a network can split into islands, each with its own frame and
  its own frequency. `Network.frame_follows_slack` restores the old wiring,
  where the slack unit itself is the frame and has to be a synchronous
  machine or an infinite bus. The two give the same modes (checked to 1e-5
  relative in `tests/test_breakers.py`) and, being only a change of
  coordinates, the same nonlinear trajectories.

A state's time derivative is written $\dot x$ throughout; every equation
below carries an explicit $\omega_b$ (`wb`, the base angular frequency,
$2\pi f_n$) factor exactly where the code has it — per-unit flux-linkage
and current states pick up this factor because $1$ per-unit flux
corresponds to $1/\omega_b$ seconds of a unit voltage integrated, the
standard per-unit machine-equation convention.

### Synchronous machine (`sm.py`)

A round-rotor synchronous machine with one field winding and 2 damper
windings (1 $d$-axis, 2 $q$-axis), its own step-up transformer, and a set of
regulators — 19 states in its default configuration, ported from
`Functions/symSG.m`.

**States** $x = [i_{gd}, i_{gq}, \psi_d, \psi_q, \psi_{fd}, \psi_{1d},
\psi_{1q}, \psi_{2q}, \Delta\omega_r, \theta, \ldots]$ — grid-side
transformer current, stator/field/damper flux linkages, speed deviation and
rotor angle, followed by whatever the fitted regulators carry (below).

#### Selectable regulators

The exciter, the stabiliser and the governor are chosen per unit through
`DerUnit.exciter`, `DerUnit.pss` and `DerUnit.governor`. They decide the
machine's states *and* its parameter set: each model brings its own
parameters under its own names, so what is typed is what that model's own
diagram shows, and an override under a name the fitted model does not have
is rejected rather than ignored.

| Slot | Model | States | Parameters |
|---|---|---|---|
| Exciter | `g2elin` (default) | $e_1, e_2, e_{fd}, e_3$ | $T_r, K_a, T_a, K_e, T_e, K_{fd}, T_{fd}$ |
| | `kundur` | $e_1, e_{tgr}$ | $T_R, K_A, T_A, T_B$ |
| PSS | `g2elin` (default) | $\Delta\omega_1, v_1, v_2, v_{pss}$ | $K_{PSS}, T_{LP}, T_{HP}, T_{1n}, T_{1d}, T_{2n}, T_{2d}$ |
| | `kundur` | $v_1, v_2, v_{pss}$ | $K_{STAB}, T_W, T_1, T_2, T_3, T_4$ |
| | `none` | — | — |
| Governor | `g2elin` (default) | $P_m$ | $m_p, T_G$ |
| | `none` | — | — |

The `kundur` pair is Kundur Fig. E12.9: a thyristor exciter — a terminal
voltage transducer, a gain and transient gain reduction, with no exciter lag
— and a stabiliser that is a washout and two lead-lag stages. Because a
thyristor bridge has no time constant of its own, $E_{fd}$ follows the
regulator instantaneously and is **algebraic** there rather than a state,
which is why that exciter is two states shorter rather than one.

`none` is equipment that is not fitted, not a gain turned down: the states,
the parameters and the reduction group all go. A machine with no governor
holds its mechanical power at the operating point, which is the
constant-torque assumption most textbook small-signal examples are worked
under — and note that a fleet with no governor anywhere has an unregulated
common frequency, which the modal tools treat as marginal by construction
rather than as an instability.

**Algebraic variables** $z = [\omega_r, v_{ed}, v_{eq}, i_d, i_{1d},
i_{fd}, i_q, i_{1q}, i_{2q}, v_{gd}, v_{gq}, C_m, C_e, \Delta P, E_t]$.

**Inputs**: own $u_s = [V_{ref}, P_{ref}, \omega_{ref}]$ — note that
$P_{ref}$ is the *governor's* setpoint, so with no governor fitted its
column of $B$ is exactly zero and it moves nothing; grid-coupling
$u_g = [\theta_g, v_{gd,g}, v_{gq,g}]$ (just $[v_{gd,g}, v_{gq,g}]$ for the
slack unit, which needs no $\theta_g$ input since $\theta_g \equiv \theta$
for it by definition).

**Electrical dynamics** (transformer current, stator/damper flux):

$$
\begin{aligned}
\dot i_{gd} &= \frac{\omega_b}{L_t}\left(v_{ed} - v_{gd} - R_t i_{gd} + \omega_r L_t i_{gq}\right) \\
\dot i_{gq} &= \frac{\omega_b}{L_t}\left(v_{eq} - v_{gq} - R_t i_{gq} - \omega_r L_t i_{gd}\right) \\
\dot\psi_d &= \omega_b\left(v_{ed} + R_a i_d + \omega_r \psi_q\right) \\
\dot\psi_q &= \omega_b\left(v_{eq} + R_a i_q - \omega_r \psi_d\right) \\
\dot\psi_{fd} &= \omega_b\frac{R_{fd}}{L_{ad}} e_{fd} - \omega_b R_{fd} i_{fd} \\
\dot\psi_{1d} &= -\omega_b R_{1d} i_{1d} \qquad
\dot\psi_{1q} = -\omega_b R_{1q} i_{1q} \qquad
\dot\psi_{2q} = -\omega_b R_{2q} i_{2q}
\end{aligned}
$$

**Mechanical (swing equation + first-order turbine/governor)**:

$$
\Delta\dot\omega_r = \frac{C_m - C_e - K_D\,\Delta\omega_r}{2H}
\qquad
\dot\theta = \omega_b\,\omega_r
\qquad
\dot P_m = \frac{P_{ref} - \Delta P - P_m}{T_G}
$$

**PSS** (washout $\to$ lead-lag cascade — note $\dot v_1$ feeds directly
into the $\dot v_2$ equation and $\dot v_2$ into $\dot v_{pss}$, a
transfer-function realization, not three independent first-order lags):

$$
\Delta\dot\omega_1 = \frac{\Delta\omega_r - \Delta\omega_1}{T_{LP}}
\qquad
\dot v_1 = \frac{T_{HP} K_{PSS}\,\Delta\dot\omega_1 - v_1}{T_{HP}}
$$

$$
\dot v_2 = \frac{T_{1n}\dot v_1 + v_1 - v_2}{T_{1d}}
\qquad
\dot v_{pss} = \frac{T_{2n}\dot v_2 + v_2 - v_{pss}}{T_{2d}}
$$

**AVR** (voltage transducer $\to$ PI-ish comparator $\to$ exciter $\to$
field-current-derivative feedback):

$$
\dot e_1 = \frac{E_t - e_1}{T_r}
\qquad
\dot e_2 = \frac{K_a(V_{ref} - e_1 - e_3 + v_{pss}) - e_2}{T_a}
$$

$$
\dot e_{fd} = \frac{K_e e_2 - e_{fd}}{T_e}
\qquad
\dot e_3 = \frac{K_{fd}\,\dot e_{fd} - e_3}{T_{fd}}
$$

**Algebraic constraints** — frequency reference, stator-side aux-load
Ohm's law, flux/current inversion, torque/power, and the terminal-voltage
magnitude:

$$
\omega_r = \Delta\omega_r + \omega_{ref}
\qquad
v_{ed} = (i_d - i_{gd})R_g
\qquad
v_{eq} = (i_q - i_{gq})R_g
$$

$$
\begin{bmatrix} i_d \\ i_{1d} \\ i_{fd} \end{bmatrix}
= \begin{bmatrix} -L_{ad}-L_l & L_{ad} & L_{ad} \\ -L_{ad} & L_{1d}+L_{ad} & L_{ad} \\ -L_{ad} & L_{ad} & L_{fd}+L_{ad} \end{bmatrix}^{-1}
\begin{bmatrix} \psi_d \\ \psi_{1d} \\ \psi_{fd} \end{bmatrix}
\qquad
\begin{bmatrix} i_q \\ i_{1q} \\ i_{2q} \end{bmatrix}
= \begin{bmatrix} -L_{aq}-L_l & L_{aq} & L_{aq} \\ -L_{aq} & L_{1q}+L_{aq} & L_{aq} \\ -L_{aq} & L_{aq} & L_{2q}+L_{aq} \end{bmatrix}^{-1}
\begin{bmatrix} \psi_q \\ \psi_{1q} \\ \psi_{2q} \end{bmatrix}
$$

$$
C_m = \frac{P_m}{\omega_r}
\qquad
C_e = \psi_d i_q - \psi_q i_d
\qquad
\Delta P = \frac{\omega_r - \omega_{ref}}{m_p}
\qquad
E_t = \sqrt{v_{ed}^2 + v_{eq}^2}
$$

Frame rotation into the machine's own dq frame (non-slack only — the
slack sets $v_{gd} = v_{gd,g}$, $v_{gq} = v_{gq,g}$ directly, since its
own frame *is* the global one):

$$
\begin{bmatrix} v_{gd} \\ v_{gq} \end{bmatrix} =
\begin{bmatrix} \cos(\theta_g-\theta) & -\sin(\theta_g-\theta) \\ \sin(\theta_g-\theta) & \cos(\theta_g-\theta) \end{bmatrix}
\begin{bmatrix} v_{gd,g} \\ v_{gq,g} \end{bmatrix}
$$

**Outputs**: own $[P_e, Q_e, \omega_r, \theta, \Delta\omega_r, E_t]$ with
$P_e = v_{ed}i_{gd}+v_{eq}i_{gq}$, $Q_e=-v_{ed}i_{gq}+v_{eq}i_{gd}$; grid
current, rotated back to the global frame for a non-slack unit (the
slack outputs $i_{gd}, i_{gq}$ directly, being already in that frame):

$$
\begin{bmatrix} i_{gd,g} \\ i_{gq,g} \end{bmatrix} =
\begin{bmatrix} \cos(\theta-\theta_g) & -\sin(\theta-\theta_g) \\ \sin(\theta-\theta_g) & \cos(\theta-\theta_g) \end{bmatrix}
\begin{bmatrix} i_{gd} \\ i_{gq} \end{bmatrix}
$$

Default electrical/mechanical/AVR/PSS parameter values are tabulated in
{doc}`operating_point <operating_point>`.

### Grid-forming converter (`gfm.py`)

A two-level VSC behind an LCL filter ($R_f/L_f$ series, $C_f$ shunt) and its
own step-up transformer ($R_t/L_t$), with an outer power-control law, a
voltage-control loop, a current-control loop, and a first-order DC-link
model — 15 states running droop, ported from `Functions/symGFM_types.m`.

**States** $x = [i_{sd}, i_{sq}, i_{gd}, i_{gq}, v_{dc}, i_{dc}, \ldots,
M_{VLd}, M_{VLq}, M_{CLd}, M_{CLq}]$ — converter and grid-side filter
currents, filter capacitor voltage, DC-link voltage and current, then the
outer law's own states (below), then the voltage-loop/current-loop PI
integrator states.

#### Selectable power-control laws

All five of `symGFM_types.m`'s laws are chosen per unit through
`DerUnit.controller`. They differ in exactly three things — which states the
outer loop carries, how it forms the frequency deviation $\Delta\omega$ that
the angle integrates, and how it forms the voltage reference the cascade
then tracks. Everything downstream is identical for all five.

| Law | Outer states | $\Delta\omega$ | $v_{ed}^{ref}$ | Parameters |
|---|---|---|---|---|
| `droop` (default) | $p_m, \theta, q_m$ | $m_p(P_{ref}-p_m)$ | $V_{ref} + (Q_{ref}-q_m)n_q$ | $m_p, n_q, w_f$ |
| `droop_filtered` | $p_m, \Delta\omega, \theta, q_m$ | a state | same | $m_p, n_q, w_f, w_c$ |
| `dvoc` | $p_m, \theta, q_m, v_{ed}^{ref}$ | $\eta(P_{ref}/V_{ref}^2 - p_m/(v_{ed}^{ref})^2)$ | a state | $\eta, \alpha, w_f$ |
| `vsm` | $\Delta\omega, \theta, \Phi$ | a state (swing, $J$/$D_p$) | $\omega\Phi$ | $J, D_p, K, D_q$ |
| `matching` | $v_{dc,m}, \theta$ | $K_\theta(v_{dc,m} - V_{dc,ref})$ | $V_{ref}$ | $K_\theta, w_f$ |

**A virtual synchronous machine and matching control carry no power filters
at all.** VSM measures $p$ and $q$ directly — its emulated inertia is what
smooths the response — and matching takes its frequency from the DC link,
which is what a machine's speed does physically. So $p_m$/$q_m$ are absent
from those two rather than merely retuned.

The gains are `script_generic.m`'s, and every one is written in terms of the
droop tuning ($\eta = m_p$, $J = 1/(m_p w_f)$, $K_\theta = m_p K_{pdc}$, and
so on). That is deliberate: the laws are meant to be *comparable*, tuned to
the same equivalent inertia and reactive gain, so a study that swaps one for
another sees what the law changes rather than what a different tuning
changes.

**Algebraic** $z = [m_d, m_q, \omega]$ (modulation index, own frequency).
**Inputs**: $u_s = [P_{ref}, Q_{ref}, V_{ref}, \omega_{ref}, V_{dc,ref}]$,
$u_g = [\theta_g, v_{gd,g}, v_{gq,g}]$ (a GFM is never the slack).

$$
\begin{bmatrix} v_{gd} \\ v_{gq} \end{bmatrix} =
\begin{bmatrix} \cos(\theta_g-\theta) & -\sin(\theta_g-\theta) \\ \sin(\theta_g-\theta) & \cos(\theta_g-\theta) \end{bmatrix}
\begin{bmatrix} v_{gd,g} \\ v_{gq,g} \end{bmatrix}
$$

**LCL filter + DC link**:

$$
\begin{aligned}
\dot i_{sd} &= \frac{\omega_b}{L_f}\left(m_d v_{dc} - v_{ed} - R_f i_{sd} + \omega L_f i_{sq}\right) \\
\dot i_{sq} &= \frac{\omega_b}{L_f}\left(m_q v_{dc} - v_{eq} - R_f i_{sq} - \omega L_f i_{sd}\right) \\
\dot i_{gd} &= \frac{\omega_b}{L_t}\left(v_{ed} - v_{gd} - R_t i_{gd} + \omega L_t i_{gq}\right) \\
\dot i_{gq} &= \frac{\omega_b}{L_t}\left(v_{eq} - v_{gq} - R_t i_{gq} - \omega L_t i_{gd}\right) \\
\dot v_{ed} &= \frac{\omega_b}{C_f}\left(i_{sd} - i_{gd} + \omega C_f v_{eq}\right) \\
\dot v_{eq} &= \frac{\omega_b}{C_f}\left(i_{sq} - i_{gq} - \omega C_f v_{ed}\right) \\
\dot v_{dc} &= \frac{\omega_b}{C_{dc}}\left(i_{dc} - G_{dc} v_{dc} - m_d i_{sd} - m_q i_{sq}\right) \\
\dot i_{dc} &= \frac{1}{T_{dc}}\left(\frac{P_{ref}}{V_{dc,ref}} + K_{pdc}(V_{dc,ref}-v_{dc}) - i_{dc}\right)
\end{aligned}
$$

**Droop + voltage/current loops** — $P$/$\omega$ and $Q$/$V$ droop, a
low-pass power-measurement filter ($\omega_f$), then a voltage loop
(PI + feed-forward, gains $K_{p/i,VL}$) whose output is the current
reference, then a current loop (PI + feed-forward, gains $K_{p/i,CL}$):

$$
q = -v_{ed}i_{gq}+v_{eq}i_{gd}
\qquad
p = v_{ed}i_{gd}+v_{eq}i_{gq}
\qquad
\Delta\omega = m_p(P_{ref}-p_m)
$$

$$
\dot p_m = \omega_f(p-p_m)
\qquad
\dot q_m = \omega_f(q-q_m)
\qquad
\dot\theta = \omega_b\,\omega
$$

$$
v_{ed,ref} = V_{ref}+(Q_{ref}-q_m)n_q
\qquad
v_{eq,ref}=0
$$

$$
\dot M_{VLd} = K_{iVL}(v_{ed,ref}-v_{ed})
\qquad
\dot M_{VLq} = K_{iVL}(v_{eq,ref}-v_{eq})
$$

$$
i_{sd,ref} = K_{pVL}(v_{ed,ref}-v_{ed}) + M_{VLd} + K_{ffi}i_{gd} - \omega_{ff}C_f v_{eq}
$$
$$
i_{sq,ref} = K_{pVL}(v_{eq,ref}-v_{eq}) + M_{VLq} + K_{ffi}i_{gq} + \omega_{ff}C_f v_{ed}
$$

$$
\dot M_{CLd}=K_{iCL}(i_{sd,ref}-i_{sd})
\qquad
\dot M_{CLq}=K_{iCL}(i_{sq,ref}-i_{sq})
$$

**Algebraic** (modulation index, own frequency):

$$
m_d = \frac{1}{v_{dc}}\Big(K_{pCL}(i_{sd,ref}-i_{sd}) + K_{ffv}v_{ed} - \omega_{ff}L_f i_{sq} + M_{CLd}\Big)
$$
$$
m_q = \frac{1}{v_{dc}}\Big(K_{pCL}(i_{sq,ref}-i_{sq}) + K_{ffv}v_{eq} + \omega_{ff}L_f i_{sd} + M_{CLq}\Big)
$$
$$
\omega = \Delta\omega + \omega_{ref}
$$

**Outputs**: $p$, $q$, $\omega$, $V_t=\sqrt{v_{ed}^2+v_{eq}^2}$, and the
grid current rotated back to the global frame the same way as the SM's
non-slack case above (using $\theta$ in place of the SM's own angle).

Default filter/loop-tuning values (pole-placement formulas for
$K_{p/i,VL}$, $K_{p/i,CL}$, $K_{pdc}$) are in {doc}`operating_point
<operating_point>`.

### Grid-following converter (`gfl.py`)

Structurally the same LCL-filter VSC as the GFM above, but current-
(not voltage-) controlled, phase-locked to the grid via an SRF-PLL
instead of self-clocking — 14 states, ported from `Functions/symGFL.m`.
Never the slack.

**States** $x = [i_{sd}, i_{sq}, i_{gd}, i_{gq}, v_{ed}, v_{eq}, v_{dc},
i_{dc}, M_d, M_q, M_{CLd}, M_{CLq}, M_{pll}, \theta_{pll}]$. **Algebraic**
$z=[m_d, m_q, \omega_{pll}]$. **Inputs**: $u_s=[V_{dc,ref}, Q_{ref},
I_{dc,ref}]$, $u_g=[\theta_g, v_{gd,g}, v_{gq,g}]$.

$$
\begin{bmatrix} v_{gd} \\ v_{gq} \end{bmatrix} =
\begin{bmatrix} \cos(\theta_g-\theta_{pll}) & -\sin(\theta_g-\theta_{pll}) \\ \sin(\theta_g-\theta_{pll}) & \cos(\theta_g-\theta_{pll}) \end{bmatrix}
\begin{bmatrix} v_{gd,g} \\ v_{gq,g} \end{bmatrix}
$$

**LCL filter + DC link** (identical form to the GFM's, $\omega$ replaced
by the PLL's own tracked frequency $\omega_{pll}$):

$$
\begin{aligned}
\dot i_{sd} &= \frac{\omega_b}{L_f}\left(m_d v_{dc} - v_{ed} - R_f i_{sd} + \omega_{pll} L_f i_{sq}\right) \\
\dot i_{sq} &= \frac{\omega_b}{L_f}\left(m_q v_{dc} - v_{eq} - R_f i_{sq} - \omega_{pll} L_f i_{sd}\right) \\
\dot i_{gd} &= \frac{\omega_b}{L_t}\left(v_{ed} - v_{gd} - R_t i_{gd} + \omega_{pll} L_t i_{gq}\right) \\
\dot i_{gq} &= \frac{\omega_b}{L_t}\left(v_{eq} - v_{gq} - R_t i_{gq} - \omega_{pll} L_t i_{gd}\right) \\
\dot v_{ed} &= \frac{\omega_b}{C_f}\left(i_{sd} - i_{gd} + \omega_{pll} C_f v_{eq}\right) \\
\dot v_{eq} &= \frac{\omega_b}{C_f}\left(i_{sq} - i_{gq} - \omega_{pll} C_f v_{ed}\right) \\
\dot v_{dc} &= \frac{\omega_b}{C_{dc}}\left(i_{dc} - G_{dc} v_{dc} - m_d i_{sd} - m_q i_{sq}\right) \\
\dot i_{dc} &= \frac{1}{T_{dc}}\left(I_{dc,ref}-i_{dc}\right)
\end{aligned}
$$

**DC-voltage / reactive-power / current loops** — a $d$-axis current
reference regulates the DC-link voltage directly (no separate power
loop); the $q$-axis current reference is the reactive-power loop's own
integrator state, fed straight into the current loop:

$$
i_{sd,ref} = K_{pd}(V_{dc,ref}-v_{dc}) + M_d
\qquad
\dot M_d = K_{id}(V_{dc,ref}-v_{dc})
$$

$$
q=-v_{ed}i_{gq}+v_{eq}i_{gd}
\qquad
p=v_{ed}i_{gd}+v_{eq}i_{gq}
\qquad
\dot M_q = K_{iq}(Q_{ref}-q)
$$

$$
\dot M_{CLd}=K_{iCL}(i_{sd,ref}-i_{sd})
\qquad
\dot M_{CLq}=K_{iCL}(M_q-i_{sq})
$$

**SRF-PLL** — the integral-of-$v_{eq}$ term drives $v_{eq}\to 0$, locking
$\theta_{pll}$ to the terminal-voltage vector angle:

$$
\dot M_{pll} = K_{i,pll}\,v_{eq}
\qquad
\dot\theta_{pll} = \omega_b\,\omega_{pll}
$$

**Algebraic**:

$$
m_d = \frac{1}{v_{dc}}\Big(K_{pCL}(i_{sd,ref}-i_{sd}) + K_{ffv}v_{ed} - \omega_{ff}L_f i_{sq} + M_{CLd}\Big)
$$
$$
m_q = \frac{1}{v_{dc}}\Big(K_{pCL}(M_q-i_{sq}) + K_{ffv}v_{eq} + \omega_{ff}L_f i_{sd} + M_{CLq}\Big)
$$
$$
\omega_{pll} = M_{pll} + K_{p,pll}v_{eq} + \omega_{ff}
$$

**Outputs**: $p$, $q$, $\omega_{pll}$ (recomputed at output as
$M_{pll}+K_{p,pll}v_{eq}+\omega_{ff}$ — algebraically identical to the
constraint above, kept as a separate output expression in the source),
$V_t$, and grid current rotated back to the global frame using
$\theta_{pll}$.

### Infinite bus (`ib.py`)

The simplest DER: an ideal voltage source of fixed magnitude $V_{up}$ and
frequency $\omega_{up}$ behind its own series impedance $R_{up}/L_{up}$
(the DER's own transformer) — 3 states, ported from `Functions/symIB.m`.
No control loops, no algebraic variables at all. Always the slack (never
wired as a non-slack unit — see {doc}`interconnect <interconnect>`), so
its own frame *is* the global frame and its equations never rotate
anything.

**States** $x=[i_{gd,g}, i_{gq,g}, \theta_{up}]$. **Inputs**:
$u_s=[\omega_{up}, V_{up}]$, $u_g=[v_{gd,g}, v_{gq,g}]$.

$$
\dot i_{gd,g} = \frac{\omega_b}{L_{up}}\left(V_{up} - v_{gd,g} - R_{up}i_{gd,g} + \omega_{up}L_{up}i_{gq,g}\right)
$$
$$
\dot i_{gq,g} = \frac{\omega_b}{L_{up}}\left(-v_{gq,g} - R_{up}i_{gq,g} - \omega_{up}L_{up}i_{gd,g}\right)
$$
$$
\dot\theta_{up} = \omega_b\,\omega_{up}
$$

**Outputs**: own $[P_{up}, Q_{up}] = [V_{up}i_{gd,g},\; -V_{up}i_{gq,g}]$;
grid-facing $[i_{gd,g}, i_{gq,g}, \theta_{up}, \omega_{up}]$ (echoed
directly — an IB defines the global frame, so nothing needs rotating).

### Transmission line (`line.py`)

A series $R_l$–$L_l$ branch between two buses $j$ (`from_bus`) and $k$
(`to_bus`), living entirely in the global dq frame (no rotation, unlike
every DER above) — 2 states, ported from `Functions/symLine.m`.

**States** $x=[i_{ld,g}, i_{lq,g}]$. **Inputs**:
$[\omega_g, v_{gdj,g}, v_{gqj,g}, v_{gdk,g}, v_{gqk,g}]$.

$$
\dot i_{ld,g} = \frac{\omega_b}{L_l}\left(v_{gdj,g}-v_{gdk,g} - R_l i_{ld,g} + \omega_g L_l i_{lq,g}\right)
$$
$$
\dot i_{lq,g} = \frac{\omega_b}{L_l}\left(v_{gqj,g}-v_{gqk,g} - R_l i_{lq,g} - \omega_g L_l i_{ld,g}\right)
$$

**Output**: $[i_{ld,g}, i_{lq,g}]$ itself (this current is what
{doc}`interconnect <interconnect>`'s node current-balance wiring sums
into the two endpoint buses, $+$ at $k$, $-$ at $j$). $L_l$ plays the
role of the line's reactance $x_{pu}$ directly, per-unit convention
throughout this codebase.

### Network node (`node.py`)

The shunt line-charging capacitance $C_l$ seen at one bus — its "voltage"
state *is* that bus's own dq voltage, so this component is what actually
defines bus voltage dynamically in the closed-loop system (every other
component treats its own bus voltage as an algebraic/input quantity) — 2
states, ported from `Functions/symNode.m`.

**States** $x=[v_{gd,g}, v_{gq,g}]$. **Input**: $[\omega_g, i_{shd,g},
i_{shq,g}]$ — the net shunt current injected at this bus (KCL sum of every
incident DER/line/load current, built by {doc}`interconnect
<interconnect>`'s wiring rules).

$$
\dot v_{gd,g} = \frac{\omega_b}{C_l}\left(i_{shd,g} + \omega_g C_l v_{gq,g}\right)
\qquad
\dot v_{gq,g} = \frac{\omega_b}{C_l}\left(i_{shq,g} - \omega_g C_l v_{gd,g}\right)
$$

**Output**: $[v_{gd,g}, v_{gq,g}]$ itself — every other component's
`vgd_g`/`vgq_g` input is wired to this output.

### Constant-impedance load (`load.py`)

A series $R_c$–$L_c$ branch to (implicitly) a zero-voltage reference,
i.e. a constant-impedance load equivalent at its bus — 2 states, ported
from `Functions/symLoad.m`. $R_c$/$L_c$ (`load_rx` in
{doc}`operating_point <operating_point>`) are derived from the load's
$P$/$Q$/$V$ at the power-flow operating point, not free parameters.

**States** $x=[i_{cd,g}, i_{cq,g}]$. **Input**: $[\omega_g, v_{gd,g},
v_{gq,g}]$.

$$
\dot i_{cd,g} = \frac{\omega_b}{L_c}\left(v_{gd,g} - R_c i_{cd,g} + \omega_g L_c i_{cq,g}\right)
\qquad
\dot i_{cq,g} = \frac{\omega_b}{L_c}\left(v_{gq,g} - R_c i_{cq,g} - \omega_g L_c i_{cd,g}\right)
$$

**Output**: $[i_{cd,g}, i_{cq,g}]$ — subtracted (load convention) at its
bus's current balance, the same way a line's current is added/subtracted
at its two endpoints.

### Linearization, in one place

Every component above shares one linearization step
({func}`~g2elin_core.components.base.ComponentDAE.linearize`), algebraic-
variable elimination by direct substitution into the Jacobians formed
once by `build_dae()`:

$$
F_x = \frac{\partial f}{\partial x} \quad F_z = \frac{\partial f}{\partial z} \quad F_u = \frac{\partial f}{\partial u}
\qquad
G_x = \frac{\partial g}{\partial x} \quad G_z = \frac{\partial g}{\partial z} \quad G_u = \frac{\partial g}{\partial u}
\qquad
H_x = \frac{\partial h}{\partial x} \quad H_z = \frac{\partial h}{\partial z} \quad H_u = \frac{\partial h}{\partial u}
$$

$$
A = F_x - F_z G_z^{-1} G_x
\qquad
B = F_u - F_z G_z^{-1} G_u
\qquad
C = H_x - H_z G_z^{-1} G_x
\qquad
D = H_u - H_z G_z^{-1} G_u
$$

evaluated numerically at each component's own operating point (all the
$0$-subscripted values each `*OperatingPoint` class computes). This is
exactly the small-signal linearization every `sym*.m` file in the MATLAB
source performs — see the module's own docstring above for the one
deliberate implementation difference (symbolic-vs-numeric order of
substitution and inversion, same result either way).

## Reference

```{eval-rst}
.. automodule:: g2elin_core.components.base
```

```{eval-rst}
.. automodule:: g2elin_core.components.sm
```

```{eval-rst}
.. automodule:: g2elin_core.components.gfm
```

```{eval-rst}
.. automodule:: g2elin_core.components.gfl
```

```{eval-rst}
.. automodule:: g2elin_core.components.ib
```

```{eval-rst}
.. automodule:: g2elin_core.components.frame
```

```{eval-rst}
.. automodule:: g2elin_core.components.line
```

```{eval-rst}
.. automodule:: g2elin_core.components.node
```

```{eval-rst}
.. automodule:: g2elin_core.components.load
```
