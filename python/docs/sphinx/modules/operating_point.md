# `pu_base` and `operating_point`

Two small modules that sit between the solved power flow and the symbolic
{doc}`component models <components>`: `pu_base` computes the physical base
values a per-unit system needs, and `operating_point` turns bus voltages/
angles/powers into the exact operating-point arguments each component's
`*OperatingPoint` class expects.

## `pu_base` — amplitude-invariant per-unit bases

Ported from `Functions/PU_calc.m`. Note the AC bases use an
**amplitude-invariant** (peak-quantity) convention — `Ub = sqrt(2/3)*Ul`,
`Ib = sqrt(2)*Il` — rather than the more common RMS convention, matching the
EMT-oriented per-unit system the MATLAB toolbox uses throughout. The
resulting impedance base still reduces to the familiar
`Un_kV**2 / Sn_MVA`, which is what actually matters for the R/X/B
per-unit conversions used elsewhere in this package.

```{mermaid}
flowchart LR
    IN["wb, un_kv, sn_mva"] --> FN["pu_base()"]
    FN --> OUT["PuBase<br/>ub, ib, zb, lb, cb<br/>ub_dc, ib_dc, zb_dc, cb_dc"]
```

## `operating_point` — from `PowerFlowResult` to per-component operating points

Ports `script_generic.m`'s "SM Initializations" section (computing `SM.Pt`,
`SM.Et`, flux linkages, etc. from `bus_sol`) plus the node/line/load
operating-point loops in its "Substituting in the individual state spaces"
section. This is also where two MATLAB fidelity quirks are implemented
deliberately (not accidentally): **every DER unit's transformer impedance
comes from DG#1's row** (`first_tr`, used for every unit rather than each
unit's own transformer), and **every node's shunt susceptance comes from
line #1's `b_pu`** (handled downstream in the node component, but the values
this module hands it are what make that visible).

```{mermaid}
flowchart TD
    N["Network"] --> CO["compute_operating_point()"]
    PF["PowerFlowResult"] --> CO
    CO --> SLACK["1. Slack SM first<br/>defines theta_g_rad (global angle ref)"]
    SLACK --> SM["2. Remaining SM units<br/>-> SmOperatingPoint"]
    SLACK --> GFMGFL["3. GFM / GFL units<br/>-> GfmOperatingPoint / GflOperatingPoint"]
    SLACK --> NODES["4. Plain nodes<br/>-> node_vg (rotated to global frame)"]
    SLACK --> LINES["5. Lines<br/>solve steady-state 2x2 -> line_i0"]
    SLACK --> LOADS["6. Loads<br/>constant-impedance equivalent -> load_rx"]
    SM --> RESULT["NetworkOperatingPoint"]
    GFMGFL --> RESULT
    NODES --> RESULT
    LINES --> RESULT
    LOADS --> RESULT
```

`sm_params()`, `gfm_params()` and `gfl_params()` hold the default
electrical/control/AVR/PSS parameter values transcribed directly from
`script_generic.m`'s "Defining Default Tuning" section — the same values
every preset case in the MATLAB source uses (none of the example scripts
override them), so hard-coding them here is a faithful port, not a
simplification. One exception: `gfm_params()`'s `tr_cl` (the inner current
loop's closed-loop rise-time tuning target, `Inner_current_loop.tr` in the
MATLAB source) is an optional keyword, default unchanged from the fixed
value every preset used before — added so a caller can retune just that
one loop and trace how it moves the closed-loop eigenvalues, without
duplicating the pole-placement formula outside this function (see
`notebooks/tour.ipynb`'s root-locus section).

## Physics

### Per-unit bases

Given $\omega_b = 2\pi f_n$, a bus's nominal line-to-line voltage
$U_{n}$ (kV) and the network's own MVA base $S_n$:

$$
U_b = \sqrt{\tfrac{2}{3}}\,U_n
\qquad
I_l = \frac{S_n}{\sqrt3\,U_n}
\qquad
I_b = \sqrt2\,I_l
\qquad
Z_b = \frac{U_b}{I_b}
\qquad
L_b = \frac{Z_b}{\omega_b}
\qquad
C_b = \frac{1}{Z_b\,\omega_b}
$$

The $\sqrt{2/3}$/$\sqrt2$ factors make this an **amplitude-invariant**
(peak-quantity) base — the dq-frame convention used throughout every
component's equations — rather than the more common RMS convention; the
impedance base $Z_b$ still reduces to the familiar $U_{n}^2/S_n$ (the
amplitude factors cancel), which is what the R/X/B per-unit values on
every `Line`/`Transformer` are already expressed against. A converter's
DC-link base scales off the same AC quantities (a two-level VSC's DC bus
voltage is nominally $2U_b$ at unity modulation):

$$
U_{b,dc} = 2\,U_b
\qquad
I_{b,dc} = \tfrac34\,I_b
\qquad
Z_{b,dc} = \tfrac83\,Z_b
\qquad
C_{b,dc} = \tfrac38\,C_b
$$

### Default synchronous-machine parameters

Transcribed from `script_generic.m`'s "Defining Default Tuning" section —
every WSCC/CIGRE preset in this codebase uses these values unmodified.

| Electrical (pu) | Value | | Electrical (pu) | Value |
|---|---|---|---|---|
| $R_a$ | 0.003 | | $L_{1q}$ | 0.7252 |
| $L_l$ | 0.15 | | $R_{1q}$ | 0.00619 |
| $L_{ad}$ | 1.66 | | $L_{2q}$ | 0.125 |
| $L_{aq}$ | 1.61 | | $R_{2q}$ | 0.02368 |
| $L_{fd}$ | 0.165 | | $R_g$ (aux load, $=S_n/P_L$) | 1000 |
| $R_{fd}$ | 0.0006 | | | |
| $L_{1d}$ | 0.1713 | | | |
| $R_{1d}$ | 0.0284 | | | |

| Mechanical | Value | AVR | Value | PSS | Value |
|---|---|---|---|---|---|
| $H$ (s) | 5.0 | $T_r$ (s) | 0.02 | $T_{LP}$ (s) | 0.03 |
| $K_D$ | 0.0 | $K_a$ | 300.0 | $K_{PSS}$ | 2.0 |
| $m_p$ (droop, pu) | 0.005 | $T_a$ (s) | 0.001 | $T_{HP}$ (s) | 2.0 |
| $T_G$ (s) | 0.2 | $K_e$ | 1.0 | $T_{1n}/T_{1d}$ (s) | 0.05 / 0.02 |
| | | $T_e$ (s) | 0.0001 | $T_{2n}/T_{2d}$ (s) | 3.0 / 5.4 |
| | | $K_{fd}$ | 0.001 | | |
| | | $T_{fd}$ (s) | 0.1 | | |

### GFM/GFL loop-tuning: rise-time pole placement

Every PI loop in `gfm_params()`/`gfl_params()` is tuned the same way
(`script_generic.m`'s "Inner\_current\_loop"/"Outer\_..." blocks): treat
the loop's own plant as first-order ($\frac{1}{s\tau + 1}$-shaped after
normalizing by its own DC gain), pick a target rise time $t_r$ and a
fixed damping ratio $\zeta = 0.707$, and back out the natural frequency
and PI gains that place the closed loop there:

$$
\omega_n = \frac{3}{\zeta\, t_r}
\qquad
K_p = \frac{2\zeta\omega_n \tau - 1}{g}
\qquad
K_i = \frac{\omega_n^2 \tau}{g}
$$

where $\tau$ is the plant's own time constant and $g$ its DC gain. Every
loop below is this same formula with different $(\tau, g, t_r)$:

| Loop | Plant $\tau$ | Plant gain $g$ | $t_r$ | Result |
|---|---|---|---|---|
| GFM current loop | $L_1/(\omega_b R_1)$ | $1/R_1$ | 0.1 ms | $K_{pCL}, K_{iCL}$ |
| GFM voltage loop | — ($g=\omega_b/C_f$, integrator plant) | $\omega_b/C_f$ | 15 ms | $K_{pVL}, K_{iVL}$ |
| GFM DC-link voltage | $C_{dc}/(\omega_b G_{dc})$ | $1/G_{dc}$ | 5 ms | $K_{pdc}$ (I-only: $K_{pdc}=\left(\tfrac{3\tau_{dc}}{t_{r}}-1\right)/g$) |
| GFL current loop | $L_1/(\omega_b R_1)$ | $1/R_1$ | 10 ms | $K_{pCL}, K_{iCL}$ |
| GFL DC-link voltage | $C_{dc}/(\omega_b G_{dc})$ | $-1/G_{dc}$ | 100 ms | $K_{pd}, K_{id}$ |
| GFL reactive-power loop | — (I-only) | — | 100 ms | $K_{iq} = -3/t_r$ |
| GFL PLL | — ($g=\omega_b$) | $\omega_b$ | 50 ms | $K_{p,pll}, K_{i,pll}$ |

The droop coefficients themselves are fixed, not pole-placed: $m_p =
1/K_{D,opt}$ with $K_{D,opt}=200$ (so $m_p = 0.005$ pu, matching the SM's
own $m_p$ above — deliberately, so a GFM's droop response is comparable
to an SM's governor droop at the same operating point), $n_q =
10^{-4}$ pu, and the power-measurement filter $\omega_f =
\frac{1}{2 m_p H_{1st}}$ with a first-order-equivalent inertia
$H_{1st}=3.0\,\mathrm s$.

### From a solved power flow to a DER's own operating point

Given a DER's terminal voltage phasor $\bar V_0 = V_t\,e^{j\theta_t}$ (from
the power-flow bus solution) and its gross injected complex power $\bar
S_0 = P_0 + jQ_0$ (net power-flow injection plus the unit's own
`p_cons_mw`/`q_cons_mvar`), every operating-point class starts from the
same terminal current:

$$
\bar I_0 = \left(\frac{\bar S_0}{\bar V_0}\right)^{\!*} = \frac{P_0 - jQ_0}{\bar V_0^{\,*}}
$$

**Synchronous machine** — the classical "$E_q$ behind $X_q$" angle-finding
step locates the rotor's own $d$-axis before any dq rotation is possible,
using the machine's own stator resistance $R_a$ and $q$-axis synchronous
reactance $L_q = L_l + L_{aq}$:

$$
\delta_0 = \angle\!\left(\bar V_0 + (R_a + jL_q)\bar I_0\right)
\qquad
\theta_0 = \delta_0 - \frac{\pi}{2}
$$

Rotating $\bar V_0$, $\bar I_0$, and the current split off the machine's
own small aux load ($\bar I_{c,0}=\bar V_0/R_g$, $\bar I_{g,0} = \bar
I_0-\bar I_{c,0}$) into the machine's own dq frame by $e^{-j\theta_0}$
gives $v_{ed0}, v_{eq0}, i_{d0}, i_{q0}, i_{gd0}, i_{gq0}$ directly, and
the flux linkages/field current follow from the stator voltage equations
at steady state ($\dot\psi \equiv 0$):

$$
\psi_{d0} = v_{eq0}+R_a i_{q0}
\qquad
\psi_{q0} = -v_{ed0}-R_a i_{d0}
\qquad
i_{fd0} = \frac{v_{eq0}+R_a i_{q0}+(L_{ad}+L_l)i_{d0}}{L_{ad}}
$$

$$
e_{fd0}=L_{ad}\,i_{fd0}
\qquad
\psi_{fd0}=(L_{ad}+L_{fd})i_{fd0}-L_{ad}i_{d0}
\qquad
\psi_{1d0}=L_{ad}(i_{fd0}-i_{d0})
\qquad
\psi_{1q0}=\psi_{2q0}=-L_{aq}i_{q0}
$$

`SmOperatingPoint` cross-checks $i_{d0}$/$i_{q0}$/$i_{fd0}$ against an
*independent* second derivation — solving the same flux/current matrix
relation {doc}`components <components>`'s algebraic equations use — and
raises immediately if they disagree by more than float tolerance; this
is what caught the double-rotation bug documented in
{doc}`components <components>`'s own page.

**GFM/GFL** — no rotor angle to locate first; the unit's own frame angle
*is* the terminal-voltage angle directly, $\theta_0=\theta_t$. The filter
capacitor's charging current and the converter-side (behind $R_f/L_f$)
voltage follow directly:

$$
\bar E_{g,0}=V_t e^{j\theta_0}
\qquad
\bar I_{g,0}=\left(\frac{\bar S_0}{\bar E_{g,0}}\right)^{\!*}
\qquad
\bar I_{c,0}=j\omega C_f\,\bar E_{g,0}
\qquad
\bar I_{s,0}=\bar I_{g,0}+\bar I_{c,0}
\qquad
\bar V_{m,0}=\bar E_{g,0}+\bar I_{s,0}(R_f+jL_f)
$$

then every dq quantity is $e^{-j\theta_0}$ of its own phasor above,
exactly the same physics for both unit types (a GFL's PI integrator
states $M_d, M_q, M_{CLd}, M_{CLq}, M_{pll}$ are then back-solved from
each loop's own steady-state condition $\dot M \equiv 0$, since
`GFL_subs.m` reads them from upstream data columns this port has no
equivalent source for — see `GflOperatingPoint`'s own docstring).

**Grid-side node/line/load operating points** are simpler still: a plain
bus's dq voltage is just its power-flow $(V,\theta)$ rotated into the
*global* frame by $e^{-j\theta_g}$ (no per-unit rotor angle involved); a
line's current solves its own steady-state $2\times2$ Ohm's-law system
directly from its two end voltages; a load's equivalent series
impedance is $Z=V^2/S$ at angle $\arccos(P/S)$ (inductive convention,
$Q\ge0$).

## Reference

```{eval-rst}
.. automodule:: g2elin_core.pu_base
```

```{eval-rst}
.. automodule:: g2elin_core.operating_point
```
