# `pipeline` and `modal` — linear small-signal analysis end to end

`pipeline.linearize_network()` is the one function that chains everything
in {doc}`components <components>`, {doc}`operating_point <operating_point>`
and {doc}`interconnect <interconnect>` into a single closed-loop
state-space; `modal.analyze()` turns that state-space into eigenvalues,
damping/frequency, and participation factors.

## `linearize_network()`

```{mermaid}
flowchart TD
    N["Network"] --> OP["compute_operating_point()"]
    PF["PowerFlowResult"] --> OP
    OP --> SM["linearize_sm() per SM unit"]
    OP --> GFM["linearize_gfm() per GFM unit"]
    OP --> GFL["linearize_gfl() per GFL unit"]
    OP --> IB["linearize_ib() for an IB slack"]
    OP --> NODE["linearize_node() per plain bus<br/>(uses line #1's b_pu — MATLAB quirk)"]
    OP --> LINE["linearize_line() per line"]
    OP --> LOAD["linearize_load() per load"]
    SM --> AN["assemble_network()"]
    GFM --> AN
    GFL --> AN
    IB --> AN
    NODE --> AN
    LINE --> AN
    LOAD --> AN
    AN --> SYS["AssembledSystem<br/>(A, B, C, D)"]
```

Raises `NotImplementedError` if any DER unit isn't SM/GFM/GFL/IB, and the
slack isn't SM or IB specifically ({doc}`components.ib <components>` is
always the slack — see its own docstring). Wiring the infinite bus in
(`interconnect.assemble`'s `"ib_slack"` port kind,
`interconnect.network_assembly`'s `_SLACK_KIND_BY_UNIT_TYPE`) was built
alongside the SMIB presets ({doc}`network <network>`) — it existed
earlier but had no real preset to validate against until then.

## `modal.analyze()`

Ported from the numeric core of `Functions/modal_analysis.m` — eigenvalues,
the damping/frequency table, and participation factors.

```{mermaid}
flowchart LR
    A["AssembledSystem.A"] --> EIG["scipy.linalg.eig(A, left=True, right=True)<br/>(single call — right & left eigenvectors<br/>come back consistently paired)"]
    EIG --> SORT["sort by ascending real part"]
    SORT --> NORM["normalize: W @ V = I"]
    NORM --> PART["participation[i,j] =<br/>|V[i,j]|*|W[j,i]| / sum_k(|V[k,j]|*|W[j,k]|)"]
    PART --> RES["ModalAnalysisResult<br/>eigenvalues, eigenvectors, participation"]
    RES --> TAB["summary_table()<br/>freq/damping + top-3 participating states per mode"]
```

**Bug fixed here, not incidental**: an earlier version computed right
eigenvectors from `eig(A)` and left eigenvectors from a *separate*
`eig(A.conj().T)` call, assuming the two came back in matching order —
nothing guarantees that, and mismatched pairing silently produced
near-zero or NaN participation factors on cases with closely-spaced
eigenvalues (first caught on the CIGRE preset, which has two
near-identical GFM units). The current single-call `scipy.linalg.eig(A,
left=True, right=True)` sidesteps the whole class of bug.

## `modal.toolbox` — the rest of `modal_analysis.m`

`modal/analysis.py`'s own docstring used to flag its plotting and
sensitivity-tensor features as "not ported yet." Read in full, that MATLAB
file's toolbox turns out to be four more pieces of pure linear algebra on
the same eigendecomposition `analyze()` already computes (plus one thin
`scipy.signal` wrapper) — nothing here is a new model or a new derivation:

```{mermaid}
flowchart TD
    MODAL["ModalAnalysisResult<br/>(eigenvalues, left/right eigenvectors, participation)"] --> SENS["eigenvalue_sensitivity(mode)<br/>d(lambda_i)/d(A[k,j]) = w_i[k]*v_i[j]<br/>masked to A's structurally-nonzero entries"]
    MODAL --> SHAPE["mode_shape(mode)<br/>top-5 states' right-eigenvector phase angle<br/>(magnitude fixed at 1 -- a shape plot, not a magnitude one)"]
    MODAL --> FREE["free_response(state, offset, t)<br/>x(t) = V @ (c * exp(eig*t)), c = W @ x0<br/>closed-form -- no ODE solve"]
    SYS["AssembledSystem<br/>(A, B, C, D)"] --> STEP["step_response(input, output, amplitude, t)<br/>SISO reduction + scipy.signal.step"]
```

Each is cross-checked against an independent reference in
`tests/test_modal_toolbox.py`, not just "runs and is finite": `free_response`
against `scipy.linalg.expm(A*t) @ x0` directly, `step_response` against
`scipy.signal.step` on the same SISO-reduced system built by hand.

## Physics: modal analysis

Every mode of the closed-loop linear system $\dot x = A_{tot}x$ is an
eigenpair $(\lambda_i, v_i)$ of $A_{tot}$: $A_{tot}v_i = \lambda_i v_i$.
Complex $A_{tot}$ eigenvalues always come in conjugate pairs
$\lambda=\sigma\pm j\omega_d$ for a real system, each pair describing one
oscillatory mode:

$$
f_n = \frac{|\lambda|}{2\pi}
\qquad
\zeta = \frac{-\operatorname{Re}(\lambda)}{|\lambda|}
\qquad
f_d = f_n\sqrt{1-\zeta^2}
$$

— undamped natural frequency, damping ratio, and damped (actually
observed oscillation) frequency; $\zeta>0$ (equivalently
$\operatorname{Re}(\lambda)<0$) is the stability condition for that mode.
`analyze()` computes both the right eigenvectors $V$ ($A_{tot}V = V\Lambda$)
and left eigenvectors $W$ ($WA_{tot}=\Lambda W$) from one
`scipy.linalg.eig(..., left=True, right=True)` call (see the module's own
"bug fixed here" note below for why not two separate calls), normalized
so $WV=I$. The **participation factor** of state $i$ in mode $j$ —
Kundur's standard measure of how much a given state actually
participates in a given oscillatory mode, dimensionless and normalized
so each mode's column sums to 1 — is

$$
p_{ij} = \frac{|V_{ij}|\,|W_{ji}|}{\sum_k |V_{kj}|\,|W_{jk}|}
$$

### `modal.toolbox` — sensitivity, mode shape, free/step response

**Eigenvalue sensitivity** to a structural $A_{tot}$ entry (Kundur eq.
12.75-family result): perturbing element $A_{kj}$ shifts eigenvalue
$\lambda_i$ by

$$
\frac{\partial \lambda_i}{\partial A_{kj}} = w_i[k]\,v_i[j]
$$

masked to only the entries where $A_{tot}$ is structurally nonzero (a
sensitivity value at a zero entry isn't a sensitivity to any real
physical parameter).

**Mode shape**: the relative *phase* of the top-participating states'
own right-eigenvector components for mode $i$ — magnitude fixed at 1 by
construction, since this is a shape (relative timing/direction of
oscillation between states) not a magnitude comparison:

$$
\text{shape}_k = \angle V_{k,i}\,, \quad k \in \text{top participating states}
$$

**Free (initial-condition) response** — a closed-form modal expansion,
no ODE integration needed, valid because the system is linear:

$$
x(t) = V\,\big(c \odot e^{\Lambda t}\big)\,,\qquad c = W x_0
$$

($x_0$ a unit perturbation at one chosen state, $\odot$ elementwise
product against each mode's own $e^{\lambda_i t}$ — this is literally
the general solution of a linear ODE in modal coordinates).

**Step response** between one chosen exogenous input and output is a
SISO reduction of the full MIMO $(A,B,C,D)$ — pick input $i$'s column of
$B$ and output $j$'s row of $C$/$D$ — then `scipy.signal.step` on that
reduced system, exact (not approximate) since the system is linear and
the reduction is just picking out one input/output pair:

$$
\dot x = A x + B_{:,i}\,u \qquad y = C_{j,:}x + D_{j,i}\,u
$$

## Reference

```{eval-rst}
.. automodule:: g2elin_core.pipeline
```

```{eval-rst}
.. automodule:: g2elin_core.modal.analysis
```

```{eval-rst}
.. automodule:: g2elin_core.modal.toolbox
```
