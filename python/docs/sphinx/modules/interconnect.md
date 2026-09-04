# `interconnect` — wiring components into one closed-loop system

Functionally ports `Functions/assoc_matrices.m` plus the concatenation/
closed-loop formula in `script_generic.m`'s Section IX — but instead of
replicating ~300 lines of hand-rolled row/column index arithmetic (verified
once against a MATLAB run, near-impossible to verify blind), this builds
the same F/G/K/L interconnection from an explicit list of named-port wiring
rules, checked against `assoc_matrices.m`'s Parts I–VII by inspection (each
rule is documented in the module's own docstring: reference-angle
propagation, DG terminal voltage, node current balance, line/load voltage
feed).

## `Block` / `Wiring` / `Topology`

```{mermaid}
classDiagram
    class PortSpec {
        <<protocol>>
        n_us n_ug n_out_s n_out_g : int
    }
    class Block {
        name kind : str
        comp : PortSpec
        state_off input_off output_off : int
        input_row(port) int
        output_col(port) int
    }
    class Wiring {
        input_block : Block
        input_port : str
        terms : list~(Block, str, float)~
    }
    class Topology {
        blocks : list~Block~
        n_x n_u n_y : int
        G : ndarray
        F : ndarray
        input_names : list~str~
    }
    class AssembledSystem {
        A B C D : ndarray
        state_names input_names output_names
    }
    Block --> PortSpec
    Wiring --> Block : reads/writes ports of
    Topology --> Block : offsets
    Topology ..> AssembledSystem : assemble() builds this from Topology + block A/B/C/D
```

`PortSpec` is a `Protocol` — just the four port-count attributes
(`n_us`, `n_ug`, `n_out_s`, `n_out_g`). Both `LinearComponent` (from
{doc}`components <components>`) and the nonlinear wrapper in
{doc}`timedomain.emt <timedomain>` satisfy it structurally, which is what
lets `compute_topology()` be shared code between the linear closed-loop
assembly here and the nonlinear coupled-Newton assembly in `timedomain/emt.py`
— topology (*which port connects to which*) is identical either way; only
what runs on top of it (linear algebra vs. a nonlinear solve) differs.

## From wiring rules to a closed-loop state-space

```{mermaid}
flowchart TD
    BW["build_blocks_and_wiring()<br/>(network_assembly.py)"] --> TOPO["compute_topology(blocks, wiring)"]
    TOPO --> G["G: full input = F @ u_exo + G @ y<br/>(interconnection selector matrix)"]
    TOPO --> F["F: selects each block's own (us) inputs"]
    G --> ASM["assemble(blocks, wiring)"]
    F --> ASM
    BLKDIAG["block-diagonal A_ol/B_ol/C_ol/D_ol<br/>(one per-component LinearComponent, concatenated)"] --> ASM
    ASM --> CLOSED["closed-loop elimination<br/>E_ol = inv(I - D_ol@G)<br/>A_tot = A_ol + B_ol@G@E_ol@C_ol<br/>B_tot, C_tot, D_tot similarly"]
    CLOSED --> OUT["AssembledSystem<br/>(A, B, C, D, state/input/output names)"]
```

`network_assembly.build_blocks_and_wiring()` is the network-specific part:
given a `Network` and its already-linearized components, it decides *which*
`Block`s exist (one per SM/GFM/GFL/IB/node/line/load) and *which* `Wiring`
rules connect them, then hands that to the kind-agnostic `assemble()`
above. It raises `NotImplementedError` if the slack unit is anything other
than a synchronous machine or an infinite bus — `_SLACK_KIND_BY_UNIT_TYPE`
picks `"sm_slack"` or `"ib_slack"` accordingly; every other wiring rule in
that module is already generic over "the slack block" by role, not by
kind, so adding the IB slack kind needed no other change there.

**Block names are per-type counters**, not the DER's raw id: `"SM_1"`,
`"GFM_1"`, `"GFM_2"`, `"GFL_1"` rather than `"DG_1"`/`"DG_2"`/`"DG_3"`/
`"DG_4"` — counted in `network.der_units` declaration order (every preset
lists them in ascending id order already). Every state/input/output name
downstream is built by interpolating a block's name (e.g. `f"dw_r_{{{block.name}}}"`
→ `"dw_r_{SM_2}"`), so this one change is what makes a name say what kind
of unit it is without cross-referencing the network definition — a single
surgical change here propagates everywhere (`AssembledSystem`'s and
`NonlinearNetworkModel`'s `state_names`/`input_names`/`output_names`
alike), since both are built by the same `build_blocks_and_wiring()`.

## Physics: closed-loop elimination

Every component contributes its own open-loop linear state-space (its
own $A_i,B_i,C_i,D_i$ from {doc}`components <components>`'s
linearization). Stacking all $n$ blocks block-diagonally gives the
**open-loop** system

$$
A_{ol}=\operatorname{blkdiag}(A_1,\dots,A_n)
\qquad
B_{ol}=\operatorname{blkdiag}(B_1,\dots,B_n)
\qquad
C_{ol}=\operatorname{blkdiag}(C_1,\dots,C_n)
\qquad
D_{ol}=\operatorname{blkdiag}(D_1,\dots,D_n)
$$

with the full input vector $u$ split into each block's own exogenous
inputs $u_{exo}$ (selected by $F$, e.g. an SM's $P_{ref}$) and the
interconnection signals wired between blocks' own outputs $y$ (selected
by $G$, built directly from the wiring rules below):

$$
u = F\,u_{exo} + G\,y
\qquad
y = C_{ol}x + D_{ol}u
$$

Substituting the first equation into the second and solving for $y$
(hence $u$) in terms of $x$ and $u_{exo}$ alone gives the algebraic loop's
closed form — the same elimination {doc}`components <components>` uses
per-component, applied once more at the network level:

$$
E = \left(I - D_{ol}G\right)^{-1}
$$

$$
A_{tot} = A_{ol} + B_{ol}\,G\,E\,C_{ol}
\qquad
B_{tot} = \left(B_{ol}\,G\,E\,D_{ol} + B_{ol}\right)F
\qquad
C_{tot} = L\,E\,C_{ol}
\qquad
D_{tot} = L\,E\,D_{ol}\,F
$$

($L$ selects each block's own named ("$s$") outputs — e.g. an SM's
$P_e$/$Q_e$/$\omega_r$ — out of the full output vector $y$, the same way
$F$ selects named inputs out of $u$.) $A_{tot}$ is exactly what
{doc}`modal analysis <pipeline>` eigen-decomposes.

**The wiring rules themselves** — what actually populates $G$ — are
physical KVL/KCL statements in the dq frame, not free choices:

- **Reference-frame propagation**: every non-slack DER's own $\theta_g$
  input, and every node/line/load's own $\omega_g$ input, equals the
  slack's own $\theta$/$\omega_r$ output — one shared rotating reference
  frame for the whole network, defined by whichever unit is the slack.
- **DER terminal voltage**: a DER's $v_{gd,g}/v_{gq,g}$ input equals its
  own transformer's `hv_bus` node's voltage output — a DER sees exactly
  the voltage of the grid bus its transformer connects to.
- **Node current balance (KCL)**: a node's injected shunt current
  $i_{sh}$ equals the signed sum of every incident branch's current —
  $+1$ for a DER whose transformer lands here, $+1$ for an incoming
  line ($k=$ this bus), $-1$ for an outgoing line ($j=$ this bus), $-1$
  for a load:
  $$
  i_{sh,d} = \sum_{\text{DER here}} i_{gd} \;+\!\!\sum_{\text{lines }k=\text{here}}\!\! i_{ld} \;-\!\!\sum_{\text{lines }j=\text{here}}\!\! i_{ld} \;-\!\!\sum_{\text{loads here}} i_{cd}
  $$
  (and identically for the $q$-axis) — this is the same statement as
  $Y_{bus}$'s row sum in power flow, just carried per-branch instead of
  pre-summed into an admittance matrix, since each branch here is its
  own dynamic state, not a static impedance.
- **Line/load voltage**: a line's two endpoint voltage inputs, and a
  load's own voltage input, equal their respective node's voltage
  output directly (KVL — a branch sees the potential difference of the
  nodes it connects).

## Reference

```{eval-rst}
.. automodule:: g2elin_core.interconnect.assemble
```

```{eval-rst}
.. automodule:: g2elin_core.interconnect.network_assembly
```
