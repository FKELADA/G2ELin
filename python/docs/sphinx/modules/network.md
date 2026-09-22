# `network` — schema and presets

The typed network model that replaces the MATLAB `Y_network`/`Y_line`/
`Y_TR`/`Y_DER` numeric-column convention with named, validated
[pydantic](https://docs.pydantic.dev/) fields. Every other module in
`g2elin_core` takes a `Network` (or the `PowerFlowResult` computed from one)
as its starting point — see the {doc}`architecture diagram </index>`.

## Schema

```{mermaid}
classDiagram
    class Network {
        +name: str
        +f_hz: float
        +sn_mva: float
        +buses: list~Bus~
        +lines: list~Line~
        +transformers: list~Transformer~
        +loads: list~Load~
        +shunts: list~Shunt~
        +der_units: list~DerUnit~
        +bus(bus_id) Bus
    }
    class Bus {
        +id: int
        +name: str
        +vn_kv: float
    }
    class Line {
        +from_bus: int
        +to_bus: int
        +r_pu: float
        +x_pu: float
        +b_pu: float
    }
    class Transformer {
        +hv_bus: int
        +lv_bus: int
        +r_pu: float
        +x_pu: float
        +sn_mva: float
    }
    class Load {
        +bus: int
        +p_mw: float
        +q_mvar: float
    }
    class DerUnit {
        +bus: int
        +unit_type: UnitType
        +bus_type: BusType
        +v_set_pu: float
        +p_set_mw: float
        +controller: GfmController?
    }
    class BusType {
        <<enumeration>>
        SLACK
        PV
        PQ
    }
    class UnitType {
        <<enumeration>>
        NONE
        INFINITE_BUS
        GFM
        GFL
        SYNCHRONOUS_MACHINE
    }
    class GfmController {
        <<enumeration>>
        DROOP
        DROOP_FILTERED
        DVOC
        VSM
        MATCHING
    }

    Network "1" *-- "many" Bus
    Network "1" *-- "many" Line
    Network "1" *-- "many" Transformer
    Network "1" *-- "many" Load
    Network "1" *-- "many" DerUnit
    DerUnit --> UnitType : unit_type
    DerUnit --> BusType : bus_type
    DerUnit --> GfmController : controller (GFM only)
```

Note that load-flow bus type (`BusType`: slack/PV/PQ) and unit type
(`UnitType`: infinite bus/GFM/GFL/SM/none) live on `DerUnit`, not on `Bus`
itself — matching `Y_DER`'s columns in the original MATLAB convention. A
plain network node with no DER attached is just a `Bus`.

Validators (not shown above — see the source) enforce that every bus
reference in a `Line`/`Transformer`/`Load`/`DerUnit` resolves to a real
`Bus`, and that exactly one bus is `SLACK`. These are structural
guarantees the original MATLAB scripts relied on by convention, not
schema — encoding them here means a malformed network fails fast at
construction time rather than deep inside a solver.

### Transformers: two roles, one element

Whether a transformer is a **unit step-up** or a **branch** is decided by its
LV bus, not by a flag (`breakers.unit_transformers` /
`breakers.branch_transformer_indices`):

| | unit step-up | branch transformer |
| --- | --- | --- |
| LV bus | carries a DER unit | an ordinary grid bus |
| dynamic model | inside the unit's own model, as its `Rt`/`Lt` -- the classic machine behind a transformer impedance | a block of its own, the same RL branch a line uses |
| its LV bus | never appears; it lives inside the unit | an ordinary node |
| `tap_ratio` / `shift_degree` | not modelled (the unit's own ratio) | modelled |
| opening its breaker | takes the unit with it | drops just the branch, as for a line |

A branch transformer used to be absent from the dynamic model entirely while
carrying power in the power flow, which left the operating point the model was
linearized about not being an equilibrium at the two buses it joined -- ~450
pu/s of drift at exactly those buses, with no validation error.

The ratio needs no new component maths. A transformer branch is a series
impedance behind an ideal ratio $a\,e^{j\varphi}$, so the series part sees

$$v' = \frac{v_{hv}}{a\,e^{j\varphi}},\qquad i_{hv} = \frac{i\,e^{j\varphi}}{a}$$

which in dq is a $1/a$ scaling and a rotation by $-\varphi$ -- plain
coefficients on the wiring the interconnection already supports, and the same
ones in the linear and nonlinear paths. Power is conserved across the ideal
part: $v_{hv} i_{hv} = a v' \cdot i/a = v' i$. Two transformers claiming the
same unit's LV bus are rejected, since that would decide by list order which
grid bus the unit injects into.

### Bus shunts

A bus's own dynamic model (`components/node.py`) *is* a capacitance to
ground, $C\,dv/dt = i - j\omega C v$, which makes the two kinds of shunt
behave very differently:

* a **capacitor bank** (`q_mvar < 0`) is simply more of that capacitance. It
  adds no states and no block at all -- its current is already in the node
  equation.
* a **reactor** (`q_mvar > 0`) cannot be folded in: a negative capacitance
  would flip the sign of the node's dynamics, and netting the two against each
  other would only hold at exactly nominal frequency. It gets the line's RL
  branch with its far end wired to an empty sum, which is ground.

Its `q_mvar` is the nameplate ($Q = V^2/X$); the resistance comes from `r_pu`,
or from $X/R = 50$ when that is not given. The power flow is handed the P and
Q of that same R-X pair, so both sides of the tool agree on what the device
draws.

### Each bus's capacitance

`breakers.node_capacitances` gives every bus half the charging of the lines
meeting at it (a pi-model puts half at each end) plus its capacitor banks.

`Network.nodes_share_first_line_b` restores the MATLAB toolbox's convention,
where every bus instead borrows the *first line's* charging whatever is
connected to it. The ported presets set it, so they still reproduce that
tool's numbers. Measured, the borrowed value was 0.95x to 1.61x off across
WSCC's buses and 0.06x to 1.28x across CIGRE's, and it was the whole reason
the operating point was not an equilibrium: the predicted drift
$\omega_b |v| (\text{ratio} - 1)$ matched the measured one at every bus.
Giving each bus its own takes WSCC's drift from 178 to 0.085 pu/s, and makes
EMT *faster* (7.9 s -> 1.1 s for 50 ms), because the solver no longer chases a
transient that only existed because the starting point was wrong.

A bus whose lines declare no charging at all -- distribution-feeder data
routinely omits it -- is refused rather than divided by, naming
`Network.min_node_b_pu` as the way to say what to use instead.

## Presets

21 networks in 5 families. The first two are transcribed line-by-line
from a complete MATLAB script, with exact numeric values and source-line
comments tying each field back to the `.m` file it came from — each
shares one raw topology across every DER-mix variant in its family
(`_wscc9_topology()`/`_cigre_raw_topology()`), with a thin per-variant
function supplying just that variant's own DER rows:

- `wscc9_3sm()` and 7 more DER-mix variants — WSCC 9-bus
  (`Functions/WSCC_raw.m` + `WSCC/script_WSCC.m`'s `switch model.name`,
  one `case` per variant): `wscc9_2sm_1gfl`, `wscc9_1sm_2gfl`,
  `wscc9_1sm_1gfm_1gfl`, `wscc9_1sm_2gfm`, `wscc9_2sm_1gfm`,
  `wscc9_1gfm_2gfl`, `wscc9_3gfm`. The last two have a **GFM as the
  slack** — power flow works (pandapower's slack handling is agnostic to
  unit type), but modal/EMT raise `NotImplementedError` (HTTP 501),
  since {doc}`interconnect <interconnect>` only wires up an SM/IB slack
  so far.
- `cigre_islanded_1sm_2gfm_1gfl()` and 3 more — CIGRE MV benchmark
  feeder, islanded (`Functions/CIGRE_raw.m` +
  `CIGRE/script_CIGRE_Islanded.m`'s own `switch model.name`):
  `cigre_islanded_1sm_1gfm_1gfl`, `cigre_islanded_1sm_3gfm_1gfl`,
  `cigre_islanded_2sm_2gfm_2gfl` (this last one the only preset with two
  synchronous machines). `cigre_islanded_1sm_1gfm_1gfl()` — the smallest
  network here with one of each major DER type — is the subject of
  {doc}`../cigre_walkthrough` and its companion notebook,
  {doc}`../_notebooks/cigre_1sm_1gfm_1gfl_walkthrough`.

Two MATLAB quirks are preserved deliberately rather than silently
"fixed," for fidelity with the original tool's actual (if arguably buggy)
behavior — see the docstrings below for exactly where and why.

The other three are **SMIB** (single-machine-infinite-bus) presets —
`sm_smib()`, `gfm_smib()`, `gfl_smib()` — and are *not* a straight port
the way the two above are: `Functions/preset_networks.m` names exactly
these three cases but its `switch` body for each is empty, so the
topology/line/load values are real (transcribed from `Functions/SMIB_raw.m`,
confirmed to share `CIGRE_raw.m`'s own base-value convention numerically)
while the DER attachment is a documented reconstruction following this
project's own `network_form.m` convention. See `sm_smib()`'s own
docstring for the full account. These are also the first presets to
exercise `g2elin_core.components.ib` (the infinite-bus component, built
early in this project but never wired into the interconnection until
these presets needed a real slack for it — see {doc}`interconnect
<interconnect>`).

### Kundur two-area

`kundur_two_area()` and `kundur_two_area_classic()` are the first presets
that are **not** ports of the MATLAB toolbox: two areas of two 900 MVA
machines joined by a weak 220 km double-circuit tie, from Kundur's *Power
System Stability and Control*, Example 12.6 — the textbook case for
inter-area oscillations and for what a stabilizer is for. They use this
tool's own conventions throughout (each bus its own capacitance, each unit
its own transformer), and they are the first to need three of the things
added for them: unit ratings, capacitor banks, and per-bus capacitance.

The book publishes the machine as `Xd`/`Xq`/`Xl` with transient and
subtransient reactances and open-circuit time constants, while
{doc}`components <components>` is written in mutual and leakage inductances
with explicit field and damper windings. `kundur_machine_params()` converts
between them; `tests/test_kundur.py` inverts the conversion and checks all
eight published quantities come back.

Checked against the book's own operating point: the tie carries 400.4 MW
(published 400), and the rotor angles come out 9.77°, 27.08° and 37.27°
behind G1 against the published 9.7°, 27.0° and 37.2°.

There are two of them because the mode *frequencies* depend on controls this
tool models differently:

| | inter-area | local (area 1 / area 2) |
| --- | --- | --- |
| `kundur_two_area` (this tool's AVR + PSS + governor) | 0.71 Hz, +12.9% | 1.23 / 1.26 Hz |
| `kundur_two_area_classic` (the book's: fast exciter, no PSS, constant torque) | **0.61 Hz, −2.5%** | 1.12 / 1.15 Hz |
| Kundur, published | ~0.55 Hz, negative | ~1.1 Hz |

The `classic` variant reproduces what the example exists to show: without a
stabilizer the inter-area mode is *unstable*. The governor turned out to
move it most — its 0.5% droop stiffens the mode from 0.61 to 0.71 Hz and adds
one of its own near 0.25 Hz — so `classic` sets a droop large enough that it
does not respond, the SM model having no governor off switch. The residual
offset from the published frequency is the exciter and stabilizer structure,
which is this tool's and not the book's; the tests assert the mode structure
and the sign of the damping rather than frequencies the model cannot be
expected to reproduce exactly.

## Topology layout

`network.topology.compute_topology_layout()` computes a 2D graph layout
(Kamada-Kawai, over the bus graph — lines and transformers as edges) for
the web UI's Network and Power Flow tabs. A DER unit or load isn't a
separate geometric node in this schema — a DER sits on its own real `Bus`
(`network_form.m` always gives it one), and a load just references an
existing bus id — so the graph to lay out is exactly `network.buses` +
`network.lines`/`network.transformers`; DER/load presence is metadata
decorating an existing node, not a node of its own.

`BusNode`/`TopologyEdge` also carry each element's electrical/control
parameters, not just layout — a `Line`'s `r_pu`/`x_pu`/`b_pu`/`length_km`,
a `Transformer`'s `r_pu`/`x_pu`/`sn_mva`, and (for a DER-hosting bus) the
`DerUnit`'s own dispatch fields plus, for SM/GFM/GFL, the *derived*
electrical/control-parameter dict from `operating_point.sm_params()`/
`gfm_params()`/`gfl_params()` — the same values those functions hand to
the linearizer, not a duplicate of the static schema. This needs no power
flow (those three functions are pure functions of network-level static
data — `sn_mva`/`f_hz`/the DER's own bus `vn_kv`/the first transformer's
impedance, the same "every unit uses DG#1's transformer" quirk
{doc}`components <components>` documents), so `/topology` stays cheap.
The web UI renders this as a click-to-inspect detail panel on each
bus/line/transformer in the diagram.

## Validation

pydantic's own validators on `Network` guarantee two structural properties
unconditionally (every bus reference resolves, exactly one slack DER).
Every hand-crafted preset in this codebase also satisfies a further
invariant by construction, following `network_form.m`'s convention: each
DER unit sits on its own private bus, reached from the rest of the network
through exactly one transformer, and that bus is never also a Line
endpoint, a Load's bus, or another Transformer's `hv_bus`. A network
assembled by hand (the web UI's network editor or drag-and-drop builder)
can easily violate any part of that — e.g. a DER wired with "Wire: Line"
instead of "Wire: Transformer" — and `operating_point.compute_operating_point()`
plus `interconnect.network_assembly.build_blocks_and_wiring()` both assume
it unconditionally, previously crashing with a bare `KeyError`/`IndexError`
deep inside instead of a clear message when it didn't hold.

`network.validation.validate_network(network) -> list[NetworkIssue]` checks
this (plus duplicate bus/DER ids, self-loop lines/transformers, and overall
graph connectivity) in one pass and returns *every* issue found, not just
the first — each tagged `"error"` (blocks the affected capability) or
`"warning"` (e.g. a GFM/GFL slack: power flow still works, modal/EMT
don't support it yet) and which of `powerflow`/`modal`/`emt` it
affects. `compute_operating_point()` calls this itself and raises one
combined `ValueError` listing every `"error"`-severity issue if any exist
(-> HTTP 422, not a 500) — so this validation runs regardless of whether a
caller checks proactively. The web UI's network editor also calls
`POST /api/network/validate` directly (see {doc}`../api_and_web`) as a
cheap, no-power-flow-needed pre-flight check, showing every issue in the
editor's own panel rather than discovering them one crash at a time.

For all of this exercised hands-on — a `random_network()` builder that
respects every convention above by construction, several deliberately
broken variants showing exactly what `validate_network()` catches
(including cases these unconditional pydantic validators alone can't,
like a DER wired to another DER's own bus), and a 30-seed batch proving
the builder is robust rather than one lucky draw — see
{doc}`../_notebooks/random_network`.

## Reference

```{eval-rst}
.. automodule:: g2elin_core.network.schema
```

```{eval-rst}
.. automodule:: g2elin_core.network.presets
```

```{eval-rst}
.. automodule:: g2elin_core.network.topology
```

```{eval-rst}
.. automodule:: g2elin_core.network.validation
```
