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

## Presets

15 networks in 4 families. The first two are transcribed line-by-line
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
