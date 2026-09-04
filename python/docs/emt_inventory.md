# G2ELib_V1.slx inventory (P6 investigation)

Generated with [`tools/slx_inspect.py`](../tools/slx_inspect.py). Regenerate with:

```
python tools/slx_inspect.py inventory
```

## Why this exists

P6 ("detailed EMT models" — switching-level converter physics, abc-frame/
unbalanced faults, distributed line effects; see the README's P6 section
for why that's a separate phase from the EMT time-domain simulation
already implemented) is the one feature from the original migration plan
with **no symbolic MATLAB source** to port for the physics it would add —
that physics exists only as Simulink block diagrams in `G2ELib_V1.slx`.
Before committing to (or scoping down) a translation effort, the diagrams
needed to actually be read. This document is the result of that reading,
not a finished port: it is what a P6 implementation would have to
translate, and the one piece of it that has actually been cross-checked
end-to-end.

The key fact that makes this readable at all: `.slx` is a **zip archive of
XML**, not an opaque binary. `simulink/systems/system_<SID>.xml` holds one
file per subsystem — every block's type, name, and parameters as `<P>`
elements, and every wire as `<Line Src=... Dst=...>` (with `<Branch>` for
fan-out). No MATLAB or Simulink installation is needed to read it, only a zip
reader and an XML parser — see `tools/slx_inspect.py`.

## Block-type histogram (whole library, 1135 blocks / 55 subsystems)

| count | block type | count | block type |
|---:|---|---:|---|
| 159 | From | 8 | Saturate |
| 106 | Gain | 8 | Step |
| 101 | Goto | 7 | BusSelector |
| 98 | Reference | 6 | MultiPortSwitch |
| 84 | Sum | 6 | Ground |
| 80 | Inport | 4 | Sqrt |
| 61 | Constant | 4 | StateSpace |
| 59 | Product | 3 | Switch |
| 54 | SubSystem | 3 | Scope |
| 51 | Outport | 3 | Fcn |
| 48 | PMIOPort | 3 | Trigonometry |
| 35 | Demux | 2 | FirstOrderHold |
| 21 | Terminator | 2 | TransferFcn |
| 20 | Mux | | |
| 18 | Integrator | | |
| 17 | Selector | | |
| 15 | BusCreator | | |
| 14 | TransportDelay | | |
| 13 | PMComponent | | |

**Why this matters for scoping**: only `PMComponent` (13) + `Ground` (6) +
`PMIOPort` (48) — 67 of 1135 blocks (~6%) — are Simscape/SimPowerSystems
*physical* (acausal circuit) blocks. Everything else is ordinary Simulink
signal-flow (gains, sums, products, integrators, transfer functions,
state-space blocks, PLLs built from primitives) with parameter values stored
as plain text in the XML. The physical blocks are concentrated in a few
subsystems (Transformer, Load, Upstream Network) whose underlying physics
(RL/RLC branches) is standard and already implemented in this codebase's EMT
component models — so the *hard* 6% is also the *already-solved* 6%, and the
bulk of the work is patiently transcribing readable control-logic diagrams,
not reverse-engineering opaque physics.

`Reference` (98) blocks are links to library blocks (e.g. MathWorks'
"Synchronous Machine pu Fundamental") rather than inline definitions — those
still need to be resolved against Simulink's own block library semantics,
which is not fully documented in the `.slx` alone. `From`/`Goto` (260
combined) implement tag-based signal routing rather than direct wires, which
the current inventory tool does not yet resolve into direct edges.

## Subsystem hierarchy (from `system_root`)

```
- DC/AC Converters  [system_933, 1 blocks]
  - Voltage Source Converters (VSC)  [system_936, 2 blocks]
    - Grid-Following VSC (GFL)  [system_4388, 7 blocks]
      - Init_GFL_i  [system_8394, 6 blocks]
      - VSC - GFL - Ideal  [system_4391, 146 blocks]
        - Current Loop  [system_5071, 26 blocks]
        - Filter  [system_4473, 13 blocks]
        - P-control  [system_5056, 8 blocks]
        - PLL  [system_5057, 20 blocks]
        - Power Measurement  [system_4660, 26 blocks]
        - Q-control  [system_5059, 6 blocks]
        - Subsystem  [system_4716, 8 blocks]
        - Subsystem1  [system_4725, 8 blocks]
        - VIM SU  [system_5066, 29 blocks]
    - Grid-Forming VSC (GFM)  [system_939, 8 blocks]
      - Init_GFM_i  [system_8399, 6 blocks]
      - VSC - GFM - Generic - Ideal  [system_1571, 244 blocks]
        - Current Loop  [system_1585, 26 blocks]
        - DC Voltage Control   [system_1612, 11 blocks]
        - Damping Resistor  [system_1624, 6 blocks]
          - Subsystem  [system_1628, 5 blocks]
        - Filter  [system_1644, 13 blocks]
        - Matching SU (SI)  [system_1742, 15 blocks]
        - P-f Droop SU (PU)  [system_1762, 13 blocks]
        - PI_dc_voltage controller  [system_2219, 6 blocks]
        - PLL  [system_1778, 20 blocks]
        - Power Measurement  [system_2308, 26 blocks]
        - Q-V Droop MMU (PU)  [system_1846, 13 blocks]
        - Subsystem  [system_1871, 8 blocks]
        - Subsystem1  [system_1880, 8 blocks]
        - Transient Virtual Impedance  [system_1904, 21 blocks]
        - VIM SU  [system_1938, 29 blocks]
        - VSM AVR (PU)  [system_1969, 25 blocks]
          - Park Transform2  [system_7251, 9 blocks]
          - emf_ calculation2  [system_7252, 16 blocks]
        - VSM SU (PU)  [system_2012, 17 blocks]
        - Virtual Impedance  [system_2033, 16 blocks]
        - Voltage Loop  [system_2052, 27 blocks]
        - dVOC - P-f (PU)   [system_2079, 33 blocks]
          - APC  [system_2086, 13 blocks]
          - P-AVR  [system_2113, 8 blocks]
          - RPC  [system_2122, 13 blocks]
- Load Models  [system_670, 2 blocks]
  - Load_i  [system_7422, 14 blocks]
- Miscellaneous Network Elements  [system_223, 5 blocks]
  - Bus_i  [system_7413, 13 blocks]
  - Transformer_i  [system_7415, 20 blocks]
  - Upstream Network  [system_8427, 6 blocks]
- Synchronous Machine Models  [system_215, 7 blocks]
  - Simulink - Synchronous Machine - Full Controls  [system_5216, 81 blocks]
    - First-Order Filter1  [system_5486, 3 blocks]
    - First-Order Filter2  [system_5487, 3 blocks]
    - First-Order Filter3  [system_5488, 3 blocks]
    - Synchronous Machine pu Fundamental  [system_5462, 13 blocks]
```

The GFM subsystem (`system_1571`, 244 blocks) is the largest single unit —
larger than GFL (146) and the synchronous-machine block (81) combined —
because it implements *five* interchangeable outer-loop control types in one
diagram (Droop, VSM, dVOC, Matching, plus their shared inner current/voltage
loops and virtual-impedance/PLL infrastructure), matching `symGFM_types.m`'s
multi-branch structure on the MATLAB side (this codebase currently ports only
the `'Droop'` branch — see [`components/gfm.py`](../src/g2elin_core/components/gfm.py)).

## Worked example: independent cross-validation of `gfm.py`'s Droop Q-V control

`system_1846`, "Q-V Droop MMU (PU)" (13 blocks, nested under GFM →
`VSC - GFM - Generic - Ideal`), computes the reactive-power droop reference
voltage. Its full wiring, read directly from the XML:

| SID | block | params |
|---|---|---|
| 1847 | Inport `Q` | |
| 1848 | Inport `Q_set` | |
| 1849 | Inport `Eg_set_dq0` | |
| 2336 | Reference `State-Space (with initial outputs)1` | low-pass filter: `A=-q_filter_wf, B=1, C=q_filter_wf, Y0=GFM_Qref` |
| 1856 | Sum `Sum` | `+-` : `Q_set - qm` |
| 1852 | Gain `Gain` | `Droop_nq` |
| 1851 | Demux | splits `Eg_set_dq0` into d/q |
| 1853 | Math `Hypot` | `hypot(d, q)` = `\|Eg_set_dq0\|` |
| 1857 | Sum `Sum1` | `++` : `\|Eg_set_dq0\| + (Q_set-qm)*Droop_nq` |
| 1855 | Selector | picks component 3 (the 0-sequence) of `Eg_set_dq0` |
| 1854 | Product | `Eg_set_dq0 * Sum1_out / \|Eg_set_dq0\|` (rescales the vector to the new magnitude) |
| 1859 | Outport `Eg_dq0_ref` | |

Tracing the wires (`Q` → State-Space filter → `qm`; `Sum(Q_set, -qm)` →
`Gain(Droop_nq)` → `Sum1(+, +)` with `hypot(Eg_set_dq0)`) gives exactly:

```
qm      = lowpass(Q, wf=q_filter_wf, y0=GFM_Qref)
ve_ref  = hypot(Eg_set_dq0)
ved_ref = ve_ref + (Q_set - qm) * Droop_nq
```

This is **exactly** the equation already ported in
[`components/gfm.py`](../src/g2elin_core/components/gfm.py) from
`Functions/symGFM_types.m`'s `'Droop'` case
(`ved_ref = ve_ref + (q_ref - qm)*nq`, `symSG`-style low-pass on `qm`). The
two were derived independently — one from the MATLAB symbolic source, one by
hand-tracing this XML — and agree exactly. This is the first place in the
whole migration where the ported physics has been checked against the
*actual Simulink diagram* rather than only against the MATLAB script that
(per the project's own account) built it.

## What this does and doesn't establish

**Established**: the `.slx` is fully readable without MATLAB; the control
logic is mostly ordinary signal-flow blocks with literal parameters, not
opaque compiled physics; the diagrams and the already-ported nonlinear
component equations agree, at least for one subsystem.

**Not established, and not attempted here**: a full port. Translating all
1135 blocks across 55 subsystems into verified Python EMT code would still
require (a) resolving `Reference`/library-linked blocks (a MathWorks
Synchronous Machine block among them) against Simulink's own semantics, (b)
resolving `From`/`Goto` tag-routed signals (260 blocks) into explicit wiring,
(c) modeling the Simscape/SimPowerSystems acausal network blocks (13
`PMComponent` + 6 `Ground` + 48 `PMIOPort`) as differential-algebraic
equations rather than signal flow, and (d) some way to validate the result's
*behavior*, not just its structure — there is no MATLAB/Simulink installation
available in this environment to run the original model and diff against.
None of that is done. This document is a scoping and verification artifact,
not a partial implementation.
