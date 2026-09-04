# P6: detailed EMT models — planned, not started

{doc}`P4's EMT simulation <modules/timedomain>` already integrates the
network's actual nonlinear time-domain equations — this phase is about
specific higher-fidelity physics that formulation structurally can't
represent: switching-level converter behavior (PWM/modulation, not an
averaged control law), a three-phase/abc-frame or explicitly unbalanced
representation (faults, negative-sequence currents), and distributed/
traveling-wave line models. See the README's P6 section for the full
rationale for why this is a separate, later phase rather than a
prerequisite for calling P4 "EMT simulation."

Unlike every other feature in this project, that physics has **no symbolic
MATLAB source to port** — it exists only as Simulink block diagrams in
`G2ELib_V1.slx`, with no `sym*.m`-equivalent nonlinear equations anywhere
in the MATLAB repo. This page documents what reading that `.slx` file
directly has established — real investigative groundwork, not an
implementation; nothing here has been built yet. The full source document
is included below; it was generated with the dev-time inspector at
`python/tools/slx_inspect.py`.

## Subsystem hierarchy (from `system_root`)

```{mermaid}
graph TD
    ROOT["system_root"]
    ROOT --> DCAC["DC/AC Converters"]
    DCAC --> VSC["Voltage Source Converters (VSC)"]
    VSC --> GFL["Grid-Following VSC (GFL)"]
    VSC --> GFM["Grid-Forming VSC (GFM)"]

    GFL --> GFLINIT["Init_GFL_i"]
    GFL --> GFLIDEAL["VSC - GFL - Ideal<br/>146 blocks"]
    GFLIDEAL --> GFLCL["Current Loop"]
    GFLIDEAL --> GFLFILT["Filter"]
    GFLIDEAL --> GFLPC["P-control"]
    GFLIDEAL --> GFLPLL["PLL"]
    GFLIDEAL --> GFLPM["Power Measurement"]
    GFLIDEAL --> GFLQC["Q-control"]
    GFLIDEAL --> GFLVIM["VIM SU"]

    GFM --> GFMINIT["Init_GFM_i"]
    GFM --> GFMIDEAL["VSC - GFM - Generic - Ideal<br/>244 blocks (largest single unit)"]
    GFMIDEAL --> GFMCL["Current Loop"]
    GFMIDEAL --> GFMDCV["DC Voltage Control"]
    GFMIDEAL --> GFMDR["Damping Resistor"]
    GFMIDEAL --> GFMFILT["Filter"]
    GFMIDEAL --> GFMMATCH["Matching SU"]
    GFMIDEAL --> GFMPFDROOP["P-f Droop SU"]
    GFMIDEAL --> GFMPI["PI_dc_voltage controller"]
    GFMIDEAL --> GFMPLL["PLL"]
    GFMIDEAL --> GFMPM["Power Measurement"]
    GFMIDEAL --> GFMQV["Q-V Droop MMU<br/>verified against gfm.py"]
    GFMIDEAL --> GFMTVI["Transient Virtual Impedance"]
    GFMIDEAL --> GFMVIM["VIM SU"]
    GFMIDEAL --> GFMAVR["VSM AVR"]
    GFMIDEAL --> GFMVSMSU["VSM SU"]
    GFMIDEAL --> GFMVI["Virtual Impedance"]
    GFMIDEAL --> GFMVL["Voltage Loop"]
    GFMIDEAL --> GFMDVOC["dVOC - P-f"]

    ROOT --> LM["Load Models"]
    LM --> LOADI["Load_i"]
    ROOT --> MNE["Misc Network Elements"]
    MNE --> BUSI["Bus_i"]
    MNE --> TRI["Transformer_i"]
    MNE --> UPN["Upstream Network"]
    ROOT --> SMM["Synchronous Machine Models"]
    SMM --> SMFULL["Simulink SM Full Controls<br/>81 blocks"]
    SMFULL --> SMPU["Synchronous Machine pu Fundamental<br/>(MathWorks library block)"]

    classDef verified fill:#0ca30c,color:#fff,stroke:#087a08;
    class GFMQV verified
```

The GFM subsystem (244 blocks) is the largest single unit — larger than
GFL (146) and the synchronous-machine block (81) combined — because it
implements *five* interchangeable outer-loop control types in one diagram
(Droop, VSM, dVOC, Matching, plus shared inner current/voltage loops and
virtual-impedance/PLL infrastructure), matching `symGFM_types.m`'s
multi-branch structure. This codebase currently ports only the `'Droop'`
branch (`components.gfm`) — the green node above.

## Block-type histogram (whole library, 1135 blocks / 55 subsystems)

```{mermaid}
%%{init: {"theme": "base"} }%%
xychart-beta
    title "Block types by count (top 15 of 32)"
    x-axis ["From", "Gain", "Goto", "Reference", "Sum", "Inport", "Constant", "Product", "SubSystem", "Outport", "PMIOPort", "Demux", "Terminator", "Mux", "Integrator"]
    y-axis "Count" 0 --> 170
    bar [159, 106, 101, 98, 84, 80, 61, 59, 54, 51, 48, 35, 21, 20, 18]
```

Only `PMComponent` (13) + `Ground` (6) + `PMIOPort` (48) — 67 of 1135
blocks, ~6% — are Simscape/SimPowerSystems *physical* (acausal circuit)
blocks, concentrated in the Transformer/Load/Upstream-Network subsystems
whose underlying physics (RL/RLC branches) is standard and already
implemented in this codebase's (P4) EMT component models. Everything else is
ordinary Simulink signal-flow with literal parameter values in the XML —
readable, just voluminous.

```{include} ../emt_inventory.md
:heading-offset: 1
```
