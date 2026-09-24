# Cross-validation against independent references

**September 2026.** What this report is for: everything G2ELin had been
checked against until now was itself. The model-order reduction added in this
cycle is validated against the tool's own full-order model, the nonlinear
model against the tool's own linearisation, the presets against the tool's
own power flow. All of that establishes *self-consistency* — that the pieces
agree with each other — and none of it establishes that the underlying
physics is right. A systematic error in the machine model would pass every
one of those tests.

This is the first check against references that were not produced by this
code.

> **A dated record.** The numbers below are what the presets gave in
> September 2026 and are left as they were. `kundur_two_area_classic` has
> since been rebuilt on the Kundur exciter model with the book's own values
> and an explicit "no governor", so a re-run will not reproduce them
> exactly. What the report concluded about *why* the models differ still
> holds; the figures are a snapshot, not a current baseline.

---

## Summary

| Reference | What it tests | Result |
|---|---|---|
| MATLAB toolbox's own state matrices | every component's linearisation | **7 of 7 match to ≤3×10⁻¹⁶** |
| ANDES 2.0, Kundur two-area | electromechanical modes, whole-network | **all three modes within 1.5 %** once the cases are made the same case |
| Kundur's book | published operating point | already tested; **eigenvalues not verified** (see Limitations) |

Two things this does **not** establish, stated up front because they matter
more than the numbers above: the grid-forming and grid-following converter
models have been checked only at component level, never against another tool
at network level; and the IEEE 14/39/118 presets cannot be used for dynamic
validation at all (they carry one invented machine, see Limitations).

---

## 1. Component linearisation, against the MATLAB toolbox

### What the reference is

`matlab/Symbolic/A_*.txt` holds the state matrices the original MATLAB
toolbox generated, entry by entry, as symbolic expressions. These are a
genuinely independent reference: different code, different language, and the
equations written out a second time by hand. They were sitting in the
repository unused.

The Python port and the MATLAB source use identical symbol names, identical
state ordering and identical parameter names, so the two matrices are
directly comparable without any translation.

### Method

`A = Fx − Fz·Gz⁻¹·Gx` is a *formula*. Both tools compute that formula; the
comparison evaluates both at the same randomly chosen point in parameter and
equilibrium space, five points per component. That tests the whole expression
rather than one operating point — a disagreement anywhere in it shows up. The
point does not need to be a physical equilibrium, only a place where `Gz`
inverts.

Reproduce with `python tools/compare_matlab_symbolic.py`.

### Result

| Component | Variant | States | Non-zero entries | Worst relative error |
|---|---|---|---|---|
| Synchronous machine | slack | 19 | 80 | 2.4 × 10⁻¹⁶ |
| Grid-forming converter (droop) | non-slack | 15 | 77 | 3.1 × 10⁻¹⁶ |
| Grid-following converter | non-slack | 14 | 59 | 1.9 × 10⁻¹⁶ |
| Infinite bus | slack | 3 | 4 | 0 |
| Line | — | 2 | 4 | 0 |
| Load | — | 2 | 4 | 0 |
| Node | — | 2 | 2 | 0 |

Machine precision throughout. Every component's linearisation reproduces the
original toolbox exactly.

### What this does and does not prove

It proves the **port is faithful** — that the migration from MATLAB to Python
introduced no error in any component's equations or in the algebraic
elimination. That is worth having: it is the foundation everything else in
the tool is built on, and it is now checked rather than assumed.

It does **not** prove the equations are right. Both tools implement the same
physics from the same source, so a modelling error in the original would be
reproduced faithfully and this check would still pass. For that, section 2.

### A detail worth recording

The reference files carry no record of which slack variant they were
generated for. The match itself identifies it: the machine file is the slack
variant, the converter files are non-slack — which is how `script_generic.m`
builds a case. A wrong variant differs in whole columns, not in the last
digit, so there is no ambiguity.

---

## 2. Electromechanical modes, against ANDES

### What the reference is

[ANDES](https://github.com/CURENT/andes) 2.0, an open-source power system
simulator with its own eigenvalue routine, its own machine models in standard
time-constant form (GENROU), its own exciter and governor models, and its own
Kundur two-area case. No shared ancestry with this code.

### Most of the work was making the two cases the same case

This is the part worth reading. Two tools disagreeing usually means they were
asked different questions. Every difference below was found by inspecting the
case data, not assumed:

**What matched.** ANDES's `kundur_full` carries the book's machine data
exactly — H = 6.5 / 6.175 s on the 900 MVA base, X_d = 1.8, X_q = 1.7,
X'_d = 0.3, X'_q = 0.55, X''_d = 0.25, and all four open-circuit time
constants (8.0, 0.03, 0.4, 0.05 s). Both tools use constant-impedance loads
for eigenvalue analysis.

**What did not.**

- **The network is different.** ANDES has three circuits on the 7–8 tie where
  the book has two, and doubles the 5–6, 6–7 and 9–10 lines. It has **no
  shunt capacitor banks**; the book specifies 200 MVAr at bus 7 and 350 MVAr
  at bus 9.
- **The controls are different.** ANDES drives its machines with EXDC2 and
  TGOV1. This tool's `kundur_two_area_classic` preset uses a fast exciter with
  no stabiliser, which is the whole point of Example 12.6.
- **`kundur_gencls.dyr` is not the book.** That case uses H = 13 s (double)
  and X'_d = 0.25 (against 0.3). It was not used.

### Result

Each row adds one matched element. Reproduce with `python tools/compare_andes.py`.

| Case | Inter-area | Local 1 | Local 2 |
|---|---|---|---|
| **ANDES `kundur_full`** | 0.6473 Hz / 3.43 % | 1.1120 Hz / 8.66 % | 1.1459 Hz / 8.86 % |
| G2ELin, preset as shipped | 0.6146 Hz / −2.67 % | 1.1260 Hz / 5.33 % | 1.1580 Hz / 5.45 % |
| | −5.04 % / −6.10 pp | +1.26 % / −3.32 pp | +1.05 % / −3.40 pp |
| + ANDES's exciter and governor | 0.5929 Hz / 8.40 % | 1.1246 Hz / 8.89 % | 1.1608 Hz / 8.72 % |
| | −8.40 % / +4.97 pp | +1.14 % / **+0.23 pp** | +1.30 % / **−0.14 pp** |
| **+ ANDES's third tie circuit** | **0.6376 Hz / 7.28 %** | **1.1270 Hz / 8.78 %** | **1.1626 Hz / 8.64 %** |
| | **−1.49 %** / +3.85 pp | **+1.35 %** / +0.12 pp | **+1.46 %** / −0.21 pp |

**All three electromechanical modes agree to within 1.5 % in frequency
between two independent implementations.** The local modes agree to within a
quarter of a percentage point in damping.

### Reading the intermediate steps

The steps are as informative as the endpoint, because each isolates one thing:

- **The exciter governs the damping.** Matching it moved the local modes from
  −3.4 pp to within ±0.25 pp. The tool's own preset gives the inter-area mode
  *negative* damping where ANDES's EXDC2 gives +3.4 % — which is not a
  disagreement but the textbook result: a fast exciter destabilises the
  inter-area mode, a slow DC exciter does not. That is what Example 12.6 is
  for.
- **The tie circuit governs the inter-area frequency.** Adding ANDES's third
  circuit moved it 0.5929 → 0.6376 Hz, closing an 8.4 % gap to 1.5 %. A
  stronger tie raises the frequency at which the two areas swing against each
  other; the local modes barely moved (+0.2 %), because they do not cross the
  tie. The selectivity is the evidence.
- **The shunt banks pull the other way.** Removing them drops the inter-area
  mode back to 0.5901 Hz. ANDES not having them is the reason the third
  circuit alone does not close the gap completely.

### What remains unexplained

The inter-area **damping** still differs by 3.9 pp (7.28 % against 3.43 %)
after everything above. The most likely cause is the governor: TGOV1 is a
lead-lag with T2 = 2.1 s and T3 = 7.0 s, while this tool's governor is
first-order, so only the droop and T1 = 0.49 s carry across. Those time
constants are comparable to the inter-area mode's 1.6 s period and far from
the local modes' 0.9 s — which fits the pattern exactly, since the local
modes' damping agrees to 0.2 pp.

This was **not** confirmed. An attempt to remove the governor from both tools
failed: ANDES's initialisation did not converge with TGOV1 out of service, so
that comparison was discarded rather than reported. Confirming it needs
either a lead-lag governor in this tool or a first-order one in ANDES.

---

## 3. A side result: the reduction validated independently

The comparison above was run on the **reduced** model — quasi-stationary
network, 6th-order machines — because that is the model class ANDES uses. So
the 1.5 % agreement is evidence for the model-order reduction work as well as
for the underlying models.

Separately, on the same network, moving from the full model to
quasi-stationary + 6th order shifts the modes by:

| | Full order | RMS (order6 + quasi-stationary) | Shift |
|---|---|---|---|
| Inter-area | 0.6140 Hz | 0.6146 Hz | +0.10 % |
| Local 1 | 1.1244 Hz | 1.1260 Hz | +0.14 % |
| Local 2 | 1.1561 Hz | 1.1580 Hz | +0.16 % |

The reduction is, on this network, worth about a tenth of a percent — well
inside the 1.5 % agreement with ANDES, so it is not what limits it.

---

## Limitations

These are the reasons not to over-read the results above.

**The IEEE presets cannot be used for dynamic validation.** IEEE 14, 39 and
118 were imported from pandapower, which carries no dynamic data. Every
machine in each case therefore shares **one** parameter set — a generic
machine, not the published IEEE dynamic data:

```
ieee14: 5 units, 1 distinct parameter set
ieee39: 10 units, 1 distinct parameter set
```

Their eigenvalues are a property of an invented machine. They are fine for
topology, power flow, and performance benchmarking — which is what they were
imported for — and meaningless for comparison against published modal
results. Kundur is currently the *only* preset carrying real, per-machine,
published dynamic data.

**Kundur's published eigenvalues were not verified.** The intention was to
check against the book's own tabulated mode frequencies and damping ratios.
The web search needed to confirm those published values hit a session rate
limit, and quoting them from memory in a validation report would defeat its
purpose. The existing test suite checks the book's *operating point* (the
400 MW tie flow, the rotor-angle differences) and checks the modes only
against loose ranges. **Tightening that against the book's printed table
remains the single highest-value item outstanding**, and it needs someone
with the book open.

**The converters have no network-level reference.** GFM and GFL match the
MATLAB toolbox at component level (section 1), and nothing beyond that. ANDES
has converter models (REGCA1, REGCP1 and similar) but they are
positive-sequence grid-following models of a different structure, not the
averaged dq-frame GFM/GFL implemented here, so a comparison would be
measuring the difference between two model *definitions*. For a converter
reference, the repository's own Simulink models are the better path — see
below.

**One network, one operating point.** Everything in section 2 is the Kundur
two-area system at one dispatch. Agreement there does not guarantee agreement
on a network with converters, a different loading, or islanded operation.

**The remaining inter-area damping difference is diagnosed, not proven.**

---

## What to do next

In order of value per unit of effort:

1. **Kundur's book eigenvalues.** Cheapest and most defensible. The preset,
   the data and the test file already exist; this is transcribing a printed
   table and tightening an existing test. Needs the book.

2. **The Simulink models already in this repository.** `matlab/CIGRE/` and
   `Library_V5_2024/` contain the EMT models these presets were ported from.
   They are the natural reference for the *full-order* model and for the
   converters — the two things section 2 could not reach. Needs MATLAB.

3. **A converter reference.** Either the Simulink models above, or building a
   matched case in a tool with averaged converter models. Worth doing: the
   converters are what this tool is *for*, and they are the least validated
   part of it.

4. **A second ANDES network.** IEEE-39 with real dynamic data attached would
   test whether the 1.5 % agreement holds on a larger system. This requires
   sourcing published IEEE-39 machine data and entering it — which would also
   fix the limitation above, and is the only way to make the IEEE presets
   useful for anything dynamic.

---

## Reproducing

```bash
cd python
python tools/compare_matlab_symbolic.py   # section 1, no extra dependencies
pip install andes
python tools/compare_andes.py             # section 2
```

ANDES generates its model code on first use. In a sandbox that blocks
multiprocessing, prepare it single-threaded first:

```python
import andes; andes.prepare(quick=True, ncpu=1)
```
