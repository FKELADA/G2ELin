# `reduction` — model order, and what a time-domain run means

G2ELin was built to model every element in full: the network's
electromagnetic transients, a machine's stator flux, a converter's filter
and inner loops. That is the right model for studying converter-network
interaction, and the wrong one for a 118-bus electromechanical study — not
because it is inaccurate, but because those fast dynamics are what makes a
model stiff, and a stiff model of a large network is a model nobody waits
for.

Model order reduction is the switch between the two. It is the same network,
the same operating point and the same equations; what changes is which
states are still integrated.

## The one idea

Every element's model is built as a DAE: differential equations for its
states, algebraic constraints, outputs. Reducing the order means **moving an
equation from the first group to the second**:

$$\frac{dx}{dt} = f(x, z, u) \qquad\longrightarrow\qquad 0 = f(x, z, u)$$

Nothing else changes. The linearisation
({meth}`~g2elin_core.components.base.ComponentDAE.linearize`) already
eliminates algebraic variables, the interconnection's closed-loop formula
already resolves the algebraic loops this creates — with a quasi-stationary
network that elimination *is* the admittance-matrix solve a phasor tool does
— and the nonlinear solver already solves `g(x, z, u) = 0` coupled across
the whole network. Even the initialisation is untouched: the operating point
is looked up by symbol, so a state that moves from `x` to `z` simply arrives
in a different vector.

## Two reductions, not one

```{warning}
`algebraic` and `frozen` are **not** two words for "drop this state". Using
one where the other belongs is the most common way a hand-rolled reduced
model goes quietly wrong.
```

| Mode | Operation | Valid for | Means |
|---|---|---|---|
| `dynamic` | `dx/dt = f` | — | the state is integrated |
| `algebraic` | `0 = f` | states **faster** than the study | the state responds instantly |
| `frozen` | `x ≡ x₀` | states **slower** than the study | the state does not move |

`algebraic` is residualization, the singular-perturbation limit: correct for
stator flux, filter currents, network branches, inner control loops.
`frozen` is truncation: correct for the field flux in the classical machine
model, where `E'` is held constant.

The field winding is the case that catches people. Its time constant
(`T'do`) is 5–10 s — a *slow* state. Setting `dψ_fd/dt = 0` gives an
*instantaneous* field winding, the exact opposite of the classical model.
That is why the 2nd-order level freezes it and every other level keeps it
dynamic.

## The catalogue

The unit of control is a **state group** — one physical approximation, such
as "stator flux" or "inner current loop" — not an individual state. A
**level** is a named preset over those groups. Anything a level expresses,
explicit group settings can express too; levels exist so the familiar models
have their familiar names.

### Network

| Level | Meaning |
|---|---|
| `full` | every passive element integrates its own `L di/dt` and `C dv/dt` |
| `quasi_stationary` | every passive element is algebraic — the phasor network |

Groups: `nodes`, `lines`, `loads`, `shunts`, `transformers`, so a mixed
network (dynamic lines, algebraic loads) is expressible.

`network_frequency` decides whether the `ω·L` and `ω·C` terms follow the
reference frame's own speed (the default, what an EMT model does) or are
pinned to nominal, which is what phasor tools do.

### The catalogue is per *unit*, not per type

A machine's AVR, stabiliser and governor, and a converter's power-control
law, are chosen per unit, and they decide which state groups that unit even
has: a machine with the Kundur exciter has an `avr` group of two states
rather than four, one with no governor has no `governor` group at all, and a
converter running VSM has `frequency` and `flux` where droop has
`power_filter`.

So `element(kind)` answers about the *type* — the default models — while
`Network.unit_element(der)` answers about one unit, which is what every
per-unit caller goes through. The **group ids are shared** wherever the
physics is, so a saved level or per-group override stays meaningful when a
model is swapped; only the states behind the group change.

### Synchronous machine

The controller states (governor, PSS, AVR) stay dynamic in every named
level: a reduced-order machine in a stability study keeps its full controls.
The totals below are for the default regulators — a different exciter or a
missing governor shifts every row by the same amount.

| Level | Machine states kept | Total |
|---|---|---|
| `full` | order-8 plus the step-up transformer current | 19 |
| `order8` | `ψ_d ψ_q ψ_fd ψ_1d ψ_1q ψ_2q Δω θ` | 17 |
| `order6` | stator flux algebraic (Sauer–Pai / Anderson–Fouad) | 15 |
| `order5` | 6th order without the second q-axis damper | 14 |
| `order4` | two-axis (`E'd`, `E'q`) | 13 |
| `order3` | flux decay / one-axis | 12 |
| `order2` | classical: swing only, field flux **frozen** | 11 |

### Grid-forming converter

| Level | Dropped | Total |
|---|---|---|
| `full` | — | 15 |
| `no_trafo` | transformer current | 13 |
| `no_filter` | + LC filter | 9 |
| `no_inner` | + inner current loop | 7 |
| `no_voltage` | + cascaded voltage loop | 5 |
| `droop` | + DC link → droop law and power filters only | 3 |

### Grid-following converter

| Level | Dropped | Total |
|---|---|---|
| `full` | — | 14 |
| `no_trafo` | transformer current | 12 |
| `no_filter` | + LCL filter | 8 |
| `no_inner` | + inner current loop | 6 |
| `no_dc` | + DC link | 4 |
| `pll` | + outer loops → a current source behind its PLL | 2 |

`θ`, `θ_pll` and the swing equation can never be removed: they define the
frames everything else is written in.

## What the choice makes the simulation

The web interface calls this the model *class*, and derives it rather than
asking for it:

| Class | When | A time-domain run is |
|---|---|---|
| **EMT** | network dynamic **and** units keep their fast electrical states | an electromagnetic-transient simulation |
| **RMS** | neither does | an electromechanical (phasor) simulation |
| **Mixed** | one but not the other | valid, and neither of the above |

This is why the interface's time-domain page is no longer called "EMT
simulation": what it produces depends on the model order, and both are
first-class.

## Pairing rules

A reduced model is not a pick-and-mix. Three combinations are buildable but
physically inconsistent, and {mod}`g2elin_core.network.validation` warns
about each:

1. **Quasi-stationary network + a unit that still integrates its transformer
   current.** The unit keeps an electromagnetic transient the grid it feeds
   no longer has.
2. **Quasi-stationary network + dynamic stator flux** (or a dynamic
   converter filter). The standard pairing drops both together; keeping only
   the machine's produces modes with no physical counterpart.
3. **Quasi-stationary network + `nodes_share_first_line_b`.** A
   quasi-stationary bus equation *is* the power flow's bus equation, so a
   shunt susceptance the power flow never saw puts the two out of agreement.
   Measured on WSCC-9: bus voltages settle 0.03 pu away from the power
   flow's and the units start ~3000× further from equilibrium. A dynamic
   network hides this as a transient that decays in microseconds; a
   quasi-stationary one cannot.

They are warnings, not errors. Asking for a mixed-timescale model on purpose
— to find out what one element's fast dynamics contribute — is a legitimate
thing to do, and none of these makes the model unsolvable.

## The same categories classify modes

Every state group also declares what *kind* of phenomenon it is, and
{mod}`g2elin_core.modal.classify` uses that to label each mode by which kind
of state does most of the participating in it. The eigenvalue map draws one
marker shape per kind.

| Kind | What is in it |
|---|---|
| **Synchronisation** | machine swing, grid-forming droop angle and power filters, grid-following PLL |
| **Control** | AVR, exciter, PSS, governor, converter voltage and current loops |
| **Unit electrical** | stator and rotor flux, damper windings, filters, DC link, step-up transformer current |
| **Network** | bus voltages and branch currents |
| **Mixed** | no kind holds half of it — often the interesting ones |
| **Reference angle** | the model's own free angles: coordinates, not dynamics |

The split is by *mechanism*, not by which box the state sits in. A converter
has no rotor, but its droop law and a machine's swing equation do the same
job, so they share a kind and "which are the synchronisation modes?" has one
answer across a mixed fleet.

```{note}
Synchronisation wins on a much smaller share than the others — 20% rather
than 50%. That is the field's own definition: a mode is electromechanical
when the rotor states participate *significantly*, not when they participate
most. Every machine's AVR and exciter take part in its local mode too, and
there are three times as many of those states, so a plain largest-share rule
hands rotor modes to the controllers riding on them. Measured on WSCC-9: the
1.63 Hz local mode is 32% rotor and 58% control.

The same rule correctly leaves the 0.66 and 0.76 Hz modes as control — those
are 1–2% rotor, field-flux and exciter modes that merely happen to sit in the
frequency band where an inter-area mode would be. Frequency alone would
misclassify them; participation does not.
```

## Is it safe? The adequacy check

The usual answer to "can I use a reduced model here?" is "it depends on your
case", which is no help in front of a specific case. {mod}`g2elin_core.modal.adequacy`
gives a specific answer, and it can, because the tool has both halves: the
full-order model, and participation factors.

It does two things:

**The screen.** Linearise at full order. For every state the reduction would
remove, measure its participation in the modes inside the band of interest.
A state that only ever participates in modes far outside that band can go —
that is the timescale separation singular perturbation needs, *verified on
this network* rather than assumed. A state carrying weight in a slow mode is
one whose removal will move that mode.

**The measurement.** Build the reduced model too and compare the spectra
mode by mode in the band, pairing them by proximity. So the report says both
"this looked risky" and "here is what it actually cost".

**A mode that disappears is not automatically a loss.** Dropping the damper
windings deletes their own time constants; that is the reduction doing what
it was asked. So each vanished mode is weighed by how much of it the removed
states accounted for: more than half, and it was theirs to take
(`expected_loss`); less, and the reduction has taken out something the
surviving states were part of. Without that distinction every real reduction
looked unsafe — on WSCC-9's classical model, 9 modes flagged as lost where
only 3 genuinely were — and a verdict that always says no tells you nothing.

Verdicts: `safe` (every mode reproduced), `check` (modes move — fine for
screening, not for controller design), `unsafe` (a mode the surviving states
were part of has no counterpart at all).

On WSCC-9 the check accepts a quasi-stationary network with 6th-order
machines and rejects the classical model, naming `ψ_2q` at 58 % participation
in a 4.7 Hz mode as the reason. That is a real property of that network, not
a rule of thumb.

## Where the settings live

On the {class}`~g2elin_core.network.schema.Network`, not on the analysis
request. The level is part of what a saved case *is*; it travels with an
exported network; and the API's model cache is keyed on the network's own
JSON, so changing a level invalidates every cached model for free.

```python
from g2elin_core.network.schema import ModelOptions

network.models = ModelOptions(
    network_level="quasi_stationary",
    sm_level="order6",
    gfm_level="droop",
    sm_states={"avr": "frozen"},      # per-group override, on top of the level
)
network.der_units[2].level = "order2"  # and per unit, on top of that
```

## Performance

Reduction removes states, but what decides how long a run takes is the
stiffness. On WSCC-9, a quasi-stationary network with 6th-order machines
goes from 88 states to 46 — and from a fastest mode at 8.2 × 10⁶ rad/s to
6.3 × 10³, a factor of 1300. That ratio is the number to watch, and it is
why the adequacy report prints it.

```{note}
Reduction moves states *into* the algebraic vector, so the coupled Newton
solve gets larger, not smaller. It is only a win because that solve is
sparse: the Jacobian of a 118-bus model is 0.055 % non-zero, and solving it
sparsely rather than densely took one right-hand-side evaluation there from
15.9 s to 0.037 s. See {mod}`g2elin_core.timedomain.emt`.
```

## Further reading

The framework is singular perturbation theory; see Kundur §5.3 and §13 for
the machine-side pairing rules, IEEE Std 1110 for the model-order naming,
and Milano *et al.*, *Foundations and Challenges of Low-Inertia Systems*
(PSCC 2018) for why the quasi-stationary assumption is under pressure in
converter-dominated grids. The adequacy check exists because that last
question — *when* is RMS not good enough — is usually answered only after
the fact, by an EMT run that disagrees.

```{eval-rst}
.. automodule:: g2elin_core.reduction
   :members:
   :undoc-members:

.. automodule:: g2elin_core.modal.adequacy
   :members:
   :undoc-members:
```
