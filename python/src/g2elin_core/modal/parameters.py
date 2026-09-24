"""Which *physical parameter* drives a mode.

The eigenvalue sensitivity ``dlambda/dA_ij`` says which entries of the state
matrix a mode is sensitive to. That is one step short of the question
actually being asked, because ``A_ij`` is not a thing anyone can change --
it is an expression in the physical parameters. What a user wants to know is
which parameter to reach for.

The two are connected by the chain rule, and the connection is exact:

    dlambda/dp = sum_ij (dlambda/dA_ij) (dA_ij/dp)

The first factor is the sensitivity the modal toolbox already computes (as a
complex number here, not its magnitude -- the phase is what separates a
parameter that moves a mode left from one that moves it up). The second comes
from perturbing the parameter.

**Perturbing it cheaply is the whole difficulty.** Done the obvious way --
change the parameter, relinearise the network, subtract -- a 118-bus scan
takes nearly two hours: 1674 parameters, twice each, at 0.4 s a
linearisation. Two things make it four seconds instead, and both rest on the
same observation: *a unit's parameter only touches that unit*.

Its operating point is the unit's own, and everything an operating point is
built from besides the parameters is fixed by the power flow, so only that
unit is rebuilt rather than all 54
(:class:`~g2elin_core.components.base.RebuiltWithParams`).

And its four matrices are one block of the assembled system, so the assembly
-- which inverts an ``n_y``-square matrix, 0.34 s of that 0.4 s -- does not
have to be redone at all. The change is carried through it in closed form
instead (:class:`~g2elin_core.interconnect.assemble.PerturbationProjector`),
from weights built once before the scan.

What is left per parameter is one component linearisation and four small
inner products, well under a millisecond.

There is a single exception, and it is the reason the slack machine gets
special treatment below: in the MATLAB-compatible frame that machine's rotor
angle is the reference every other unit's angle is measured from, so the
three constants that set it really do reach outside its own block.

**It also answers the question symbolically, without any symbolic algebra.**
A parameter appears in ``A_ij`` exactly when ``dA_ij/dp`` is non-zero, and
that derivative is already being computed. So "which parameters live in this
entry" comes out of the same pass -- and it comes out for reduced-order
models too, where no printed symbolic matrix exists.

One consequence of working on the assembled model: a parameter is reported
as the *component* sees it, on the network's base, with a unit's own rating
(``DerUnit.sn_mva``) already rebased in. Sensitivities are given per 1% of
the value, which means the same thing on either base, so the choice does not
reach the answer -- but reading a value on one base and perturbing it on the
other would, and that is a mistake this module made once.

The worked example, on any machine: the rotor-speed row of ``A`` carries a
factor ``1/(2H)`` in every one of its entries, because the swing equation
divides by ``2H``. So an electromechanical mode's largest sensitivities land
in that row, and this module reports ``H`` as the parameter behind them.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from g2elin_core.components.gfl import linearize_gfl
from g2elin_core.components.gfm import linearize_gfm
from g2elin_core.components.sm import linearize_sm
from g2elin_core.interconnect import PerturbationProjector, assembly_parts
from g2elin_core.interconnect.network_assembly import build_blocks_and_wiring
from g2elin_core.modal.analysis import ModalAnalysisResult
from g2elin_core.network.schema import Network
from g2elin_core.operating_point import (
    compute_operating_point,
    overridable_param_keys,
    rebase_params,
)
from g2elin_core.pipeline import linear_components, linearize_network
from g2elin_core.powerflow import PowerFlowResult

# Relative step for the parameter perturbation. Large enough to stay clear of
# round-off in A, small enough that the second-order term is negligible; a
# central difference makes the error O(h^2).
RELATIVE_STEP = 1e-5
# Below this share of the largest |dA/dp| in an entry, a parameter is taken
# not to appear in it. Guards against a derivative that is really round-off.
PRESENCE_FLOOR = 1e-8


@dataclass(frozen=True)
class ParameterEffect:
    """What one parameter does to one mode."""

    unit: int
    unit_label: str
    parameter: str
    value: float
    d_lambda: complex          # dlambda/dp, in the parameter's own units
    # A 1% increase in the parameter, expressed the way the mode is read.
    d_freq_hz: float
    d_damping_pct: float

    @property
    def magnitude(self) -> float:
        """How much a 1% change moves the mode, as one number for ranking:
        the eigenvalue's displacement relative to its own size."""
        return float(abs(self.d_lambda) * abs(self.value) / 100.0)

    def as_dict(self) -> dict:
        return {
            "unit": self.unit, "unit_label": self.unit_label, "parameter": self.parameter,
            "value": self.value, "d_lambda_real": self.d_lambda.real,
            "d_lambda_imag": self.d_lambda.imag, "d_freq_hz": self.d_freq_hz,
            "d_damping_pct": self.d_damping_pct, "magnitude": self.magnitude,
        }


@dataclass(frozen=True)
class EntryParameters:
    """One high-sensitivity entry of ``A``, and the parameters inside it.

    This is the bridge the whole module is for: ``d(dw_r)/d(psi_d)`` means
    nothing on its own, but knowing that the entry is built out of ``H`` and
    the machine's inductances says what to change.
    """

    row_state: str
    col_state: str
    sensitivity: float
    parameters: list[str]       # "SM_1.H", strongest first

    def as_dict(self) -> dict:
        return {"row_state": self.row_state, "col_state": self.col_state,
                "sensitivity": self.sensitivity, "parameters": self.parameters}


@dataclass
class ParameterSensitivityResult:
    mode: int
    eigenvalue: complex
    effects: list[ParameterEffect]       # every parameter, strongest first
    entries: list[EntryParameters]       # the top |dlambda/dA| entries
    notes: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        return {
            "mode": self.mode,
            "eigenvalue_real": self.eigenvalue.real, "eigenvalue_imag": self.eigenvalue.imag,
            "effects": [e.as_dict() for e in self.effects],
            "entries": [e.as_dict() for e in self.entries],
            "notes": self.notes,
        }


def _block_index(blocks, unit_id: int, network: Network) -> int | None:
    """Where one DER unit's block sits in the assembled ordering."""
    from g2elin_core.network.breakers import block_labels

    name = block_labels(network).der.get(unit_id)
    if name is None:
        return None
    for i, b in enumerate(blocks):
        if b.name == name:
            return i
    return None


def parameter_sensitivity(
    network: Network,
    result: PowerFlowResult,
    modal: ModalAnalysisResult,
    mode: int,
    *,
    units: list[int] | None = None,
    parameters: list[str] | None = None,
    n_entries: int = 8,
    participation_floor: float = 0.0,
) -> ParameterSensitivityResult:
    """Rank unit parameters by how much each moves one mode.

    ``units`` and ``parameters`` narrow the scan explicitly.
    ``participation_floor`` narrows it by relevance: a unit none of whose
    states take part in the mode cannot move it much, and on a large network
    most units are in that position. Callers that want everything pass 0.

    The power flow is held fixed, so these are derivatives at a *fixed
    operating point*. For a control parameter that is the whole story. For an
    impedance it is not: changing ``Lt`` would also move the power flow, and
    only the dynamic half of that effect is reported here.
    """
    lam = complex(modal.eigenvalues[mode])
    w = modal.left_eigenvectors[mode, :]
    v = modal.right_eigenvectors[:, mode]
    # The complex sensitivity, not its magnitude: dlambda/dA_ij = w_i[k] v_i[j].
    sensitivity = np.outer(w, v) / (w @ v)

    op = compute_operating_point(network, result)
    ops = {**op.sm_ops, **op.gfm_ops, **op.gfl_ops}
    components = linear_components(network, result, op=op)
    blocks, wiring = build_blocks_and_wiring(network, **components)
    parts = assembly_parts(blocks, wiring)
    if sensitivity.shape != parts.A_tot.shape:
        raise ValueError(
            f"this modal result has {sensitivity.shape[0]} states but the network linearises to "
            f"{parts.A_tot.shape[0]}; they are not the same model"
        )
    projector = PerturbationProjector(parts, sensitivity)

    # The entries the mode is most sensitive to, decided before the scan: a
    # parameter's presence is then only ever evaluated at those few positions
    # instead of over the whole matrix.
    magnitude = np.abs(sensitivity) * (parts.A_tot != 0)
    top_cells: list[tuple[int, int]] = []
    for flat in np.argsort(magnitude, axis=None)[::-1][:n_entries]:
        i, j = np.unravel_index(int(flat), magnitude.shape)
        if magnitude[i, j] > 0:
            top_cells.append((int(i), int(j)))
    reach = _EntryReach(parts, top_cells)

    wanted_units = set(units) if units is not None else None
    effects: list[ParameterEffect] = []
    presence: dict[str, np.ndarray] = {}
    notes: list[str] = []
    skipped_units = 0

    for der in network.der_units:
        kind = der.unit_type.value
        if wanted_units is not None and der.id not in wanted_units:
            continue
        try:
            keys = sorted(overridable_param_keys(kind))
        except Exception:
            continue                      # an infinite bus has none
        if not keys:
            continue
        index = _block_index(blocks, der.id, network)
        if index is None:
            continue
        xs, _, _ = parts.block_slices(index)
        share = float(np.abs(modal.participation[xs, mode]).sum())
        if participation_floor and share < participation_floor:
            skipped_units += 1
            continue

        unit_op = ops.get(der.id)
        linearize = _LINEARIZERS.get(kind)
        if unit_op is None or linearize is None:
            continue
        # The parameters as the component actually sees them: on the
        # *network's* base, the unit's own overrides already rebased in. A
        # sensitivity is reported per 1% of the value, and a relative change
        # means the same thing on either base, so which one this is does not
        # reach the answer -- but reading p0 from one base and perturbing on
        # the other would, which is why the scan stays on this one.
        current = unit_op.p
        modes = network.unit_modes(der)
        label = f"{kind.upper()}_{der.id}"
        # The one unit whose parameters reach past its own block. In the
        # MATLAB-compatible frame the slack machine's rotor angle *is* the
        # reference every other unit's angle is measured from, so the three
        # constants that set it re-reference the whole network and have to go
        # the slow way round. With a frame component of its own (the default)
        # the reference angle is gauge and nothing of the sort happens.
        carries_frame = (
            network.frame_follows_slack and kind == "sm" and der.bus_type.value == "slack"
        )
        for name in keys:
            if parameters is not None and name not in parameters:
                continue
            p0 = float(current.get(name, 0.0))
            if p0 == 0.0:
                continue                  # no scale to perturb relative to
            step = abs(p0) * RELATIVE_STEP
            try:
                if carries_frame and name in ANGLE_SETTING_PARAMS:
                    dA_tot = _relinearised_dA(network, result, der, name, p0, step)
                    d_lambda = complex(np.sum(sensitivity * dA_tot))
                    reached = np.array([abs(dA_tot[i, j]) for i, j in top_cells])
                else:
                    up = linearize(unit_op.with_params({**current, name: p0 + step}), modes=modes)
                    down = linearize(unit_op.with_params({**current, name: p0 - step}), modes=modes)
                    dA = (up.A - down.A) / (2 * step)
                    dB = (up.B - down.B) / (2 * step)
                    dC = (up.C - down.C) / (2 * step)
                    dD = (up.D - down.D) / (2 * step)
                    if not all(np.all(np.isfinite(m)) for m in (dA, dB, dC, dD)):
                        continue
                    if not any(np.any(m) for m in (dA, dB, dC, dD)):
                        continue
                    d_lambda = projector.project(index, dA, dB, dC, dD)
                    reached = reach.at(index, dA, dB, dC, dD)
            except Exception as exc:
                notes.append(f"{label}.{name} skipped: {type(exc).__name__}: {exc}")
                continue
            presence[f"{label}.{name}"] = reached

            delta = d_lambda * (p0 / 100.0)
            effects.append(ParameterEffect(
                unit=der.id, unit_label=label, parameter=name, value=p0, d_lambda=d_lambda,
                d_freq_hz=_d_frequency(lam, delta), d_damping_pct=_d_damping(lam, delta),
            ))

    effects.sort(key=lambda e: -e.magnitude)
    if skipped_units:
        notes.append(
            f"{skipped_units} unit(s) skipped: none of their states take part in this mode "
            f"(participation below {participation_floor:g}). Lower the floor to include them."
        )

    entries: list[EntryParameters] = []
    for cell, (i, j) in enumerate(top_cells):
        present = {k: float(m[cell]) for k, m in presence.items() if m[cell] > 0}
        if present:
            top = max(present.values())
            present = {k: val for k, val in present.items() if val >= top * PRESENCE_FLOOR}
        entries.append(EntryParameters(
            row_state=modal.state_names[i], col_state=modal.state_names[j],
            sensitivity=float(magnitude[i, j]),
            parameters=[k for k, _ in sorted(present.items(), key=lambda kv: -kv[1])],
        ))

    return ParameterSensitivityResult(
        mode=mode, eigenvalue=lam, effects=effects, entries=entries, notes=notes,
    )


_LINEARIZERS = {"sm": linearize_sm, "gfm": linearize_gfm, "gfl": linearize_gfl}

#: A synchronous machine's initial rotor angle is
#: ``phase(v + (Ra + j(Ll + Laq)) i)``, so these three constants -- and only
#: these -- move it. They matter because in the MATLAB-compatible frame that
#: angle, for the slack machine, is the network's reference.
ANGLE_SETTING_PARAMS = frozenset({"Ra", "Ll", "Laq"})


def _relinearised_dA(network: Network, result: PowerFlowResult, der, name: str,
                     p0: float, step: float) -> np.ndarray:
    """``dA_tot/dp`` from relinearising the whole network, twice.

    The slow path, for the one case the projector cannot cover: a parameter
    whose effect is not confined to its own unit's block. Two full
    linearisations, so it is reserved for the three parameters that need it.
    """
    def A_at(value: float) -> np.ndarray:
        # An override is read per unit of the unit's own rating, while `value`
        # is on the network's -- see operating_point.rebase_params.
        raw = value if der.sn_mva is None else rebase_params(
            {name: value}, from_mva=network.sn_mva, to_mva=der.sn_mva
        )[name]
        units = [d.model_copy(update={"params": {**d.params, name: raw}})
                 if d.id == der.id else d for d in network.der_units]
        return linearize_network(network.model_copy(update={"der_units": units}), result).A

    return (A_at(p0 + step) - A_at(p0 - step)) / (2 * step)


class _EntryReach:
    """How far a change in one component reaches, at a few chosen entries.

    ``dA_tot = dA_c + dB_c P + Q dC_c + Q dD_c P`` is an ``n_x``-square
    matrix, and forming one per parameter is both slow and, on a 118-bus
    model, 27 MB a time. Only a handful of its entries are ever read -- the
    ones the mode is most sensitive to -- so only those are evaluated, which
    needs nothing bigger than the component's own blocks.
    """

    def __init__(self, parts, cells: list[tuple[int, int]]) -> None:
        self.parts = parts
        self.cells = cells
        topo = parts.topology
        GE = topo.G @ parts.E_ol
        self.P = GE @ parts.C_ol          # (n_u, n_x)
        self.Q = parts.B_ol @ GE          # (n_x, n_y)

    def at(self, index: int, dA, dB, dC, dD) -> np.ndarray:
        """``|dA_tot|`` at each chosen entry, for a change confined to
        block ``index``."""
        xs, us, ys = self.parts.block_slices(index)
        x0, x1 = xs.start, xs.stop
        out = np.zeros(len(self.cells))
        for k, (i, j) in enumerate(self.cells):
            value = 0.0
            in_rows, in_cols = x0 <= i < x1, x0 <= j < x1
            if dA.size and in_rows and in_cols:
                value += dA[i - x0, j - x0]
            if dB.size and in_rows:
                value += dB[i - x0, :] @ self.P[us, j]
            if dC.size and in_cols:
                value += self.Q[i, ys] @ dC[:, j - x0]
            if dD.size:
                value += self.Q[i, ys] @ dD @ self.P[us, j]
            out[k] = abs(value)
        return out


def _d_frequency(lam: complex, delta: complex) -> float:
    """Change in undamped natural frequency (Hz) for the given eigenvalue shift."""
    if abs(lam) == 0:
        return 0.0
    return float((abs(lam + delta) - abs(lam)) / (2 * np.pi))


def _d_damping(lam: complex, delta: complex) -> float:
    """Change in damping ratio (percentage points) for the given shift."""
    def zeta(x: complex) -> float:
        return -x.real / abs(x) * 100 if abs(x) > 0 else 0.0
    return float(zeta(lam + delta) - zeta(lam))
