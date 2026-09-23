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
parameter that moves a mode left from one that moves it up). The second is
obtained by perturbing one parameter and re-linearising, which is cheap now
that a linearisation is milliseconds rather than seconds.

**It also answers the question symbolically, without any symbolic algebra.**
A parameter appears in ``A_ij`` exactly when ``dA_ij/dp`` is non-zero, and
that derivative is already being computed. So "which parameters live in this
entry" comes out of the same pass -- and it comes out for reduced-order
models too, where no printed symbolic matrix exists.

The worked example, on any machine: the rotor-speed row of ``A`` carries a
factor ``1/(2H)`` in every one of its entries, because the swing equation
divides by ``2H``. So an electromechanical mode's largest sensitivities land
in that row, and this module reports ``H`` as the parameter behind them.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from g2elin_core.modal.analysis import ModalAnalysisResult
from g2elin_core.network.schema import Network
from g2elin_core.operating_point import overridable_param_keys, unit_params
from g2elin_core.pipeline import linearize_network
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


def _perturbed(network: Network, unit_id: int, name: str, value: float) -> Network:
    units = [
        d.model_copy(update={"params": {**d.params, name: value}}) if d.id == unit_id else d
        for d in network.der_units
    ]
    return network.model_copy(update={"der_units": units})


def parameter_sensitivity(
    network: Network,
    result: PowerFlowResult,
    modal: ModalAnalysisResult,
    mode: int,
    *,
    units: list[int] | None = None,
    parameters: list[str] | None = None,
    n_entries: int = 8,
) -> ParameterSensitivityResult:
    """Rank every unit parameter by how much it moves one mode.

    ``units`` and ``parameters`` narrow the scan; by default every parameter
    of every unit that has them is tried. One linearisation per parameter per
    direction, no eigendecomposition beyond the one already done -- tracking
    an eigenvalue across a perturbation is fragile, and the chain rule makes
    it unnecessary.

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

    base = linearize_network(network, result)
    n = base.A.shape[0]
    if sensitivity.shape != (n, n):
        raise ValueError(
            f"this modal result has {sensitivity.shape[0]} states but the network linearises to {n}; "
            "they are not the same model"
        )

    wanted_units = set(units) if units is not None else None
    effects: list[ParameterEffect] = []
    # |dA/dp| per entry, kept so an entry can say which parameters built it.
    presence: dict[str, np.ndarray] = {}
    notes: list[str] = []

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
        current = unit_params(network, der)
        label = f"{kind.upper()}_{der.id}"
        for name in keys:
            if parameters is not None and name not in parameters:
                continue
            p0 = float(current.get(name, 0.0))
            if p0 == 0.0:
                continue                  # no scale to perturb relative to
            step = abs(p0) * RELATIVE_STEP
            try:
                up = linearize_network(_perturbed(network, der.id, name, p0 + step), result).A
                down = linearize_network(_perturbed(network, der.id, name, p0 - step), result).A
            except Exception as exc:
                notes.append(f"{label}.{name} skipped: {type(exc).__name__}: {exc}")
                continue
            dA = (up - down) / (2 * step)
            if not np.all(np.isfinite(dA)) or not np.any(dA):
                continue
            d_lambda = complex(np.sum(sensitivity * dA))
            presence[f"{label}.{name}"] = np.abs(dA)

            # What a 1% increase does, read the way a mode is read.
            delta = d_lambda * (p0 / 100.0)
            effects.append(ParameterEffect(
                unit=der.id, unit_label=label, parameter=name, value=p0, d_lambda=d_lambda,
                d_freq_hz=_d_frequency(lam, delta), d_damping_pct=_d_damping(lam, delta),
            ))

    effects.sort(key=lambda e: -e.magnitude)

    # The entries the mode is most sensitive to, and what is inside them.
    magnitude = np.abs(sensitivity) * (base.A != 0)
    entries: list[EntryParameters] = []
    for flat in np.argsort(magnitude, axis=None)[::-1][:n_entries]:
        i, j = np.unravel_index(int(flat), magnitude.shape)
        if magnitude[i, j] <= 0:
            continue
        present = {k: float(m[i, j]) for k, m in presence.items() if m[i, j] > 0}
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
