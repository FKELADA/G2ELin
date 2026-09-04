"""The rest of ``Functions/modal_analysis.m``'s toolbox: eigenvalue
sensitivity, mode shapes, closed-form free-motion response, and MIMO step
response. ``modal/analysis.py``'s own docstring already flagged these as
"not ported yet" -- this is that follow-up, read directly from the MATLAB
source (not guessed) to keep the math and the plotting choices faithful.

Everything here is pure linear algebra on the eigendecomposition
:func:`~g2elin_core.modal.analysis.analyze` already computes, or (for step
response) a thin :mod:`scipy.signal` wrapper over the same ``A``/``B``/
``C``/``D`` :func:`~g2elin_core.pipeline.linearize_network` already
produces -- no new model type, no new derivation.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import scipy.signal

from g2elin_core.interconnect import AssembledSystem
from g2elin_core.modal.analysis import ModalAnalysisResult


@dataclass(frozen=True)
class SensitivityEntry:
    row_state: str
    col_state: str
    value: float


@dataclass(frozen=True)
class SensitivityResult:
    mode: int
    matrix: np.ndarray  # |sensitivity|, masked to A's structurally-nonzero entries, shape (n, n)
    top: list[SensitivityEntry]  # top-8 entries, matching modal_analysis.m's max_sens table column


def eigenvalue_sensitivity(modal: ModalAnalysisResult, A: np.ndarray, mode: int, n_top: int = 8) -> SensitivityResult:
    """``sensMatrix(k,j,i) = leftV(i,k)*rightV(j,i) / (leftV(i,:)*rightV(:,i))``
    for the chosen mode ``i`` (``modal_analysis.m`` Part 3) -- the standard
    eigenvalue-to-matrix-element sensitivity (Kundur eq. 12.75-ish):
    ``d(lambda_i)/d(A[k,j]) = w_i[k] * v_i[j]``, ``w_i``/``v_i`` the left/
    right eigenvectors of mode ``i``. Masked to where ``A`` is actually
    nonzero (``modal_analysis.m``'s ``Anum_ones``) -- a sensitivity value at
    a structurally-zero ``A`` entry isn't a "sensitivity to a real
    parameter," it's noise from an entry that was never a parameter to begin
    with.
    """
    w_i = modal.left_eigenvectors[mode, :]
    v_i = modal.right_eigenvectors[:, mode]
    denom = w_i @ v_i  # ~1 by construction (W@V=I), kept explicit rather than assumed
    sens = np.abs(np.outer(w_i, v_i) / denom)
    mask = A != 0
    sens = sens * mask

    flat_idx = np.argsort(sens, axis=None)[::-1][:n_top]
    top = []
    for idx in flat_idx:
        k, j = np.unravel_index(idx, sens.shape)
        if sens[k, j] <= 0:
            continue
        top.append(SensitivityEntry(row_state=modal.state_names[k], col_state=modal.state_names[j], value=float(sens[k, j])))
    return SensitivityResult(mode=mode, matrix=sens, top=top)


@dataclass(frozen=True)
class ModeShapeResult:
    mode: int
    states: list[str]
    angles_deg: list[float]  # magnitude is always 1 by construction (modal_analysis.m's ModeMatrix)


def mode_shape(modal: ModalAnalysisResult, mode: int, n_top: int = 5) -> ModeShapeResult:
    """Phase angle of the top-``n_top`` participating states' right
    eigenvector components, for polar plotting (``modal_analysis.m`` Part:
    "Plotting the mode shape"). Magnitude is fixed at 1 in the original
    (``ModeMatrix(k,1:2,i) = [1 angle(rightV(k,i))]``) -- it's a *shape*
    plot (relative phase between states' oscillation), not a magnitude one.
    """
    part_col = np.abs(modal.participation[:, mode])
    top_idx = np.argsort(part_col)[::-1][:n_top]
    v = modal.right_eigenvectors[:, mode]
    return ModeShapeResult(
        mode=mode,
        states=[modal.state_names[i] for i in top_idx],
        angles_deg=[math.degrees(float(np.angle(v[i]))) for i in top_idx],
    )


def free_response(modal: ModalAnalysisResult, perturb_index: int, offset: float, t: np.ndarray) -> np.ndarray:
    """Closed-form modal-expansion response to an initial-condition
    perturbation, no ODE integration needed (``modal_analysis.m``'s "Free
    Motion" section): ``x(t) = V @ (c * exp(eig*t))``, ``c = W @ x0`` the
    modal coordinates of the perturbation ``x0`` (a unit vector at
    ``perturb_index`` scaled by ``offset``). Returns ``x(t)`` for *every*
    state (shape ``(n_states, len(t))``) -- cheap (this is an O(n^2) matrix
    product per timestep for these state counts, not a solve), so there's
    no reason to compute only one state's row and force a second call for
    another.
    """
    n = len(modal.state_names)
    x0 = np.zeros(n, dtype=complex)
    x0[perturb_index] = offset
    c = modal.left_eigenvectors @ x0  # modal coordinates, shape (n,)
    exp_terms = np.exp(np.outer(modal.eigenvalues, t))  # (n_modes, n_t)
    x_t = modal.right_eigenvectors @ (c[:, None] * exp_terms)  # (n_states, n_t)
    return np.real(x_t)  # imaginary parts cancel across conjugate-pair modes up to float noise


def step_response(
    system: AssembledSystem, input_name: str, output_name: str, amplitude: float, t: np.ndarray
) -> np.ndarray:
    """MIMO step response between one chosen input/output pair
    (``modal_analysis.m``'s "step_resp" section), via a SISO reduction of
    the already-linearized ``(A, B, C, D)`` -- pick the one column of ``B``
    and row of ``C``/``D`` the chosen input/output correspond to, scale the
    unit-step response by ``amplitude`` afterward (linear system, so this
    is exact, not an approximation).
    """
    i = system.input_names.index(input_name)
    j = system.output_names.index(output_name)
    b_col = system.B[:, i : i + 1]
    c_row = system.C[j : j + 1, :]
    d_val = system.D[j : j + 1, i : i + 1]
    sys = scipy.signal.StateSpace(system.A, b_col, c_row, d_val)
    _, y = scipy.signal.step(sys, T=t)
    return amplitude * y
