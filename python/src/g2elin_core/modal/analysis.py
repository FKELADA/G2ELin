"""Eigenvalue / participation-factor modal analysis, ported from the numeric
core of ``Functions/modal_analysis.m`` (plotting and the sensitivity tensor
aren't ported yet — this covers the eigenvalues, damping/frequency table,
and the participation-factor matrix).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
import scipy.linalg


@dataclass
class ModalAnalysisResult:
    eigenvalues: np.ndarray  # sorted ascending by real part, matching modal_analysis.m
    right_eigenvectors: np.ndarray  # columns
    left_eigenvectors: np.ndarray  # rows (already normalized so left @ right = I)
    participation: np.ndarray  # [state, mode], each column sums to 1
    state_names: list[str]

    def summary_table(self) -> pd.DataFrame:
        """One row per mode: eigenvalue, frequency/damping, and the top
        3 participating states — mirrors the table ``modal_analysis.m`` prints.
        """
        rows = []
        for i, lam in enumerate(self.eigenvalues):
            undamped_hz = abs(lam) / (2 * np.pi)
            damp_ratio = -lam.real / abs(lam) if abs(lam) > 0 else 0.0
            damped_hz = undamped_hz * np.sqrt(max(0.0, 1 - damp_ratio**2))
            part_col = np.abs(self.participation[:, i])
            top_idx = np.argsort(part_col)[::-1][:3]
            rows.append(
                {
                    "mode": i,
                    "real": lam.real,
                    "imag": lam.imag,
                    "undamped_hz": undamped_hz,
                    "damped_hz": damped_hz,
                    "damping_pct": damp_ratio * 100,
                    "state1": self.state_names[top_idx[0]],
                    "part1_pct": part_col[top_idx[0]] * 100,
                    "state2": self.state_names[top_idx[1]] if len(top_idx) > 1 else "",
                    "part2_pct": part_col[top_idx[1]] * 100 if len(top_idx) > 1 else 0.0,
                    "state3": self.state_names[top_idx[2]] if len(top_idx) > 2 else "",
                    "part3_pct": part_col[top_idx[2]] * 100 if len(top_idx) > 2 else 0.0,
                }
            )
        return pd.DataFrame(rows)


# A mode this close to the origin, made only of angle states, is one of the
# model's free reference angles rather than a physical mode -- see
# reference_angle_modes(). The physical modes of every case in this
# repository sit at least four orders of magnitude further out.
REFERENCE_MODE_TOL = 1e-4
_ANGLE_PARTICIPATION = 0.9


def reference_angle_modes(result: "ModalAnalysisResult") -> list[int]:
    """The modes that are only the model's own reference angles.

    Nothing pins the absolute position of the dq frame, so the model always
    has one marginal direction: turn every angle by the same amount and
    nothing physical changes. With a frame of its own
    (``components/frame.py``) there is a second one, since the frame and the
    unit it follows turn at the same speed and their difference is therefore
    conserved -- a redundancy of the coordinates, not a mode of the system.

    Both show up as eigenvalues at the origin (numerically a hair either
    side of it, which is why a plain ``max(Re) < 0`` test on the raw
    spectrum is the wrong question to ask) and both are made *entirely* of
    angle states, which is how they are told apart from a genuinely slow
    mode sitting near the origin.
    """
    angle = [i for i, n in enumerate(result.state_names) if n.startswith("theta")]
    out = []
    for j, lam in enumerate(result.eigenvalues):
        if abs(lam) <= REFERENCE_MODE_TOL and sum(result.participation[i, j] for i in angle) > _ANGLE_PARTICIPATION:
            out.append(j)
    return out


def analyze(A: np.ndarray, state_names: list[str]) -> ModalAnalysisResult:
    """Eigen-decompose ``A`` and compute the participation-factor matrix.

    Right eigenvectors ``V`` and left eigenvectors ``W`` are computed by a
    single ``scipy.linalg.eig(A, left=True, right=True)`` call so they come
    back consistently paired and ordered — computing them from two separate
    ``eig()`` calls (once on ``A``, once on ``A.conj().T``) and assuming a
    shared ordering is a real bug: nothing guarantees two independent LAPACK
    calls return eigenvalues in matching order, and mismatched pairing shows
    up as near-zero (or exactly zero) ``W @ V`` normalization coefficients.
    They're normalized so ``W @ V = I`` (``modal_analysis.m``'s
    ``coef = diag(W'*V)`` normalization), then the participation factor of
    state ``i`` in mode ``j`` is ``|V[i,j]| * |W[j,i]| / sum_k(|V[k,j]|*|W[j,k]|)``.
    """
    eigvals, W_raw, V = scipy.linalg.eig(A, left=True, right=True)
    order = np.argsort(eigvals.real)
    eigvals = eigvals[order]
    V = V[:, order]
    W_raw = W_raw[:, order]
    W = W_raw.conj().T  # rows are left eigenvectors, unnormalized

    coef = np.diag(W @ V)
    W = W / coef[:, None]

    n = A.shape[0]
    participation = np.zeros((n, n))
    for j in range(n):
        denom = np.sum(np.abs(W[j, :]) * np.abs(V[:, j]))
        for i in range(n):
            participation[i, j] = np.abs(V[i, j]) * np.abs(W[j, i]) / denom if denom else 0.0

    return ModalAnalysisResult(
        eigenvalues=eigvals,
        right_eigenvectors=V,
        left_eigenvectors=W,
        participation=participation,
        state_names=state_names,
    )
