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
