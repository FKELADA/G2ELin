"""What *kind* of thing each mode is.

An eigenvalue map of a 100-state model is a cloud of dots. The question an
engineer actually brings to it -- "which of these are the electromechanical
modes, and which are control?" -- is answered by participation factors,
which the analysis already computes: a mode belongs to whatever kind of
state does most of the participating in it.

The categories come from :mod:`g2elin_core.reduction`, where every state
group already declares one. That is deliberate reuse rather than a second
taxonomy: the model-order page and the eigenvalue map then agree about what
a stator flux *is*, and adding a state group anywhere classifies its modes
everywhere without another table to keep in step.

A mode with no clear owner is not forced into a category. "Mixed" is a real
and often interesting answer -- a control loop interacting with the network
is exactly the kind of thing worth seeing marked as such, rather than
rounded to whichever side happens to hold 34%.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from g2elin_core import reduction
from g2elin_core.modal.analysis import ModalAnalysisResult, reference_angle_modes

# A category has to account for this much of a mode's participation to own
# it. Below that the mode is genuinely shared and is reported as mixed.
DOMINANCE = 0.5
# ...except synchronisation, which wins on a much smaller share. This is not
# a fudge, it is the field's own definition: a mode is electromechanical when
# the rotor angle and speed states participate *significantly*, not when they
# participate most. The AVR and exciter of every machine take part in its
# local mode too, and there are three times as many of those states, so a
# plain "largest share" rule hands a rotor mode to the controllers that ride
# on it. Measured on WSCC-9: the 1.63 Hz local mode is 32% rotor and 58%
# control.
#
# The same rule correctly leaves the 0.66 and 0.76 Hz modes as control: those
# are 1-2% rotor -- exciter and field-flux modes, not electromechanical ones.
SYNCHRONISATION_SHARE = 0.2


@dataclass(frozen=True)
class ModeCategory:
    """One mode's classification."""

    mode: int
    category: str            # a key of reduction.CATEGORIES
    share: float             # how much of the mode that category accounts for, 0..1
    shares: dict[str, float]  # every category's share, for a tooltip

    def as_dict(self) -> dict:
        return {"mode": self.mode, "category": self.category,
                "share": self.share, "shares": self.shares}


def state_categories(state_names: list[str]) -> list[str]:
    """Each state's category, in the model's own state order."""
    return [reduction.category_of_state(n) for n in state_names]


def classify_modes(result: ModalAnalysisResult) -> list[ModeCategory]:
    """Classify every mode of a solved modal analysis.

    The model's free reference angles are singled out first: they are
    coordinates rather than dynamics (see
    :func:`~g2elin_core.modal.analysis.reference_angle_modes`), so calling
    them "synchronisation modes" because they are made of angle states would
    be true and useless.
    """
    categories = state_categories(result.state_names)
    rows = {c: [] for c in reduction.CATEGORIES}
    for i, c in enumerate(categories):
        rows.setdefault(c, []).append(i)
    reference = set(reference_angle_modes(result))

    out: list[ModeCategory] = []
    for j in range(len(result.eigenvalues)):
        column = np.abs(result.participation[:, j])
        total = float(column.sum()) or 1.0
        shares = {c: float(column[idx].sum()) / total for c, idx in rows.items() if idx}
        if j in reference:
            out.append(ModeCategory(j, reduction.REFERENCE, 1.0, shares))
            continue
        synchronising = shares.get(reduction.SYNCHRONISATION, 0.0)
        if synchronising >= SYNCHRONISATION_SHARE:
            out.append(ModeCategory(j, reduction.SYNCHRONISATION, synchronising, shares))
            continue
        best = max(shares, key=lambda c: shares[c]) if shares else reduction.MIXED
        if not shares or shares[best] < DOMINANCE:
            out.append(ModeCategory(j, reduction.MIXED, shares.get(best, 0.0), shares))
        else:
            out.append(ModeCategory(j, best, shares[best], shares))
    return out
