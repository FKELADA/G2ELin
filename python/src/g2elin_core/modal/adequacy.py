"""Is this model-order reduction safe for *this* network?

The usual answer is "it depends on your case", which is useless in front of
a specific case. This module gives a specific answer, and it can, because
the tool already has both halves: the full-order model, and participation
factors.

**The method.** Linearize the network at full order once. For every state
the chosen reduction would remove, look at how much it participates in the
modes that survive -- specifically, in the slow modes the reduced model is
supposed to reproduce. A state that is only ever involved in modes far
outside the band of interest can be removed without consequence: that is
the timescale separation singular-perturbation theory needs, verified on
the actual case rather than assumed. A state carrying real participation in
a slow mode is a state whose removal will move that mode, and the reduction
is not safe.

**And then it checks.** The participation screen predicts; the second half
of the report measures. It builds the reduced model too and compares the
two spectra mode by mode in the band of interest, so the report says both
"this looked risky" and "here is what it actually cost you".

The band is a parameter because it is the study that defines it:
electromechanical oscillations live below ~2-3 Hz, converter-driven
interactions reach tens of Hz, and a sub-synchronous resonance study cares
about a band a classical stability study would happily discard.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from g2elin_core import reduction
from g2elin_core.modal.analysis import analyze, reference_angle_modes
from g2elin_core.network.schema import Network
from g2elin_core.pipeline import linearize_network
from g2elin_core.powerflow import PowerFlowResult

# Modes at or below this frequency are what a reduced model is normally
# asked to reproduce: the electromechanical band, with headroom for the
# slower converter control modes.
DEFAULT_BAND_HZ = 5.0
# Participation at or above this, in a mode inside the band, makes a state
# one the reduction should not be dropping. 5% is low on purpose -- this is
# a screen, and a false alarm costs the user a glance while a missed one
# costs them a wrong answer.
RISK_PARTICIPATION = 0.05
# An unmatched mode this much made of removed states is one the reduction was
# asked to remove, not one it lost by accident.
EXPECTED_LOSS_SHARE = 0.5


@dataclass
class StateRisk:
    """One removed state that carries weight in a mode worth keeping."""

    state: str
    group: str  # the reduction group it belongs to, e.g. "sm.stator_flux"
    mode_hz: float
    mode_damping_pct: float
    participation: float

    def as_dict(self) -> dict:
        return {
            "state": self.state,
            "group": self.group,
            "mode_hz": self.mode_hz,
            "mode_damping_pct": self.mode_damping_pct,
            "participation": self.participation,
        }


@dataclass
class ModePair:
    """One mode of the full model and its match in the reduced one."""

    full_hz: float
    full_damping_pct: float
    reduced_hz: float | None
    reduced_damping_pct: float | None
    # How much of this mode the removed states accounted for, 0..1. A mode
    # made entirely of states the reduction deletes is *supposed* to
    # disappear with them -- see :attr:`expected_loss`.
    removed_share: float = 0.0

    @property
    def matched(self) -> bool:
        return self.reduced_hz is not None

    @property
    def expected_loss(self) -> bool:
        """An unmatched mode that was the removed states' own.

        Dropping the damper windings deletes their time constants; that is
        the reduction working, not failing. Counting those as lost dynamics
        made every real reduction look unsafe, which would have made the
        verdict worthless."""
        return not self.matched and self.removed_share >= EXPECTED_LOSS_SHARE

    @property
    def d_hz(self) -> float:
        return abs(self.reduced_hz - self.full_hz) if self.matched else float("nan")

    @property
    def d_damping_pct(self) -> float:
        return abs(self.reduced_damping_pct - self.full_damping_pct) if self.matched else float("nan")

    def as_dict(self) -> dict:
        return {
            "full_hz": self.full_hz,
            "full_damping_pct": self.full_damping_pct,
            "reduced_hz": self.reduced_hz,
            "reduced_damping_pct": self.reduced_damping_pct,
            "matched": self.matched,
            "expected_loss": self.expected_loss,
            "removed_share": self.removed_share,
            "d_hz": None if not self.matched else self.d_hz,
            "d_damping_pct": None if not self.matched else self.d_damping_pct,
        }


@dataclass
class AdequacyReport:
    band_hz: float
    n_states_full: int
    n_states_reduced: int
    removed_states: list[str]
    risks: list[StateRisk]
    modes: list[ModePair]
    max_d_hz: float
    max_d_damping_pct: float
    unmatched: int
    stiffness_full: float
    stiffness_reduced: float
    verdict: str  # "safe" | "check" | "unsafe" | "full_order"
    notes: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        return {
            "band_hz": self.band_hz,
            "n_states_full": self.n_states_full,
            "n_states_reduced": self.n_states_reduced,
            "removed_states": self.removed_states,
            "risks": [r.as_dict() for r in self.risks],
            "modes": [m.as_dict() for m in self.modes],
            "max_d_hz": self.max_d_hz,
            "max_d_damping_pct": self.max_d_damping_pct,
            "unmatched": self.unmatched,
            "stiffness_full": self.stiffness_full,
            "stiffness_reduced": self.stiffness_reduced,
            "verdict": self.verdict,
            "notes": self.notes,
        }


def _full_order(network: Network) -> Network:
    """The same network with every dynamic kept -- the reference to judge a
    reduction against."""
    from g2elin_core.network.schema import ModelOptions

    units = [
        d.model_copy(update={"level": None, "states": {}}) if d.level or d.states else d
        for d in network.der_units
    ]
    return network.model_copy(update={"models": ModelOptions(), "der_units": units})


def _group_of(state_name: str, network: Network) -> str:
    """Which reduction group a full-model state name belongs to, as
    ``"<kind>.<group>"`` -- the catalogue's own lookup, which the eigenvalue
    map's mode classification shares."""
    return reduction.group_of_state(state_name) or "?"


def _band_mode_indices(result, band_hz: float) -> list[int]:
    """Indices of the modes at or below ``band_hz``, one per complex pair and
    skipping the model's own free reference angles (coordinates, not
    dynamics).

    Taking one of each conjugate pair matters beyond tidiness: a mode and its
    conjugate have identical participation factors, so counting both listed
    every finding twice.
    """
    skip = set(reference_angle_modes(result))
    return [
        i for i, lam in enumerate(result.eigenvalues)
        if i not in skip and lam.imag >= 0 and abs(lam) / (2 * np.pi) <= band_hz
    ]


def _modes_in_band(result, band_hz: float) -> list[tuple[float, float, int]]:
    """``(frequency Hz, damping %, mode index)`` of every mode in the band."""
    out = []
    for i in _band_mode_indices(result, band_hz):
        lam = result.eigenvalues[i]
        hz = abs(lam) / (2 * np.pi)
        damping = -lam.real / abs(lam) * 100 if abs(lam) > 0 else 100.0
        out.append((hz, damping, i))
    return sorted(out)


def _pair_up(
    full: list[tuple[float, float, int]],
    reduced: list[tuple[float, float, int]],
    removed_share: dict[int, float] | None = None,
) -> list[ModePair]:
    """Match each full-model mode to its nearest unused reduced-model mode.

    Nearest-in-frequency, greedily, closest pair first: a reduction moves
    modes a little, so identity is by proximity, and matching the most
    confident pairs first keeps one badly-shifted mode from stealing
    another's partner.
    """
    remaining = list(range(len(reduced)))
    pairs: dict[int, int] = {}
    candidates = sorted(
        ((abs(f[0] - r[0]), i, j) for i, f in enumerate(full) for j, r in enumerate(reduced)),
        key=lambda c: c[0],
    )
    used_full: set[int] = set()
    for _, i, j in candidates:
        if i in used_full or j not in remaining:
            continue
        # A mode that moved by more than half its own frequency is not the
        # same mode; leave it unmatched rather than claim a false pairing.
        if abs(full[i][0] - reduced[j][0]) > max(0.5 * full[i][0], 0.05):
            continue
        pairs[i] = j
        used_full.add(i)
        remaining.remove(j)
    share = removed_share or {}
    return [
        ModePair(
            full_hz=f[0], full_damping_pct=f[1],
            reduced_hz=reduced[pairs[i]][0] if i in pairs else None,
            reduced_damping_pct=reduced[pairs[i]][1] if i in pairs else None,
            removed_share=share.get(f[2], 0.0),
        )
        for i, f in enumerate(full)
    ]


def check_adequacy(
    network: Network,
    result: PowerFlowResult,
    *,
    band_hz: float = DEFAULT_BAND_HZ,
    risk_participation: float = RISK_PARTICIPATION,
) -> AdequacyReport:
    """Judge this network's chosen model levels against its own full-order model."""
    full_net = _full_order(network)
    full = linearize_network(full_net, result)
    full_modal = analyze(full.A, state_names=full.state_names)

    reduced_net = network
    is_full_order = full_net.models == network.models and all(
        not d.level and not d.states for d in network.der_units
    )
    reduced = full if is_full_order else linearize_network(reduced_net, result)
    # One eigendecomposition per system, reused for the band comparison and
    # the stiffness figure -- on a large case each one costs seconds.
    reduced_modal = full_modal if is_full_order else analyze(reduced.A, state_names=reduced.state_names)

    kept = set(reduced.state_names)
    removed = [n for n in full.state_names if n not in kept]

    # --- the screen: what do the removed states participate in? ---
    band_modes = _band_mode_indices(full_modal, band_hz)
    risks: list[StateRisk] = []
    name_index = {n: i for i, n in enumerate(full.state_names)}
    for state in removed:
        i = name_index[state]
        for j in band_modes:
            p = float(abs(full_modal.participation[i, j]))
            if p < risk_participation:
                continue
            lam = full_modal.eigenvalues[j]
            risks.append(StateRisk(
                state=state,
                group=_group_of(state, network),
                mode_hz=abs(lam) / (2 * np.pi),
                mode_damping_pct=(-lam.real / abs(lam) * 100) if abs(lam) > 0 else 100.0,
                participation=p,
            ))
    risks.sort(key=lambda r: -r.participation)

    # --- the measurement: what did it actually cost? ---
    # How much of each full-model mode belonged to the states being removed.
    # A mode that *is* the removed dynamics is meant to disappear with them.
    removed_rows = [name_index[n] for n in removed]
    share = {
        j: float(np.abs(full_modal.participation[removed_rows, j]).sum()) if removed_rows else 0.0
        for j in band_modes
    }
    modes = _pair_up(
        _modes_in_band(full_modal, band_hz), _modes_in_band(reduced_modal, band_hz), share
    )
    matched = [m for m in modes if m.matched]
    max_d_hz = max((m.d_hz for m in matched), default=0.0)
    max_d_damp = max((m.d_damping_pct for m in matched), default=0.0)
    expected = sum(1 for m in modes if m.expected_loss)
    unmatched = sum(1 for m in modes if not m.matched and not m.expected_loss)

    def stiffness(modal) -> float:
        ev = modal.eigenvalues
        return float(np.abs(ev).max()) if len(ev) else 0.0

    notes: list[str] = []
    if is_full_order:
        verdict = "full_order"
        notes.append("Nothing is reduced -- this network runs the full model.")
    elif unmatched:
        verdict = "unsafe"
        notes.append(
            f"{unmatched} mode(s) of the full model have no counterpart in the reduced one within "
            f"{band_hz:g} Hz, and they were not the removed states' own: the reduction has taken "
            "out dynamics this case depends on."
        )
    elif max_d_hz > 0.05 or max_d_damp > 2.0:
        verdict = "check"
        notes.append(
            f"Modes move by up to {max_d_hz * 1000:.0f} mHz and {max_d_damp:.1f} points of damping. "
            "Acceptable for a screening study, probably not for a damping-controller design."
        )
    else:
        verdict = "safe"
        notes.append(
            f"Every mode below {band_hz:g} Hz is reproduced to within {max_d_hz * 1000:.0f} mHz and "
            f"{max_d_damp:.2f} points of damping."
        )
    if expected:
        notes.append(
            f"{expected} further mode(s) disappeared, but each was made of the states being "
            "removed -- those are the reduction doing what was asked of it, not a loss."
        )
    if risks and verdict in ("safe", "check"):
        worst = risks[0]
        notes.append(
            f"Screen flagged {len({r.state for r in risks})} removed state(s) with participation in "
            f"the band -- worst is {worst.state} at {worst.participation * 100:.0f}% in the "
            f"{worst.mode_hz:.2f} Hz mode. The spectra still agree, so this is a caution, not a defect."
        )
    sf, sr = stiffness(full_modal), stiffness(reduced_modal)
    if not is_full_order:
        if sr > 0 and sf / sr >= 2:
            notes.append(
                f"Fastest mode drops from {sf:.2e} to {sr:.2e} rad/s ({sf / sr:.0f}x less stiff) "
                f"for {len(removed)} fewer states -- that ratio, not the state count, is what "
                "decides how long a time-domain run takes."
            )
        else:
            notes.append(
                f"{len(removed)} fewer states, but the fastest mode barely moves "
                f"({sf:.2e} -> {sr:.2e} rad/s): something still in the model is as fast as what "
                "was removed, so a time-domain run will not get much quicker."
            )

    return AdequacyReport(
        band_hz=band_hz,
        n_states_full=full.A.shape[0],
        n_states_reduced=reduced.A.shape[0],
        removed_states=removed,
        risks=risks,
        modes=modes,
        max_d_hz=max_d_hz,
        max_d_damping_pct=max_d_damp,
        unmatched=unmatched,
        stiffness_full=sf,
        stiffness_reduced=sr,
        verdict=verdict,
        notes=notes,
    )
