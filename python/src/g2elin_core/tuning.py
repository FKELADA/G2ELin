"""Control-loop tuning of the GFM/GFL converter models: response time and
damping <-> controller gains, for each loop.

The default gains in :func:`~g2elin_core.operating_point.gfm_params` /
:func:`~g2elin_core.operating_point.gfl_params` come from pole placement
with ``wn = 3/(zeta*t_r)`` (``script_generic.m``). This module exposes the
same formulas per loop, in both directions, so a loop can be specified by
its response time ``t_r`` and damping ``zeta`` instead of raw gains -- used
by the root-locus sweep (``g2elin_api.sweep``) and mirrored by the web UI's
loop tuner (``web/js/unit_params.js``).

Loop kinds:

- ``pi2``: PI loop on a first-order plant ``K/(1 + tau s)``:
  ``Kp = (2 zeta wn tau - 1)/K``, ``Ki = wn^2 tau / K``
- ``pi0``: PI loop on an integrator of gain ``k``: ``Kp = 2 zeta wn / k``, ``Ki = wn^2 / k``
- ``p1``: proportional loop on a first-order plant: ``Kp = (3 tau / t_r - 1) G``
- ``i1``: integral-only loop: ``Ki = -3 / t_r``
- ``droop``: power loop with droop ``mp`` and emulated inertia ``H``:
  ``wf = 1/(2 mp H)`` (filter time constant ``Tf = 1/wf``)
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable


@dataclass(frozen=True)
class Loop:
    id: str
    name: str
    kind: str  # "pi2" | "pi0" | "p1" | "i1" | "droop"
    kp: str | None = None
    ki: str | None = None
    # Plant constants from the unit's (effective) parameter dict.
    K: Callable[[dict], float] | None = None
    tau: Callable[[dict], float] | None = None
    k: Callable[[dict], float] | None = None
    G: Callable[[dict], float] | None = None


_CURRENT = dict(kind="pi2", kp="KpCL", ki="KiCL", K=lambda p: 1 / p["Rf"], tau=lambda p: p["Lf"] / (p["wb"] * p["Rf"]))

LOOPS: dict[str, list[Loop]] = {
    "gfm": [
        Loop("cl", "Current loop (inner)", **_CURRENT),
        Loop("vl", "Voltage loop (outer)", "pi0", "KpVL", "KiVL", k=lambda p: p["wb"] / p["Cf"]),
        Loop("dc", "DC-link voltage loop", "p1", "Kpdc", G=lambda p: p["Gdc"], tau=lambda p: p["Cdc"] / (p["wb"] * p["Gdc"])),
        Loop("droop", "Droop & emulated inertia (power loop)", "droop"),
    ],
    "gfl": [
        Loop("cl", "Current loop (inner)", **_CURRENT),
        Loop("pll", "Phase-locked loop", "pi0", "Kppll", "Kipll", k=lambda p: p["wb"]),
        Loop("dcv", "DC-voltage loop (outer)", "pi2", "Kpd", "Kid", K=lambda p: -1 / p["Gdc"],
             tau=lambda p: p["Cdc"] / (p["wb"] * p["Gdc"])),
        Loop("q", "Reactive-power loop (outer)", "i1", ki="Kiq"),
    ],
}


def loop(unit_type: str, loop_id: str) -> Loop:
    for lp in LOOPS.get(unit_type, []):
        if lp.id == loop_id:
            return lp
    raise KeyError(f"a {unit_type} unit has no control loop {loop_id!r}")


def tuning_of(lp: Loop, p: dict) -> dict:
    """Current response time ``tr`` [s] and damping ``zeta`` of loop ``lp``
    under parameters ``p`` (``zeta`` only for PI loops); for the droop loop,
    ``H`` [s] and ``Tf`` [s]. Raises ``ValueError`` when the gains don't
    correspond to a positively damped tuning."""
    if lp.kind == "pi2":
        K, tau = lp.K(p), lp.tau(p)
        wn2 = p[lp.ki] * K / tau
        if wn2 <= 0:
            raise ValueError(f"{lp.name}: {lp.kp}/{lp.ki} don't correspond to a tuning")
        wn = math.sqrt(wn2)
        zeta = (p[lp.kp] * K + 1) / (2 * wn * tau)
    elif lp.kind == "pi0":
        k = lp.k(p)
        wn2 = p[lp.ki] * k
        if wn2 <= 0:
            raise ValueError(f"{lp.name}: {lp.kp}/{lp.ki} don't correspond to a tuning")
        wn = math.sqrt(wn2)
        zeta = p[lp.kp] * k / (2 * wn)
    elif lp.kind == "p1":
        d = p[lp.kp] / lp.G(p) + 1
        if d <= 0:
            raise ValueError(f"{lp.name}: {lp.kp} doesn't correspond to a tuning")
        return {"tr": 3 * lp.tau(p) / d}
    elif lp.kind == "i1":
        if p[lp.ki] >= 0:
            raise ValueError(f"{lp.name}: {lp.ki} doesn't correspond to a tuning")
        return {"tr": -3 / p[lp.ki]}
    else:  # droop
        return {"H": 1 / (2 * p["mp"] * p["wf"]), "Tf": 1 / p["wf"]}
    if zeta <= 0:
        raise ValueError(f"{lp.name}: the gains give a non-positive damping")
    return {"tr": 3 / (zeta * wn), "zeta": zeta}


def gains_for(lp: Loop, p: dict, **target: float) -> dict[str, float]:
    """Parameter values giving loop ``lp`` the requested tuning. Anything not
    in ``target`` keeps its current value under ``p``: e.g.
    ``gains_for(lp, p, tr=0.02)`` changes the response time at unchanged
    damping. Targets: ``tr``/``zeta`` (PI and first-order loops), ``H`` or
    ``Tf`` (droop loop; ``mp`` is kept)."""
    if lp.kind == "droop":
        if "Tf" in target:
            return {"wf": 1 / target["Tf"]}
        if "H" in target:
            return {"wf": 1 / (2 * p["mp"] * target["H"])}
        raise ValueError("droop loop: give H or Tf")
    now = tuning_of(lp, p)
    tr = target.get("tr", now["tr"])
    zeta = target.get("zeta", now.get("zeta", 1.0))
    if tr <= 0 or zeta <= 0:
        raise ValueError("response time and damping must be positive")
    wn = 3 / (zeta * tr)
    if lp.kind == "pi2":
        K, tau = lp.K(p), lp.tau(p)
        return {lp.kp: (2 * zeta * wn * tau - 1) / K, lp.ki: wn * wn * tau / K}
    if lp.kind == "pi0":
        k = lp.k(p)
        return {lp.kp: 2 * zeta * wn / k, lp.ki: wn * wn / k}
    if lp.kind == "p1":
        return {lp.kp: (3 * lp.tau(p) / tr - 1) * lp.G(p)}
    return {lp.ki: -3 / tr}  # i1


# Sweepable tuning quantities: name -> (target keyword, scale to SI).
TUNING_QUANTITIES = {"tr_ms": ("tr", 1e-3), "zeta": ("zeta", 1.0), "H": ("H", 1.0), "Tf_ms": ("Tf", 1e-3)}


def quantities_for(lp: Loop) -> list[str]:
    if lp.kind in ("pi2", "pi0"):
        return ["tr_ms", "zeta"]
    if lp.kind in ("p1", "i1"):
        return ["tr_ms"]
    return ["H", "Tf_ms"]
