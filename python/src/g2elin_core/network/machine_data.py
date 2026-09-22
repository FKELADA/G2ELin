"""Turning published synchronous-machine data into this model's parameters.

Machine data is published the way a test engineer measures it -- the
synchronous, transient and subtransient reactances, and the open-circuit time
constants that go with them -- while ``components/sm.py`` is written the way
the machine is built, in mutual and leakage inductances with explicit field and
damper windings. :func:`flux_linkage_params` converts the first into the
second.

The relations are the standard ones (Kundur, ch. 4, unsaturated):

* ``Lad = Xd - Xl`` and ``Laq = Xq - Xl`` -- the mutual inductance is what is
  left of the synchronous reactance once the leakage is taken out;
* ``X' = Xl + Lad||Lfd`` and ``X'' = Xl + Lad||Lfd||L1d`` -- each winding that
  can hold its flux appears in parallel, so the machine looks stiffer the
  faster you look at it;
* ``T'do = (Lad + Lfd) / (wb * Rfd)`` -- each winding's open-circuit time
  constant is its own inductance over its own resistance.

Inverting them, as ``tests/test_kundur.py`` does, gives the published numbers
back exactly.

Everything here is per unit of the *machine's* own rating, which is what
``DerUnit.sn_mva`` exists to declare.
"""

from __future__ import annotations

import math


def _parallel(*values: float) -> float:
    return 1.0 / sum(1.0 / v for v in values)


def flux_linkage_params(
    *,
    xd: float, xq: float, xl: float, ra: float,
    xdp: float, xqp: float, xdpp: float, xqpp: float,
    td0p: float, tq0p: float, td0pp: float, tq0pp: float,
    h: float, kd: float = 0.0, f_hz: float = 60.0,
) -> dict[str, float]:
    """Published machine data -> the parameters ``components/sm.py`` wants.

    All reactances and time constants are per unit of the machine's own rating
    and in seconds; the result is too.
    """
    wb = 2 * math.pi * f_hz
    lad, laq = xd - xl, xq - xl
    lfd = lad * (xdp - xl) / (lad - (xdp - xl))
    l1d = 1.0 / (1.0 / (xdpp - xl) - 1.0 / lad - 1.0 / lfd)
    l1q = laq * (xqp - xl) / (laq - (xqp - xl))
    l2q = 1.0 / (1.0 / (xqpp - xl) - 1.0 / laq - 1.0 / l1q)
    return {
        "Ra": ra, "Ll": xl, "Lad": lad, "Laq": laq,
        "Lfd": lfd, "Rfd": (lad + lfd) / (wb * td0p),
        "L1d": l1d, "R1d": (l1d + _parallel(lad, lfd)) / (wb * td0pp),
        "L1q": l1q, "R1q": (laq + l1q) / (wb * tq0p),
        "L2q": l2q, "R2q": (l2q + _parallel(laq, l1q)) / (wb * tq0pp),
        "H": h, "KD": kd,
    }


# A conventional large round-rotor machine, per unit of its own rating. Used
# for networks whose source data is a *power-flow* case and carries no dynamic
# data at all (the IEEE cases as pandapower ships them). It is a placeholder
# with plausible numbers, not the machine any particular case was studied with:
# results that depend on it are qualitative until real data replaces it.
GENERIC_MACHINE = dict(
    xd=1.8, xq=1.75, xl=0.15, ra=0.003,
    xdp=0.30, xqp=0.55, xdpp=0.22, xqpp=0.22,
    td0p=7.0, tq0p=0.5, td0pp=0.035, tq0pp=0.07,
    h=4.0,
)


def generic_machine_params(f_hz: float = 60.0, h: float | None = None) -> dict[str, float]:
    """:data:`GENERIC_MACHINE` as model parameters, per unit of the machine."""
    spec = dict(GENERIC_MACHINE)
    if h is not None:
        spec["h"] = h
    return flux_linkage_params(**spec, f_hz=f_hz)
