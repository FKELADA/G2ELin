"""Per-unit base-value calculator, ported from ``Functions/PU_calc.m``.

Note the AC bases use an amplitude-invariant (peak-quantity) convention —
``Ub = sqrt(2/3)*Ul``, ``Ib = sqrt(2)*Il`` — rather than the more common
RMS convention, matching the EMT-oriented per-unit system the MATLAB
toolbox uses throughout. The resulting impedance base ``Zb = Ub/Ib`` still
reduces to the familiar ``Un_kV^2 / Sn_MVA`` (the amplitude factors
cancel), which is what actually matters for the R/X/B per-unit conversions
used elsewhere in this package.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class PuBase:
    ub: float
    ib: float
    zb: float
    lb: float
    cb: float
    ub_dc: float
    ib_dc: float
    zb_dc: float
    cb_dc: float


def pu_base(*, wb: float, un_kv: float, sn_mva: float) -> PuBase:
    un = un_kv * 1e3
    sn = sn_mva * 1e6
    ub = math.sqrt(2 / 3) * un
    il = sn / (math.sqrt(3) * un)
    ib = math.sqrt(2) * il
    zb = ub / ib
    lb = zb / wb
    cb = 1 / (zb * wb)

    ub_dc = 2 * ub
    ib_dc = (3 / 4) * ib
    zb_dc = (8 / 3) * zb
    cb_dc = (3 / 8) * cb

    return PuBase(ub=ub, ib=ib, zb=zb, lb=lb, cb=cb, ub_dc=ub_dc, ib_dc=ib_dc, zb_dc=zb_dc, cb_dc=cb_dc)
