"""Assembles per-component linear state-spaces into one closed-loop system.

Functionally ports ``Functions/assoc_matrices.m`` + the concatenation/
closed-loop formula in ``script_generic.m``'s "Section IX" — but instead of
replicating ~300 lines of hand-rolled row/column index arithmetic (verified
once against a MATLAB run, near-impossible to verify blind), this builds
the same F/G/K/L interconnection from an explicit list of named-port wiring
rules, checked against ``assoc_matrices.m``'s Parts I-VII by inspection:

- Part II/IV (theta, omega): every non-slack unit's ``theta_g`` input, and
  every node/line/load's ``wg`` input, equals the slack unit's own
  ``theta``/``wr`` output — the reference frame rotates at the slack's speed.
- Part III (DG voltage): each DG's ``vgd_g``/``vgq_g`` input equals the
  voltage output of the raw node its transformer connects to.
- Part V (node current balance): a node's injected current equals the sum
  of incident DG currents (+), "to" line currents (+), "from" line
  currents (-), and load currents (-).
- Part VI (line voltages): a line's two ``vgdj_g``/``vgdk_g`` inputs equal
  its "from"/"to" node's voltage output.
- Part VII (load voltage): a load's ``vgd_g``/``vgq_g`` input equals its
  node's voltage output.

The closed-loop formula itself (``A_tot0 = A_ol + B_ol*G*E_ol*C_ol`` etc.)
is ported directly — it's a single unambiguous line, not index bookkeeping.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Protocol

import numpy as np

# Local port layouts, matching the [us | ug] input order and [s | g] output
# order each components/*.py module was written with.
_UG_PORTS = {
    "sm_slack": {"vgd_g": 0, "vgq_g": 1},
    "ib_slack": {"vgd_g": 0, "vgq_g": 1},  # matches ib_dae()'s ug_vec = [vgd_g, vgq_g]
    # The frame itself has no inputs; a unit that doesn't own it reads it.
    "frame": {"w_in": 0},
    "ib": {"theta_g": 0, "wg": 1, "vgd_g": 2, "vgq_g": 3},
    "sm": {"theta_g": 0, "vgd_g": 1, "vgq_g": 2},
    "gfm": {"theta_g": 0, "vgd_g": 1, "vgq_g": 2},  # GFM is never the slack (see components/gfm.py)
    "gfl": {"theta_g": 0, "vgd_g": 1, "vgq_g": 2},  # GFL is never the slack
    "node": {"wg": 0, "ishd_g": 1, "ishq_g": 2},
    "line": {"wg": 0, "vgdj_g": 1, "vgqj_g": 2, "vgdk_g": 3, "vgqk_g": 4},
    "load": {"wg": 0, "vgd_g": 1, "vgq_g": 2},
}
_OUTG_PORTS = {
    "sm_slack": {"igd_g": 0, "igq_g": 1, "theta": 2, "wr": 3},
    # matches ib_dae()'s output_vec = [p_up, q_up, igd_g, igq_g, theta_up, wup]
    # (p_up/q_up are the n_out_s=2 "own" outputs, skipped here); an IB is
    # always the slack (see components/ib.py), so this is the only IB kind.
    "ib_slack": {"igd_g": 0, "igq_g": 1, "theta": 2, "wr": 3},
    "frame": {"theta": 0, "wr": 1},
    "ib": {"igd_g": 0, "igq_g": 1},
    "sm": {"igd_g": 0, "igq_g": 1},
    "gfm": {"igd_g": 0, "igq_g": 1},
    "gfl": {"igd_g": 0, "igq_g": 1},
    "node": {"vgd_g": 0, "vgq_g": 1},
    "line": {"ild_g": 0, "ilq_g": 1},
    "load": {"icd_g": 0, "icq_g": 1},
}


class PortSpec(Protocol):
    """Just the sizes ``Block``'s topology needs — both ``LinearComponent``
    and the nonlinear-simulation wrapper in ``timedomain/`` satisfy this
    structurally, so the same ``Block``/``Wiring``/topology code works for
    linearized (``assemble``) and nonlinear (``timedomain.emt``) assembly.
    """

    n_us: int
    n_ug: int
    n_out_s: int
    n_out_g: int


@dataclass
class Block:
    name: str
    kind: str  # key into _UG_PORTS / _OUTG_PORTS
    comp: PortSpec
    state_off: int = 0
    input_off: int = 0
    output_off: int = 0

    def input_row(self, port: str) -> int:
        return self.input_off + self.comp.n_us + _UG_PORTS[self.kind][port]

    def output_col(self, port: str) -> int:
        ports = _OUTG_PORTS[self.kind]
        if port in ports:
            return self.output_off + self.comp.n_out_s + ports[port]
        # One of the block's own named outputs (its "s" outputs) -- how the
        # reference frame reads the speed of the machine it follows.
        return self.output_off + list(self.comp.output_names).index(port)


@dataclass
class AssembledSystem:
    A: np.ndarray
    B: np.ndarray
    C: np.ndarray
    D: np.ndarray
    state_names: list[str]
    input_names: list[str]  # exogenous (F-mapped) inputs, in block order
    output_names: list[str]  # exogenous (L-mapped) outputs, in block order


@dataclass
class Wiring:
    """One row: an input port that must equal a (signed) sum of output ports."""

    input_block: Block
    input_port: str
    terms: list[tuple[Block, str, float]] = field(default_factory=list)


@dataclass
class Topology:
    """Block offsets + the interconnection's F/G selector matrices — the
    part of assembly that's identical whether the blocks are linearized
    (``assemble``) or nonlinear (``timedomain.emt``): it only encodes *which
    port connects to which*, never any component's actual dynamics.
    """

    blocks: list[Block]
    n_x: int
    n_u: int
    n_y: int
    G: np.ndarray  # (n_u, n_y): full input vector = F @ u_exo + G @ y
    F: np.ndarray  # (n_u, n_exo_in): selects each block's exogenous (us) inputs
    input_names: list[str]  # exogenous input names, in block order (F's columns)


def compute_topology(blocks: list[Block], wiring: list[Wiring]) -> Topology:
    off_x = off_u = off_y = 0
    for b in blocks:
        b.state_off, b.input_off, b.output_off = off_x, off_u, off_y
        off_x += getattr(b.comp, "n_states", 0)
        off_u += b.comp.n_us + b.comp.n_ug
        off_y += b.comp.n_out_s + b.comp.n_out_g
    n_x, n_u, n_y = off_x, off_u, off_y

    G = np.zeros((n_u, n_y))
    for w in wiring:
        row = w.input_block.input_row(w.input_port)
        for term_block, term_port, coeff in w.terms:
            G[row, term_block.output_col(term_port)] += coeff

    # F: selects each block's exogenous (us) inputs out of the full input vector.
    n_exo_in = sum(b.comp.n_us for b in blocks)
    F = np.zeros((n_u, n_exo_in))
    input_names: list[str] = []
    exo_col = 0
    for b in blocks:
        for i in range(b.comp.n_us):
            F[b.input_off + i, exo_col] = 1.0
            input_names.append(f"{b.comp.input_names[i]}_{{{b.name}}}")
            exo_col += 1

    return Topology(blocks=blocks, n_x=n_x, n_u=n_u, n_y=n_y, G=G, F=F, input_names=input_names)


def assemble(blocks: list[Block], wiring: list[Wiring]) -> AssembledSystem:
    topo = compute_topology(blocks, wiring)
    n_x, n_u, n_y = topo.n_x, topo.n_u, topo.n_y

    A_ol = _blkdiag([b.comp.A for b in blocks], n_x, n_x)
    B_ol = _blkdiag([b.comp.B for b in blocks], n_x, n_u)
    C_ol = _blkdiag([b.comp.C for b in blocks], n_y, n_x)
    D_ol = _blkdiag([b.comp.D for b in blocks], n_y, n_u)

    # L: selects each block's exogenous (s) outputs out of the full output vector.
    n_exo_out = sum(b.comp.n_out_s for b in blocks)
    L = np.zeros((n_exo_out, n_y))
    output_names: list[str] = []
    exo_row = 0
    for b in blocks:
        for i in range(b.comp.n_out_s):
            L[exo_row, b.output_off + i] = 1.0
            output_names.append(f"{b.comp.output_names[i]}_{{{b.name}}}")
            exo_row += 1

    n_exo_in = topo.F.shape[1]
    K = np.zeros((n_exo_out, n_exo_in))

    E_ol = np.linalg.inv(np.eye(n_y) - D_ol @ topo.G)
    A_tot = A_ol + B_ol @ topo.G @ E_ol @ C_ol
    B_tot = (B_ol @ topo.G @ E_ol @ D_ol + B_ol) @ topo.F
    C_tot = L @ E_ol @ C_ol
    D_tot = L @ E_ol @ D_ol @ topo.F + K

    state_names = [f"{n}_{{{b.name}}}" for b in blocks for n in b.comp.state_names]

    return AssembledSystem(A_tot, B_tot, C_tot, D_tot, state_names, topo.input_names, output_names)


def _blkdiag(mats: list[np.ndarray], n_rows: int, n_cols: int) -> np.ndarray:
    out = np.zeros((n_rows, n_cols))
    r = c = 0
    for m in mats:
        rr, cc = m.shape
        out[r : r + rr, c : c + cc] = m
        r += rr
        c += cc
    return out
