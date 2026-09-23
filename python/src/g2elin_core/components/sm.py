"""Synchronous machine + transformer + AVR/PSS/governor, ported from
``Functions/symSG.m`` (nonlinear DAE + linearization) and ``Functions/SG_subs.m``
(numeric operating-point substitution).

Two known "use index 1 for all units" simplifications carry over unchanged
from the MATLAB source (``script_generic.m``), preserved here rather than
silently fixed, since fixing them would make results incomparable to a
MATLAB-generated reference:

- Every machine's transformer impedance (``Rt``/``Lt``) comes from the
  *first* DG's ``Y_DER`` TR_R/TR_XL columns, not its own — see
  ``script_generic.m``'s ``SG_Data`` construction (``Y_TR(1,4) Y_TR(1,5)``
  for every row `i`). **No longer the default:** each unit now uses its own
  transformer (the one the power flow uses, see
  ``operating_point.unit_transformer_rx``); set
  ``Network.units_use_first_transformer`` to reproduce the MATLAB behaviour.
- Every node's shunt susceptance (``Cl``, in :mod:`g2elin_core.components.node`)
  comes from the *first* line's B column, not the sum of its own incident
  lines' — see ``script_generic.m``'s node-substitution loop (``Cl = Y_line(1,6)``).
"""

from __future__ import annotations

import cmath
from collections.abc import Mapping
from functools import lru_cache

import numpy as np
import sympy as sp

from .base import (
    ComponentDAE, LinearComponent, NonlinearFuncs, NonlinearJacobians,
    apply_reduction, build_dae, equilibrium_subs, mode_key,
)

# Parameters
wb, Rt, Lt, Rg, Ra, Ll, Lad, Laq, Lfd, Rfd = sp.symbols("wb Rt Lt Rg Ra Ll Lad Laq Lfd Rfd")
L1d, L1q, L2q, R1d, R1q, R2q, H, KD, mp, TG = sp.symbols("L1d L1q L2q R1d R1q R2q H KD mp TG")
T_LP, T_HP, K_PSS, T1n, T1d, T2n, T2d = sp.symbols("T_LP T_HP K_PSS T1n T1d T2n T2d")
Tr, Ka, Ta, Ke, Te, Kfd, Tfd = sp.symbols("Tr Ka Ta Ke Te Kfd Tfd")

# States (PSS on)
igd, igq, phi_d, phi_q, phi_fd, phi_1d, phi_1q, phi_2q = sp.symbols(
    "igd igq phi_d phi_q phi_fd phi_1d phi_1q phi_2q"
)
Dwr, theta, Pm, Dw1, v1, v2, vpss, e1, e2, efd, e3 = sp.symbols(
    "Dwr theta Pm Dw1 v1 v2 vpss e1 e2 efd e3"
)

# Algebraic
wr, ved, veq, id_, i1d, ifd, iq, i1q, i2q, vgd, vgq, Cm, Ce, DP, Et = sp.symbols(
    "wr ved veq id_ i1d ifd iq i1q i2q vgd vgq Cm Ce DP Et"
)

# Inputs
v_ref, p_ref, w_ref, theta_g, vgd_g, vgq_g = sp.symbols("v_ref p_ref w_ref theta_g vgd_g vgq_g")

_STATE_NAMES_PSS_ON = [
    "i_gd", "i_gq", "psi_d", "psi_q", "psi_fd", "psi_1d", "psi_1q", "psi_2q",
    "dw_r", "theta", "P_m", "dw_1", "v_1", "v_2", "v_pss", "e_1", "e_2", "e_fd", "e_3",
]


@lru_cache(maxsize=16)
def sm_dae(is_slack: bool, modes: tuple[tuple[str, str], ...] = ()) -> ComponentDAE:
    """``modes`` is a :func:`~g2elin_core.components.base.mode_key` tuple
    naming the states this machine gives up -- see
    :mod:`g2elin_core.reduction` for the catalogue and the named orders it
    builds out of them. The default ``()`` is the full 19-state model.
    """
    state_vec = [igd, igq, phi_d, phi_q, phi_fd, phi_1d, phi_1q, phi_2q, Dwr, theta, Pm,
                 Dw1, v1, v2, vpss, e1, e2, efd, e3]
    alg_vec = [wr, ved, veq, id_, i1d, ifd, iq, i1q, i2q, vgd, vgq, Cm, Ce, DP, Et]
    us_vec = [v_ref, p_ref, w_ref]
    ug_vec = [vgd_g, vgq_g] if is_slack else [theta_g, vgd_g, vgq_g]
    input_vec = us_vec + ug_vec

    # --- Differential equations (physical layer + electrical + mechanical + AVR/PSS) ---
    digd = (wb / Lt) * (ved - vgd - Rt * igd + wr * Lt * igq)
    digq = (wb / Lt) * (veq - vgq - Rt * igq - wr * Lt * igd)

    dphi_d = wb * (ved + Ra * id_ + wr * phi_q)
    dphi_q = wb * (veq + Ra * iq - wr * phi_d)
    dphi_fd = (wb * Rfd / Lad) * efd - wb * Rfd * ifd
    dphi_1d = wb * (-R1d * i1d)
    dphi_1q = wb * (-R1q * i1q)
    dphi_2q = wb * (-R2q * i2q)

    dDwr = (Cm - Ce - KD * Dwr) / (2 * H)
    dtheta = wb * wr

    dPm = (p_ref - DP - Pm) / TG

    dDw1 = (1 / T_LP) * (Dwr - Dw1)
    dv1 = (1 / T_HP) * (T_HP * K_PSS * dDw1 - v1)
    dv2 = (1 / T1d) * (T1n * dv1 + v1 - v2)
    dvpss = (1 / T2d) * (T2n * dv2 + v2 - vpss)

    de1 = (Et - e1) / Tr
    de2 = (Ka * (v_ref - e1 - e3 + vpss) - e2) / Ta
    defd = (Ke * e2 - efd) / Te
    de3 = (Kfd * defd - e3) / Tfd

    diffeq_vec = [digd, digq, dphi_d, dphi_q, dphi_fd, dphi_1d, dphi_1q, dphi_2q, dDwr,
                  dtheta, dPm, dDw1, dv1, dv2, dvpss, de1, de2, defd, de3]

    # --- Algebraic constraints ---
    invd = sp.Matrix([[-Lad - Ll, Lad, Lad], [-Lad, L1d + Lad, Lad], [-Lad, Lad, Lfd + Lad]]).inv()
    invq = sp.Matrix([[-Laq - Ll, Laq, Laq], [-Laq, L1q + Laq, Laq], [-Laq, Laq, L2q + Laq]]).inv()
    currents_d = invd * sp.Matrix([phi_d, phi_1d, phi_fd])
    currents_q = invq * sp.Matrix([phi_q, phi_1q, phi_2q])

    alg = [
        wr - Dwr - w_ref,
        ved - (id_ - igd) * Rg,
        veq - (iq - igq) * Rg,
        id_ - currents_d[0],
        i1d - currents_d[1],
        ifd - currents_d[2],
        iq - currents_q[0],
        i1q - currents_q[1],
        i2q - currents_q[2],
    ]
    if is_slack:
        alg += [vgd - vgd_g, vgq - vgq_g]
    else:
        alg += [
            vgd - vgd_g * sp.cos(theta_g - theta) + vgq_g * sp.sin(theta_g - theta),
            vgq - vgd_g * sp.sin(theta_g - theta) - vgq_g * sp.cos(theta_g - theta),
        ]
    alg += [
        Cm - Pm / wr,
        Ce - (phi_d * iq - phi_q * id_),
        DP - (wr - w_ref) / mp,
        Et - sp.sqrt(ved**2 + veq**2),
    ]

    # --- Outputs ---
    if is_slack:
        outputeq_g = [igd, igq, theta, wr]
    else:
        igd_g_out = igd * sp.cos(theta - theta_g) - igq * sp.sin(theta - theta_g)
        igq_g_out = igd * sp.sin(theta - theta_g) + igq * sp.cos(theta - theta_g)
        outputeq_g = [igd_g_out, igq_g_out]

    p_e = ved * igd + veq * igq
    q_e = -ved * igq + veq * igd
    outputeq_s = [p_e, q_e, wr, theta, Dwr, Et]
    output_vec = outputeq_s + outputeq_g

    state_vec, diffeq_vec, state_names, alg_vec, alg, output_vec = apply_reduction(
        states=list(zip(state_vec, diffeq_vec, _STATE_NAMES_PSS_ON)),
        modes=dict(modes),
        alg_vec=alg_vec,
        algeq_vec=alg,
        output_vec=output_vec,
    )

    return build_dae(
        state_vec=state_vec,
        alg_vec=alg_vec,
        input_vec=input_vec,
        diffeq_vec=diffeq_vec,
        algeq_vec=alg,
        output_vec=output_vec,
        n_us=3,
        n_ug=len(ug_vec),
        n_out_s=6,
        n_out_g=len(outputeq_g),
        state_names=state_names,
        input_names=["V_ref", "P_ref", "w_ref"],
        output_names=["p_e", "q_e", "w_r", "theta", "dw_r_dot", "V_t"],
    )


@lru_cache(maxsize=16)
def sm_nonlinear_funcs(is_slack: bool, modes: tuple[tuple[str, str], ...] = ()) -> NonlinearFuncs:
    """Numeric callables for the full nonlinear SM DAE — see
    :meth:`g2elin_core.components.base.ComponentDAE.nonlinear_funcs`.
    """
    return sm_dae(is_slack, modes).nonlinear_funcs()


@lru_cache(maxsize=16)
def sm_nonlinear_jacobians(is_slack: bool, modes: tuple[tuple[str, str], ...] = ()) -> NonlinearJacobians:
    """Numeric Jacobian callables for the SM DAE's algebraic/output
    equations — see :meth:`g2elin_core.components.base.ComponentDAE.nonlinear_jacobians`.
    """
    return sm_dae(is_slack, modes).nonlinear_jacobians()


class SmOperatingPoint:
    """Everything needed to linearize one SM instance, computed from the
    solved power flow the way ``script_generic.m``'s "SM Initializations"
    section + ``SG_subs.m`` do.
    """

    def __init__(
        self,
        *,
        params: dict,
        v_terminal_pu: float,
        angle_terminal_rad: float,
        p_terminal_pu: float,
        q_terminal_pu: float,
        v_grid_pu: float,
        angle_grid_rad: float,
        p_ref_pu: float,
        theta_g_rad: float,
        is_slack: bool,
    ):
        self.p = params
        Ra_ = params["Ra"]
        Lq_ = params["Ll"] + params["Laq"]
        Lad_ = params["Lad"]
        Ll_ = params["Ll"]
        Lfd_ = params["Lfd"]

        v0 = v_terminal_pu * cmath.exp(1j * angle_terminal_rad)
        Ii = (p_terminal_pu - 1j * q_terminal_pu) / v0.conjugate()
        delta0 = cmath.phase(v0 + (Ra_ + 1j * Lq_) * Ii)
        theta0 = delta0 - np.pi / 2

        Ic0 = v0 / params["RL_pu"]
        Ig0 = Ii - Ic0
        # Rotate v0 (already v_terminal_pu*exp(j*angle_terminal_rad)) by -theta0
        # into the machine's own dq frame — NOT v0*exp(j*(angle_terminal_rad-theta0)),
        # which would double-count angle_terminal_rad. That bug was invisible for
        # the slack unit only because its bus angle is pinned to exactly 0 by the
        # power-flow slack reference, masking it for every non-slack SM.
        vdq0 = v0 * cmath.exp(-1j * theta0)
        idq0 = Ii * cmath.exp(-1j * theta0)
        igdq0 = Ig0 * cmath.exp(-1j * theta0)

        self.theta0 = theta0
        self.ved0, self.veq0 = vdq0.real, vdq0.imag
        self.id0, self.iq0 = idq0.real, idq0.imag
        self.igd0, self.igq0 = igdq0.real, igdq0.imag

        self.psi_d0 = self.veq0 + Ra_ * self.iq0
        self.psi_q0 = -self.ved0 - Ra_ * self.id0
        self.ifd0 = (self.veq0 + Ra_ * self.iq0 + (Lad_ + Ll_) * self.id0) / Lad_
        self.efd0 = Lad_ * self.ifd0
        self.psi_fd0 = (Lad_ + Lfd_) * self.ifd0 - Lad_ * self.id0
        self.psi_1d0 = Lad_ * (self.ifd0 - self.id0)
        self.psi_1q0 = -params["Laq"] * self.iq0
        self.psi_2q0 = -params["Laq"] * self.iq0

        # Grid-side voltage (algebraic vars vgd/vgq): the raw node's voltage
        # in THIS machine's own dq frame (rotate by its own theta0).
        grid_v = v_grid_pu * cmath.exp(1j * angle_grid_rad)
        vg_local = grid_v * cmath.exp(-1j * theta0)
        self.vgd0, self.vgq0 = vg_local.real, vg_local.imag

        # vgd_g/vgq_g (an input): the same raw voltage in the GLOBAL dq frame.
        # For the slack unit theta_g is defined as its own theta0, so this
        # coincides with vgd0/vgq0 above.
        vg_global = grid_v * cmath.exp(-1j * theta_g_rad) if not is_slack else vg_local
        self.vgd_g0, self.vgq_g0 = vg_global.real, vg_global.imag
        self.theta_g0 = theta_g_rad
        self.is_slack = is_slack

        self.p_ref0 = p_ref_pu
        self.w_ref0 = 1.0
        self.wr0 = 1.0
        self.Dwr0 = 0.0
        self.Pm0 = p_ref_pu
        self.Dw1_0 = 0.0
        self.v1_0 = 0.0
        self.v2_0 = 0.0
        self.vpss0 = 0.0
        self.Et0 = abs(vdq0)
        self.e1_0 = self.Et0
        self.e3_0 = 0.0
        self.v_ref0 = self.efd0 / params["Ka"] + self.e1_0
        self.e2_0 = params["Ka"] * (self.v_ref0 - self.e1_0 - self.e3_0 + self.vpss0)

        # i1d, i1q, i2q (damper-winding currents) never appear in a nonlinear
        # (product) term anywhere in the model, so they're irrelevant to the
        # *linearized* A/B/C/D matrices and were skipped above. The full
        # nonlinear algebraic residual (alg4-alg9) does need them, though —
        # solve the same invd/invq relation symSG.m uses, which also
        # cross-checks id0/iq0/ifd0 against the independently-derived values
        # above (they must match at a true equilibrium).
        Lad_, Ll_, L1d_ = params["Lad"], params["Ll"], params["L1d"]
        Lfd_, Laq_, L1q_, L2q_ = params["Lfd"], params["Laq"], params["L1q"], params["L2q"]
        invd = np.linalg.inv([[-Lad_ - Ll_, Lad_, Lad_], [-Lad_, L1d_ + Lad_, Lad_], [-Lad_, Lad_, Lfd_ + Lad_]])
        invq = np.linalg.inv([[-Laq_ - Ll_, Laq_, Laq_], [-Laq_, L1q_ + Laq_, Laq_], [-Laq_, Laq_, L2q_ + Laq_]])
        id_check, self.i1d0, ifd_check = invd @ [self.psi_d0, self.psi_1d0, self.psi_fd0]
        iq_check, self.i1q0, self.i2q0 = invq @ [self.psi_q0, self.psi_1q0, self.psi_2q0]
        for name, computed, independent in (
            ("id0", id_check, self.id0), ("iq0", iq_check, self.iq0), ("ifd0", ifd_check, self.ifd0)
        ):
            if not np.isclose(computed, independent, atol=1e-8, rtol=1e-6):
                raise AssertionError(
                    f"SM operating point inconsistent: {name} = {independent} from the power-flow "
                    f"derivation but {computed} from the flux/current matrix relation"
                )

        self.Cm0 = self.Pm0 / self.wr0
        self.Ce0 = self.psi_d0 * self.iq0 - self.psi_q0 * self.id0
        self.DP0 = (self.wr0 - self.w_ref0) / params["mp"]


def _sm_operating_subs(op: SmOperatingPoint) -> dict:
    p = op.p
    subs = {
        wb: p["wb"], Rt: p["Rt"], Lt: p["Lt"], Rg: p["RL_pu"], Ra: p["Ra"], Ll: p["Ll"],
        Lad: p["Lad"], Laq: p["Laq"], Lfd: p["Lfd"], Rfd: p["Rfd"],
        L1d: p["L1d"], L1q: p["L1q"], L2q: p["L2q"], R1d: p["R1d"], R1q: p["R1q"], R2q: p["R2q"],
        H: p["H"], KD: p["KD"], mp: p["mp"], TG: p["TG"],
        T_LP: p["T_LP"], T_HP: p["T_HP"], K_PSS: p["K_PSS"],
        T1n: p["T1n"], T1d: p["T1d"], T2n: p["T2n"], T2d: p["T2d"],
        Tr: p["Tr"], Ka: p["Ka"], Ta: p["Ta"], Ke: p["Ke"], Te: p["Te"], Kfd: p["Kfd"], Tfd: p["Tfd"],
        igd: op.igd0, igq: op.igq0, phi_d: op.psi_d0, phi_q: op.psi_q0, phi_fd: op.psi_fd0,
        phi_1d: op.psi_1d0, phi_1q: op.psi_1q0, phi_2q: op.psi_2q0, Dwr: op.Dwr0, theta: op.theta0,
        Pm: op.Pm0, Dw1: op.Dw1_0, v1: op.v1_0, v2: op.v2_0, vpss: op.vpss0,
        e1: op.e1_0, e2: op.e2_0, efd: op.efd0, e3: op.e3_0,
        wr: op.wr0, ved: op.ved0, veq: op.veq0, id_: op.id0, i1d: op.i1d0, ifd: op.ifd0,
        iq: op.iq0, i1q: op.i1q0, i2q: op.i2q0, vgd: op.vgd0, vgq: op.vgq0,
        Cm: op.Cm0, Ce: op.Ce0, DP: op.DP0, Et: op.Et0,
        v_ref: op.v_ref0, p_ref: op.p_ref0, w_ref: op.w_ref0,
        vgd_g: op.vgd_g0, vgq_g: op.vgq_g0,
    }
    if not op.is_slack:
        subs[theta_g] = op.theta_g0
    return equilibrium_subs(subs)


def linearize_sm(op: SmOperatingPoint, modes: Mapping[str, str] | None = None) -> LinearComponent:
    return sm_dae(op.is_slack, mode_key(modes)).linearize(_sm_operating_subs(op))


def sm_nonlinear_point(op: SmOperatingPoint, modes: Mapping[str, str] | None = None) -> tuple:
    """``(x0, z0, u0, p0)`` for :func:`sm_nonlinear_funcs`, at the same
    operating point :func:`linearize_sm` uses.
    """
    return sm_dae(op.is_slack, mode_key(modes)).point_from_subs(_sm_operating_subs(op))
