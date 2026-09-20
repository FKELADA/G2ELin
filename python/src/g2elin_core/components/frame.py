"""The common dq reference frame the whole network model is written in.

Every node, line and load equation uses the frame's speed ``wg`` (the
rotational terms ``wg*L*i``, ``wg*C*v``), and every unit transforms between
its own dq frame and the common one through ``theta - theta_g``. *Which*
frame that is, is a modelling choice, not physics: any frame gives the same
physical trajectories, expressed differently.

This component is the frame itself: one state, the frame angle, with no
electrical connection to anything, turning either at a fixed speed

    dtheta_f/dt = wb * w_ref,      outputs: (theta_f, w_ref)

or with a machine it follows (``driven=True``), which is the usual case:

    dtheta_f/dt = wb * w_in,       outputs: (theta_f, w_in)

Either way it feeds exactly the ``theta``/``wr`` ports a slack unit used to
feed (see ``interconnect/network_assembly.py``). The fixed variant is the
angle part an infinite bus already carries (``components/ib.py``:
``dtheta_up = wb*wup`` with ``wup`` at 1 pu) without the voltage source
behind it -- which is why every infinite-bus preset in this repository has
effectively been running on a frame of its own all along.

Why a block of its own: with the frame taken from the slack unit, that unit
could never be disconnected -- the frame would go with it. Standing apart,
every unit is just a unit: any of them can be tripped, a network can be
split into islands (one frame each, following that island's own reference
machine), and the slack keeps only its steady-state meaning, the power
flow's reference. EMT tools have no slack at all for the same reason.

**Why it follows a machine rather than simply running at nominal speed.**
Both give identical trajectories -- the nonlinear models agree to solver
tolerance, since a frame is only a choice of coordinates. But governor and
droop response leave the system at a slightly different frequency after a
disturbance, and in a frame running at exactly nominal speed the network's
dq vectors then turn slowly *forever*: there is no steady state to
linearise about, and modal analysis and the linearised overlay degrade the
longer the horizon (measured on WSCC-9: 5.6 % error at 0.3 s, 17.9 % at
1.5 s, against 4.4 % and 5.1 % for a frame that follows the system). A
frame carried by a machine settles with it, so the equilibrium exists and
small-signal analysis stays sharp.

One free angle remains: nothing pins the frame's absolute position, so a
network with no infinite bus has a zero eigenvalue -- the familiar
reference-angle mode, a property of the formulation, not an instability.
"""

from __future__ import annotations

from functools import lru_cache

import sympy as sp

from .base import ComponentDAE, LinearComponent, NonlinearFuncs, NonlinearJacobians, build_dae

wb, wref = sp.symbols("wb wref")
theta_f, w_in = sp.symbols("theta_f w_in")


@lru_cache(maxsize=2)
def frame_dae(driven: bool = True) -> ComponentDAE:
    """``driven``: the frame turns with a machine's speed (its ``w_in`` port,
    wired to that machine's own speed output); otherwise at the fixed
    ``w_ref``, which is what an island referenced to an infinite bus needs."""
    speed = w_in if driven else wref
    return build_dae(
        state_vec=[theta_f],
        alg_vec=[],
        input_vec=[w_in] if driven else [],
        diffeq_vec=[wb * speed],
        algeq_vec=[],
        output_vec=[theta_f, speed],  # outputs_g only: the ports units and branches read
        n_us=0,
        n_ug=1 if driven else 0,
        n_out_s=0,
        n_out_g=2,
        state_names=["theta"],
        input_names=[],
        output_names=[],
    )


@lru_cache(maxsize=2)
def frame_nonlinear_funcs(driven: bool = True) -> NonlinearFuncs:
    return frame_dae(driven).nonlinear_funcs()


@lru_cache(maxsize=2)
def frame_nonlinear_jacobians(driven: bool = True) -> NonlinearJacobians:
    return frame_dae(driven).nonlinear_jacobians()


def _frame_operating_subs(*, wb_val: float, theta0: float, w0: float = 1.0, driven: bool = True) -> dict:
    subs = {wb: wb_val, theta_f: theta0}
    subs[w_in if driven else wref] = w0
    return subs


def linearize_frame(*, wb_val: float, theta0: float, w0: float = 1.0, driven: bool = True) -> LinearComponent:
    """The frame as a linear block: an integrator of the speed it follows.

    ``theta0`` is where the frame starts, i.e. the operating point's own
    reference angle (``NetworkOperatingPoint.theta_g_rad``), so every unit's
    ``theta - theta_g`` starts at the value the operating point computed.
    """
    return frame_dae(driven).linearize(_frame_operating_subs(wb_val=wb_val, theta0=theta0, w0=w0, driven=driven))


def frame_nonlinear_point(*, wb_val: float, theta0: float, w0: float = 1.0, driven: bool = True) -> tuple:
    """``(x0, z0, u0, p0)`` for :func:`frame_nonlinear_funcs`."""
    return frame_dae(driven).point_from_subs(
        _frame_operating_subs(wb_val=wb_val, theta0=theta0, w0=w0, driven=driven)
    )
