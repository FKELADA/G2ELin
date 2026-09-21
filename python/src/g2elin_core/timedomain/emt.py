"""EMT (nonlinear) time-domain simulation — feature 2.2/2.3 from the
migration plan.

This is a dq-frame, averaged-converter, single-rotating-reference-frame
nonlinear time-domain model — not a switching-level, three-phase/abc,
distributed-network model. Whether that counts as "EMT" is a matter of
definition: it integrates the actual nonlinear differential-algebraic
equations (not a linearized small-signal approximation), which is the
property that matters for this project's purposes, so it's named and
scoped as the EMT simulation feature. It doesn't capture switching-level
converter physics, abc-frame/unbalanced faults, or distributed/
traveling-wave line effects.

Reuses exactly the nonlinear ``f``/``g``/``h`` callables built for
cross-validation (:mod:`g2elin_core.components.base`) and the same
interconnection topology used for linear modal analysis
(:mod:`g2elin_core.interconnect`) — the wiring between components (KVL/KCL
in the dq frame) is exact regardless of whether the components themselves
are linearized, so it doesn't need to be re-derived; only the per-component
math changes from "substitute numbers into a Jacobian" to "call a nonlinear
function and Newton-solve the coupled algebraic system".

**What's coupled and why a single global Newton solve is needed.** Each
component's own algebraic variables ``z_i`` are governed by its own
``g_i(x_i, z_i, u_i) = 0``. But a component's inputs ``u_i`` are wired to
*other* components' outputs ``y_j = h_j(x_j, z_j, u_j)`` (the same linear
``F``/``G`` selector matrices ``interconnect.compute_topology`` builds) —
and since ``h_j`` can itself depend on ``u_j``, the full input vector ``u``
and every component's ``z`` have to be solved *together*:

    0 = g(x, z, u)                     (every component's own algebraic eqs)
    u = F @ u_exo + G @ h(x, z, u)      (the interconnection, exact and linear)

at every state ``x`` the ODE integrator asks for — there's no way to solve
one component's ``z_i`` in isolation first and plug it into the next.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Iterator

import numpy as np
import scipy.integrate
import scipy.optimize

from g2elin_core.components.frame import frame_dae, frame_nonlinear_funcs, frame_nonlinear_jacobians, frame_nonlinear_point
from g2elin_core.components.gfl import (
    GflOperatingPoint, gfl_dae, gfl_nonlinear_funcs, gfl_nonlinear_jacobians, gfl_nonlinear_point,
)
from g2elin_core.components.gfm import (
    GfmOperatingPoint, gfm_dae, gfm_nonlinear_funcs, gfm_nonlinear_jacobians, gfm_nonlinear_point,
)
from g2elin_core.components.ib import ib_dae, ib_nonlinear_funcs, ib_nonlinear_jacobians, ib_nonlinear_point
from g2elin_core.components.line import line_dae, line_nonlinear_funcs, line_nonlinear_jacobians, line_nonlinear_point
from g2elin_core.components.load import load_dae, load_nonlinear_funcs, load_nonlinear_jacobians, load_nonlinear_point
from g2elin_core.components.node import node_dae, node_nonlinear_funcs, node_nonlinear_jacobians, node_nonlinear_point
from g2elin_core.components.sm import (
    SmOperatingPoint, sm_dae, sm_nonlinear_funcs, sm_nonlinear_jacobians, sm_nonlinear_point,
)
from g2elin_core.interconnect import Block, build_blocks_and_wiring, compute_topology
from g2elin_core.network.breakers import frame_references, node_capacitances
from g2elin_core.network.schema import Network
from g2elin_core.operating_point import NetworkOperatingPoint, compute_operating_point
from g2elin_core.powerflow import PowerFlowResult

Vec = np.ndarray


@dataclass
class NonlinearBlockComp:
    """A component instance's nonlinear model, bound to its own fixed
    parameters — satisfies :class:`~g2elin_core.interconnect.assemble.PortSpec`
    so it slots into the same ``Block``/``Wiring``/topology machinery
    :mod:`g2elin_core.interconnect` uses for linear assembly.
    """

    n_states: int
    n_z: int
    n_us: int
    n_ug: int
    n_out_s: int
    n_out_g: int
    f: Callable[[Vec, Vec, Vec], Vec]
    g: Callable[[Vec, Vec, Vec], Vec]
    h: Callable[[Vec, Vec, Vec], Vec]
    # Analytic Jacobians of g/h (see NonlinearNetworkModel._residual_jacobian) —
    # what makes the coupled Newton solve robust and fast instead of relying
    # on SciPy's numerically-estimated one.
    Gz: Callable[[Vec, Vec, Vec], Vec]
    Gx: Callable[[Vec, Vec, Vec], Vec]
    Gu: Callable[[Vec, Vec, Vec], Vec]
    Hz: Callable[[Vec, Vec, Vec], Vec]
    Hx: Callable[[Vec, Vec, Vec], Vec]
    Hu: Callable[[Vec, Vec, Vec], Vec]
    state_names: list[str]
    input_names: list[str]
    output_names: list[str]
    x0: Vec
    z0: Vec
    u0: Vec  # full (us, ug) operating-point input, us-prefix is the exogenous default


def _bind(dae, funcs, jacs, point) -> NonlinearBlockComp:
    x0, z0, u0, p0 = point
    return NonlinearBlockComp(
        n_states=len(dae.state_syms),
        n_z=len(dae.alg_syms),
        n_us=dae.n_us,
        n_ug=dae.n_ug,
        n_out_s=dae.n_out_s,
        n_out_g=dae.n_out_g,
        f=lambda x, z, u: funcs.f(x, z, u, p0),
        g=lambda x, z, u: funcs.g(x, z, u, p0),
        h=lambda x, z, u: funcs.h(x, z, u, p0),
        Gz=lambda x, z, u: jacs.Gz(x, z, u, p0),
        Gx=lambda x, z, u: jacs.Gx(x, z, u, p0),
        Gu=lambda x, z, u: jacs.Gu(x, z, u, p0),
        Hz=lambda x, z, u: jacs.Hz(x, z, u, p0),
        Hx=lambda x, z, u: jacs.Hx(x, z, u, p0),
        Hu=lambda x, z, u: jacs.Hu(x, z, u, p0),
        state_names=list(dae.state_names),
        input_names=list(dae.input_names),
        output_names=list(dae.output_names),
        x0=x0,
        z0=z0,
        u0=u0,
    )


def nonlinear_sm_block(op: SmOperatingPoint) -> NonlinearBlockComp:
    return _bind(
        sm_dae(op.is_slack), sm_nonlinear_funcs(op.is_slack), sm_nonlinear_jacobians(op.is_slack),
        sm_nonlinear_point(op),
    )


def nonlinear_gfm_block(op: GfmOperatingPoint) -> NonlinearBlockComp:
    return _bind(gfm_dae(), gfm_nonlinear_funcs(), gfm_nonlinear_jacobians(), gfm_nonlinear_point(op))


def nonlinear_gfl_block(op: GflOperatingPoint) -> NonlinearBlockComp:
    return _bind(gfl_dae(), gfl_nonlinear_funcs(), gfl_nonlinear_jacobians(), gfl_nonlinear_point(op))


def nonlinear_ib_block(**kwargs) -> NonlinearBlockComp:
    is_slack = kwargs.get("is_slack", True)
    return _bind(
        ib_dae(is_slack), ib_nonlinear_funcs(is_slack), ib_nonlinear_jacobians(is_slack),
        ib_nonlinear_point(**kwargs),
    )


def nonlinear_frame_block(**kwargs) -> NonlinearBlockComp:
    driven = kwargs.get("driven", True)
    return _bind(
        frame_dae(driven), frame_nonlinear_funcs(driven), frame_nonlinear_jacobians(driven),
        frame_nonlinear_point(**kwargs),
    )


def nonlinear_line_block(**kwargs) -> NonlinearBlockComp:
    return _bind(line_dae(), line_nonlinear_funcs(), line_nonlinear_jacobians(), line_nonlinear_point(**kwargs))


def nonlinear_node_block(**kwargs) -> NonlinearBlockComp:
    return _bind(node_dae(), node_nonlinear_funcs(), node_nonlinear_jacobians(), node_nonlinear_point(**kwargs))


def nonlinear_load_block(**kwargs) -> NonlinearBlockComp:
    return _bind(load_dae(), load_nonlinear_funcs(), load_nonlinear_jacobians(), load_nonlinear_point(**kwargs))


@dataclass
class NonlinearNetworkModel:
    """The nonlinear counterpart of
    :class:`~g2elin_core.interconnect.assemble.AssembledSystem`: same
    ``Block``s and ``Topology`` (built by the same
    :func:`~g2elin_core.interconnect.build_blocks_and_wiring` +
    :func:`~g2elin_core.interconnect.compute_topology`), but each block's
    ``comp`` is a :class:`NonlinearBlockComp` instead of a
    ``LinearComponent`` — so solving requires the coupled Newton step this
    class's :meth:`solve_algebraic` implements, described in the module
    docstring.
    """

    blocks: list[Block]
    topology: object  # interconnect.Topology
    z_offsets: list[int]
    n_z: int
    # What measurements (timedomain.measurements) need beyond the blocks: the
    # network the model was built from, and the rotation between the
    # power-flow angles and the model's common frame (operating point's
    # theta_g), so measured angles match the power flow's at t = 0.
    network: object = None
    theta_g0: float = 0.0
    # The operating point the blocks were built around (network events
    # rebuild single blocks from it -- see timedomain.events).
    op: object = None
    # Optional warm start for the first algebraic solve (a post-event model
    # starts from the pre-event solution, not from its blocks' own guesses).
    zu_guess: tuple | None = None
    # Optional initial state replacing the blocks' own operating-point one
    # (a post-event model starts from the pre-event state).
    x_init: Vec | None = None

    @property
    def state_names(self) -> list[str]:
        return [f"{n}_{{{b.name}}}" for b in self.blocks for n in b.comp.state_names]

    @property
    def input_names(self) -> list[str]:
        """Each block's own named ("us") inputs -- e.g. an SM's own P_ref,
        not the anonymous interconnection ("ug") port signals wired between
        blocks (see ``interconnect.assemble``'s ``_UG_PORTS``, which have no
        names of their own outside that internal wiring).
        """
        return [f"{n}_{{{b.name}}}" for b in self.blocks for n in b.comp.input_names]

    @property
    def output_names(self) -> list[str]:
        """Each block's own named ("s") outputs -- e.g. an SM's
        ``p_e``/``q_e``/``w_r``/``theta``/``dw_r_dot``/``V_t`` -- the same
        set :class:`~g2elin_core.interconnect.assemble.AssembledSystem`'s
        ``output_names`` exposes for the linear model, since both come from
        the same ``ComponentDAE.output_names``.
        """
        return [f"{n}_{{{b.name}}}" for b in self.blocks for n in b.comp.output_names]

    def initial_state(self) -> Vec:
        if self.x_init is not None:
            return self.x_init.copy()
        parts = [b.comp.x0 for b in self.blocks if b.comp.n_states]
        return np.concatenate(parts) if parts else np.zeros(0)

    def initial_algebraic_guess(self) -> tuple[Vec, Vec]:
        if self.zu_guess is not None:
            return self.zu_guess[0].copy(), self.zu_guess[1].copy()
        z_parts = [b.comp.z0 for b in self.blocks if b.comp.n_z]
        u_parts = [b.comp.u0 for b in self.blocks]
        z0 = np.concatenate(z_parts) if z_parts else np.zeros(0)
        u0 = np.concatenate(u_parts) if u_parts else np.zeros(0)
        return z0, u0

    def default_u_exo(self) -> Vec:
        parts = [b.comp.u0[: b.comp.n_us] for b in self.blocks if b.comp.n_us]
        return np.concatenate(parts) if parts else np.zeros(0)

    def _unpack_x(self, x_all: Vec) -> list[Vec]:
        return [x_all[b.state_off : b.state_off + b.comp.n_states] for b in self.blocks]

    def _unpack_z(self, z_all: Vec) -> list[Vec]:
        return [z_all[zo : zo + b.comp.n_z] for b, zo in zip(self.blocks, self.z_offsets)]

    def _unpack_u(self, u_all: Vec) -> list[Vec]:
        return [u_all[b.input_off : b.input_off + b.comp.n_us + b.comp.n_ug] for b in self.blocks]

    def _residual(self, zu_flat: Vec, x_list: list[Vec], u_exo: Vec) -> Vec:
        z_all, u_all = zu_flat[: self.n_z], zu_flat[self.n_z :]
        z_list, u_list = self._unpack_z(z_all), self._unpack_u(u_all)
        g_parts, y_parts = [], []
        for b, x_i, z_i, u_i in zip(self.blocks, x_list, z_list, u_list):
            if b.comp.n_z:
                g_parts.append(b.comp.g(x_i, z_i, u_i))
            y_parts.append(b.comp.h(x_i, z_i, u_i))
        g_all = np.concatenate(g_parts) if g_parts else np.zeros(0)
        y_all = np.concatenate(y_parts) if y_parts else np.zeros(0)
        u_target = self.topology.F @ u_exo + self.topology.G @ y_all
        return np.concatenate([g_all, u_all - u_target])

    def _residual_jacobian(self, zu_flat: Vec, x_list: list[Vec], u_exo: Vec) -> Vec:
        """Analytic ``d(residual)/d(z, u)`` — see ``solve_algebraic``'s
        docstring for why this replaced a numerically-estimated one.
        Assembled the same way the linear A/B/C/D elimination combines
        per-component Jacobians with the interconnection's ``G`` matrix,
        just kept as an explicit matrix instead of solved in closed form
        (Newton needs it fresh at every iterate, not once at an operating
        point). ``x`` doesn't vary during this solve, so ``Gx``/``Hx``
        aren't needed here — only ``Gz``/``Gu``/``Hz``/``Hu``.
        """
        z_all, u_all = zu_flat[: self.n_z], zu_flat[self.n_z :]
        z_list, u_list = self._unpack_z(z_all), self._unpack_u(u_all)
        n_z, n_u, n_y = self.n_z, self.topology.n_u, self.topology.n_y

        dg_dz = np.zeros((n_z, n_z))
        dg_du = np.zeros((n_z, n_u))
        dh_dz = np.zeros((n_y, n_z))
        dh_du = np.zeros((n_y, n_u))

        for b, x_i, z_i, u_i, zo in zip(self.blocks, x_list, z_list, u_list, self.z_offsets):
            n_zi = b.comp.n_z
            n_ui = b.comp.n_us + b.comp.n_ug
            n_oi = b.comp.n_out_s + b.comp.n_out_g
            u_slice = slice(b.input_off, b.input_off + n_ui)
            if n_zi:
                z_slice = slice(zo, zo + n_zi)
                dg_dz[z_slice, z_slice] = b.comp.Gz(x_i, z_i, u_i)
                dg_du[z_slice, u_slice] = b.comp.Gu(x_i, z_i, u_i)
                dh_dz[b.output_off : b.output_off + n_oi, z_slice] = b.comp.Hz(x_i, z_i, u_i)
            dh_du[b.output_off : b.output_off + n_oi, u_slice] = b.comp.Hu(x_i, z_i, u_i)

        G = self.topology.G
        top = np.hstack([dg_dz, dg_du])
        bottom = np.hstack([-G @ dh_dz, np.eye(n_u) - G @ dh_du])
        return np.vstack([top, bottom])

    def solve_algebraic(
        self, x_all: Vec, u_exo: Vec, z_guess: Vec | None = None, u_guess: Vec | None = None
    ) -> tuple[Vec, Vec]:
        """Newton-solves the coupled ``g(x, z, u) = 0`` / interconnection
        system for ``(z, u)`` at a fixed state ``x`` and exogenous input,
        using the analytic Jacobian from ``_residual_jacobian``.

        Uses ``method="hybr"`` (Powell hybrid, MINPACK's ``hybrj`` when a
        Jacobian is supplied) first — fast, and converges from the naive
        per-component guess for every network validated so far (WSCC-9,
        CIGRE). An earlier version of this method used SciPy's
        numerically-estimated Jacobian, and with it ``"hybr"`` converged
        fine exactly at the operating point but reliably failed ("not
        making good progress") for the nearby points an ODE integrator's
        own step-size/collocation probing evaluates; with the analytic
        Jacobian that mostly went away.

        "Mostly" is doing real work in that sentence: SMIB's initial guess
        (found when wiring the infinite bus in) sits ~0.23 pu away from the
        true coupled equilibrium in the reactive-current direction — WSCC's
        equivalent offset is ~0.08 pu and ``hybr`` handles it; SMIB's is
        outside ``hybr``'s basin of convergence from that starting point,
        even though a true equilibrium demonstrably exists nearby (``lm``
        finds it to machine precision from the identical guess). Rather
        than tune per-network tolerances, this retries with
        Levenberg-Marquardt (``"lm"``, also given the analytic Jacobian) on
        an ``hybr`` failure — slower (estimates its own step directions
        instead of taking Newton steps directly) but has not yet failed
        where ``hybr`` also failed, on any network tried.
        """
        if z_guess is None or u_guess is None:
            z_guess, u_guess = self.initial_algebraic_guess()
        x_list = self._unpack_x(x_all)
        zu0 = np.concatenate([z_guess, u_guess])
        sol = scipy.optimize.root(
            self._residual, zu0, args=(x_list, u_exo), jac=self._residual_jacobian, method="hybr",
        )
        if not sol.success:
            sol = scipy.optimize.root(
                self._residual, zu0, args=(x_list, u_exo), jac=self._residual_jacobian, method="lm",
            )
        if not sol.success:
            raise RuntimeError(f"algebraic Newton solve failed: {sol.message}")
        return sol.x[: self.n_z], sol.x[self.n_z :]

    def rhs(
        self, x_all: Vec, u_exo: Vec, z_guess: Vec | None = None, u_guess: Vec | None = None
    ) -> tuple[Vec, Vec, Vec]:
        """One nonlinear DAE evaluation: solves the algebraic system, then
        returns ``(xdot, z, u)`` — the latter two so callers can warm-start
        the next call's Newton solve from this one's solution.
        """
        z_all, u_all = self.solve_algebraic(x_all, u_exo, z_guess, u_guess)
        x_list, z_list, u_list = self._unpack_x(x_all), self._unpack_z(z_all), self._unpack_u(u_all)
        parts = [
            b.comp.f(x_i, z_i, u_i)
            for b, x_i, z_i, u_i in zip(self.blocks, x_list, z_list, u_list)
            if b.comp.n_states
        ]
        xdot = np.concatenate(parts) if parts else np.zeros(0)
        return xdot, z_all, u_all

    def _inputs_and_outputs(self, x_all: Vec, z_all: Vec, u_all: Vec) -> tuple[dict[str, float], dict[str, float]]:
        """Named inputs/outputs at one already-solved ``(x, z, u)`` triple
        -- the per-sample body shared by :meth:`recover_inputs_and_outputs`
        (called once per already-integrated sample, each needing its own
        fresh :meth:`solve_algebraic` first) and :func:`simulate_steps`
        (called once per solver step, where ``(z, u)`` are already that
        step's own cached algebraic solution -- no extra solve needed
        there, see that function's own docstring).
        """
        x_list, z_list, u_list = self._unpack_x(x_all), self._unpack_z(z_all), self._unpack_u(u_all)
        inputs: dict[str, float] = {}
        outputs: dict[str, float] = {}
        for b, x_b, z_b, u_b in zip(self.blocks, x_list, z_list, u_list):
            for k, name in enumerate(b.comp.input_names):
                inputs[f"{name}_{{{b.name}}}"] = float(u_b[k])
            if b.comp.output_names:
                y_b = b.comp.h(x_b, z_b, u_b)
                for k, name in enumerate(b.comp.output_names):
                    outputs[f"{name}_{{{b.name}}}"] = float(y_b[k])
        return inputs, outputs

    def recover_signals(
        self, t: Vec, x: Vec, u_exo_fn: Callable[[float], Vec] | None = None,
        measure: Callable[[Vec, Vec, Vec], dict[str, float]] | None = None,
    ) -> tuple[dict[str, Vec], dict[str, Vec], dict[str, Vec]]:
        """Like :meth:`recover_inputs_and_outputs`, plus, when ``measure`` is
        given (``measure(x, z, u) -> {name: value}``, see
        :mod:`g2elin_core.timedomain.measurements`), those measurements at
        each sample -- from the same per-sample algebraic solve, so they cost
        nothing extra when inputs/outputs are recovered anyway."""
        u_exo_default = self.default_u_exo()
        z_guess, u_guess = self.initial_algebraic_guess()
        inputs: dict[str, list[float]] = {name: [] for name in self.input_names}
        outputs: dict[str, list[float]] = {name: [] for name in self.output_names}
        meas: dict[str, list[float]] = {}
        for i in range(len(t)):
            x_i = x[:, i]
            u_exo = u_exo_fn(t[i]) if u_exo_fn is not None else u_exo_default
            z_i, u_i = self.solve_algebraic(x_i, u_exo, z_guess, u_guess)
            z_guess, u_guess = z_i, u_i
            step_inputs, step_outputs = self._inputs_and_outputs(x_i, z_i, u_i)
            for name, val in step_inputs.items():
                inputs[name].append(val)
            for name, val in step_outputs.items():
                outputs[name].append(val)
            if measure is not None:
                for name, val in measure(x_i, z_i, u_i).items():
                    meas.setdefault(name, []).append(val)
        return (
            {k: np.array(v) for k, v in inputs.items()},
            {k: np.array(v) for k, v in outputs.items()},
            {k: np.array(v) for k, v in meas.items()},
        )

    def recover_inputs_and_outputs(
        self, t: Vec, x: Vec, u_exo_fn: Callable[[float], Vec] | None = None
    ) -> tuple[dict[str, Vec], dict[str, Vec]]:
        """Named input/output trajectories at each already-integrated
        sample ``(t[i], x[:, i])`` -- e.g. from :func:`simulate`'s
        ``sol.t``/``sol.x``.

        ``solve_ivp`` only returns state trajectories; the algebraic
        variables ``z``/full inputs ``u`` it solves for internally are
        evaluated at the integrator's own adaptive step points, not
        necessarily at the returned (possibly interpolated) sample times.
        So this re-solves the coupled algebraic system once per returned
        sample (warm-started from the previous one, same as during
        integration) to recover a ``(z, u)`` consistent with that exact
        ``(t[i], x[:, i])``, then evaluates each block's own named outputs
        ``h(x, z, u)[:n_out_s]`` there. Costs about as much as the
        integration itself (one more Newton solve per sample) -- bounded by
        the same ``t_final``/timestep limits as ``simulate()``, not
        unbounded.
        """
        inputs, outputs, _ = self.recover_signals(t, x, u_exo_fn)
        return inputs, outputs


def build_nonlinear_network(network: Network, result: PowerFlowResult) -> NonlinearNetworkModel:
    """The nonlinear counterpart of ``pipeline.linearize_network`` — same
    operating point, same topology, nonlinear component models instead of
    linearized ones. Same DER-type/slack-type restrictions apply (see
    ``interconnect/network_assembly.py``).
    """
    if not result.converged:
        raise ValueError("power flow did not converge; can't build a model around an unsolved point")

    op: NetworkOperatingPoint = compute_operating_point(network, result)
    wb_val = 2 * np.pi * network.f_hz

    der_components: dict[int, NonlinearBlockComp] = {}
    for der_id, sm_op in op.sm_ops.items():
        der_components[der_id] = nonlinear_sm_block(sm_op)
    for der_id, gfm_op in op.gfm_ops.items():
        der_components[der_id] = nonlinear_gfm_block(gfm_op)
    for der_id, gfl_op in op.gfl_ops.items():
        der_components[der_id] = nonlinear_gfl_block(gfl_op)
    for der_id, ib_kwargs in op.ib_ops.items():
        der_components[der_id] = nonlinear_ib_block(**ib_kwargs)
    missing = {d.id for d in network.der_units} - der_components.keys()
    if missing:
        raise NotImplementedError(f"unsupported DER unit type(s) for ids {sorted(missing)}")

    # Each bus's own capacitance -- see pipeline.linearize_network, which
    # builds the linear model from the same numbers.
    node_b = node_capacitances(network)
    node_components = {
        bus_id: nonlinear_node_block(wb_val=wb_val, b_pu=node_b[bus_id], wg0=1.0, vgd_g0=vgd, vgq_g0=vgq)
        for bus_id, (vgd, vgq) in op.node_vg.items()
    }
    line_components = [
        nonlinear_line_block(
            wb_val=wb_val, r_pu=ln.r_pu, x_pu=ln.x_pu, wg0=1.0, ild_g0=i0[0], ilq_g0=i0[1],
            vgdj_g0=op.node_vg[ln.from_bus][0], vgqj_g0=op.node_vg[ln.from_bus][1],
            vgdk_g0=op.node_vg[ln.to_bus][0], vgqk_g0=op.node_vg[ln.to_bus][1],
        )
        for ln, i0 in zip(network.lines, op.line_i0)
    ]
    load_components = []
    for idx, load in enumerate(network.loads):
        r_pu, x_pu = op.load_rx[idx]
        vgd, vgq = op.node_vg[load.bus]
        load_components.append(
            nonlinear_load_block(wb_val=wb_val, r_pu=r_pu, x_pu=x_pu, wg0=1.0, vgd_g0=vgd, vgq_g0=vgq)
        )

    # One reference frame per island (components/frame.py), each following
    # the unit that island is referenced to -- an infinite bus turns at its
    # own fixed speed, a machine or grid-forming converter carries its frame
    # with it. All of them start at the operating point's reference angle.
    frame_components = None if network.frame_follows_slack else {
        ref: nonlinear_frame_block(wb_val=wb_val, theta0=op.theta_g_rad, driven=driven)
        for ref, driven in frame_references(network).items()
    }
    # Shunt reactors: an RL branch to zero volts (capacitor banks are in
    # node_b above and add no block).
    shunt_components = {
        idx: nonlinear_line_block(
            wb_val=wb_val, r_pu=rx[0], x_pu=rx[1], wg0=1.0,
            ild_g0=op.shunt_i0[idx][0], ilq_g0=op.shunt_i0[idx][1],
            vgdj_g0=op.node_vg[network.shunts[idx].bus][0], vgqj_g0=op.node_vg[network.shunts[idx].bus][1],
            vgdk_g0=0.0, vgqk_g0=0.0,
        )
        for idx, rx in op.shunt_rx.items()
    }
    blocks, wiring = build_blocks_and_wiring(
        network,
        der_components=der_components,
        node_components=node_components,
        line_components=line_components,
        load_components=load_components,
        frame_components=frame_components,
        shunt_components=shunt_components,
    )
    topology = compute_topology(blocks, wiring)

    z_offsets: list[int] = []
    off = 0
    for b in blocks:
        z_offsets.append(off)
        off += b.comp.n_z
    return NonlinearNetworkModel(
        blocks=blocks, topology=topology, z_offsets=z_offsets, n_z=off, network=network, theta_g0=op.theta_g_rad, op=op,
    )


@dataclass
class EmtSimulationResult:
    t: Vec
    x: Vec  # shape (n_states, len(t))
    state_names: list[str]
    scipy_result: object = field(repr=False)


def find_state_index(model: NonlinearNetworkModel, name_contains: str) -> int:
    """Locates a state by a (unique) substring of its name, e.g.
    ``"dw_r_{SM_2}"`` or just ``"SM_2}"`` if that's unambiguous. Raises if
    zero or more than one state matches, rather than silently guessing.
    """
    matches = [i for i, n in enumerate(model.state_names) if name_contains in n]
    if len(matches) == 0:
        raise ValueError(f"no state name contains {name_contains!r}")
    if len(matches) > 1:
        names = [model.state_names[i] for i in matches]
        raise ValueError(f"{name_contains!r} matches more than one state: {names}")
    return matches[0]


def simulate(
    model: NonlinearNetworkModel,
    t_span: tuple[float, float],
    *,
    x0: Vec | None = None,
    u_exo_fn: Callable[[float], Vec] | None = None,
    t_eval: Vec | None = None,
    method: str = "Radau",
    rtol: float = 1e-6,
    atol: float = 1e-8,
    first_step: float | None = None,
) -> EmtSimulationResult:
    """Integrates the nonlinear DAE from ``x0`` (defaults to the model's own
    operating point — pass an explicit ``x0`` to start from a perturbed
    state). ``u_exo_fn(t)`` lets a caller drive a disturbance (e.g. a P_ref step) by
    returning a modified exogenous-input vector at each t; defaults to
    holding it fixed at the operating point.

    The system's stiffest modes reach ~1e6-1e7 rad/s (see the "SM operating
    point" README note on the Rg penalty parameter), so this defaults to
    ``Radau`` (implicit, L-stable) rather than an explicit method.
    """
    if x0 is None:
        x0 = model.initial_state()
    z_guess, u_guess = model.initial_algebraic_guess()
    u_exo_default = model.default_u_exo()
    cache = {"z": z_guess, "u": u_guess}

    def rhs_fn(t: float, x: Vec) -> Vec:
        u_exo = u_exo_fn(t) if u_exo_fn is not None else u_exo_default
        xdot, z_sol, u_sol = model.rhs(x, u_exo, cache["z"], cache["u"])
        cache["z"], cache["u"] = z_sol, u_sol
        return xdot

    sol = scipy.integrate.solve_ivp(
        rhs_fn, t_span, x0, method=method, t_eval=t_eval, rtol=rtol, atol=atol,
        dense_output=False, first_step=first_step,
    )
    if not sol.success:
        raise RuntimeError(f"EMT integration failed: {sol.message}")
    return EmtSimulationResult(t=sol.t, x=sol.y, state_names=model.state_names, scipy_result=sol)


@dataclass
class EmtStep:
    """One accepted solver step from :func:`simulate_steps` -- the state
    at that step, plus this step's own named inputs/outputs read off the
    same cached algebraic solution the step itself already computed (not
    a fresh :meth:`NonlinearNetworkModel.recover_inputs_and_outputs`
    re-solve -- see that function's docstring for why that one needs to
    re-solve and this one doesn't).
    """

    t: float
    x: Vec
    inputs: dict[str, float]
    outputs: dict[str, float]
    # This step's solved algebraic variables and full input vector (for
    # measurements, which need more than the named inputs/outputs).
    z: Vec | None = None
    u: Vec | None = None


def simulate_steps(
    model: NonlinearNetworkModel,
    t_span: tuple[float, float],
    *,
    x0: Vec | None = None,
    u_exo_fn: Callable[[float], Vec] | None = None,
    method: str = "Radau",
    rtol: float = 1e-6,
    atol: float = 1e-8,
    first_step: float | None = None,
    max_step: float = np.inf,
) -> Iterator[EmtStep]:
    """The step-by-step sibling of :func:`simulate`: yields one
    :class:`EmtStep` after every individual solver step, instead of
    returning the complete trajectory at once -- for a caller that wants
    to *observe* the integration as it happens (a live-tracing UI, a
    progress indicator) rather than wait for the whole thing. Same
    warm-started coupled-algebraic-solve pattern ``simulate()`` uses
    internally, driven through scipy's low-level OOP stepper class (the
    same one ``method``'s string names ``solve_ivp`` accepts, via
    ``getattr(scipy.integrate, method)``) instead of the one-shot
    ``solve_ivp()`` convenience wrapper, which has no way to yield control
    back between steps -- it runs straight through to completion.

    Getting inputs/outputs "for free" at each step (no extra Newton solve,
    unlike :meth:`NonlinearNetworkModel.recover_inputs_and_outputs`) relies
    on evaluating them from ``cache["z"]``/``cache["u"]`` right after each
    ``stepper.step()`` call -- valid because every ``OdeSolver`` subclass
    (Radau included) evaluates the RHS at its own just-accepted solution
    point as part of normal stepping (needed for the next step's error
    estimate), so the cache is already consistent with ``stepper.t``/
    ``stepper.y`` to within the same ``rtol`` the whole integration uses,
    not a separate approximation.

    ``max_step`` caps the solver's own adaptive step size but does not
    force a literal fixed step -- there is no such mode for an implicit
    adaptive method. It also does not reliably control total runtime for a
    stiff system: capping it from unbounded down to 0.03s changed a
    representative stiff trajectory's total step count by under 1% (2036
    vs 2025 steps) while nearly doubling wall-clock time (more, smaller
    steps). Treat it as a *redraw-granularity* knob for a caller watching
    the steps arrive, not a speed control.
    """
    if x0 is None:
        x0 = model.initial_state()
    z_guess, u_guess = model.initial_algebraic_guess()
    u_exo_default = model.default_u_exo()
    cache = {"z": z_guess, "u": u_guess}

    def rhs_fn(t: float, x: Vec) -> Vec:
        u_exo = u_exo_fn(t) if u_exo_fn is not None else u_exo_default
        xdot, z_sol, u_sol = model.rhs(x, u_exo, cache["z"], cache["u"])
        cache["z"], cache["u"] = z_sol, u_sol
        return xdot

    stepper_cls = getattr(scipy.integrate, method)
    stepper = stepper_cls(
        rhs_fn, t_span[0], x0, t_span[1], max_step=max_step, rtol=rtol, atol=atol, first_step=first_step,
    )

    # The stepper's own constructor already called rhs_fn(t_span[0], x0)
    # once (every OdeSolver subclass needs the initial derivative), so
    # cache["z"]/cache["u"] already hold the *solved* (not just guessed)
    # algebraic state for x0 by the time we get here.
    inputs0, outputs0 = model._inputs_and_outputs(x0, cache["z"], cache["u"])
    yield EmtStep(t=t_span[0], x=x0.copy(), inputs=inputs0, outputs=outputs0, z=cache["z"].copy(), u=cache["u"].copy())

    while stepper.status == "running":
        stepper.step()
        inputs_i, outputs_i = model._inputs_and_outputs(stepper.y, cache["z"], cache["u"])
        yield EmtStep(
            t=stepper.t, x=stepper.y.copy(), inputs=inputs_i, outputs=outputs_i, z=cache["z"].copy(), u=cache["u"].copy()
        )

    if stepper.status == "failed":
        raise RuntimeError(f"EMT live integration failed at t={stepper.t}")
