"""Shared machinery for linearized small-signal component models.

Ports the pattern used throughout ``Functions/sym*.m``: build a full
nonlinear DAE symbolically (differential equations ``diffeqVec``, algebraic
constraints ``algeqVec``, outputs ``outputeqVec``), take Jacobians to get
``Fx, Fz, Fu, Gx, Gz, Gu, Hx, Hz, Hu``, then eliminate the algebraic
variables:

    A = Fx - Fz @ inv(Gz) @ Gx      B = Fu - Fz @ inv(Gz) @ Gu
    C = Hx - Hz @ inv(Gz) @ Gx      D = Hu - Hz @ inv(Gz) @ Gu

One deliberate deviation from the MATLAB source: the MATLAB toolbox
substitutes numeric parameter/operating-point values *before* forming
``Ai = Fx - Fz*inv(Gz)*Gx`` (symbolic matrix inversion, feasible in MATLAB
because it's cached per case). Substitution and matrix algebra commute, so
this module instead keeps ``Fx..Hu`` symbolic (cheap — just
differentiation, no inversion) and defers substitution to
:meth:`ComponentDAE.linearize`, where the ``Gz`` inverse and matrix products
run numerically in numpy. Same result, avoids symbolic inversion of a
15x15+ matrix.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field

import numpy as np
import sympy as sp

# The three things that can become of one state when a model's order is
# reduced (see :func:`apply_reduction`).
DYNAMIC = "dynamic"
ALGEBRAIC = "algebraic"
FROZEN = "frozen"


@dataclass(frozen=True)
class LinearComponent:
    """A numeric linearized state-space, split into "own" (s) and "grid" (g)
    inputs/outputs the way the MATLAB toolbox splits ``nIP_*`` / ``nOP_*``.

    ``us`` inputs are the component's own setpoints/references (e.g. a
    generator's ``P_ref``); ``ug`` inputs/outputs are the interconnection
    signals (node voltages, injected currents, the global reference angle)
    that :mod:`g2elin_core.interconnect` wires between components.
    """

    A: np.ndarray
    B: np.ndarray
    C: np.ndarray
    D: np.ndarray
    n_us: int
    n_ug: int
    n_out_s: int
    n_out_g: int
    state_names: list[str]
    input_names: list[str]  # us names only, length n_us
    output_names: list[str]  # s-output names only, length n_out_s

    @property
    def n_states(self) -> int:
        return self.A.shape[0]


@dataclass(frozen=True)
class NonlinearFuncs:
    """Numeric (numpy-backed) callables for a component's full nonlinear DAE,
    lambdified from the same symbolic equations the linearization was built
    from. All three take ``(x, z, u, p)`` — numpy arrays matching
    :attr:`ComponentDAE.state_syms` / ``alg_syms`` / ``input_syms`` /
    ``param_syms`` order — and return a 1-D numpy array.

    This is what makes EMT (nonlinear) time-domain simulation possible
    without re-deriving physics: it's the *same* ``diffeqVec``/``algeqVec``/
    ``outputeqVec`` every ``sym*.m`` file builds before linearizing, kept
    instead of discarded (see the migration plan's "the nonlinear models
    already exist" finding).
    """

    f: object  # (x, z, u, p) -> state derivatives, shape (n_states,)
    g: object  # (x, z, u, p) -> algebraic residuals, shape (n_alg,)
    h: object  # (x, z, u, p) -> outputs, shape (n_out_s + n_out_g,)


@dataclass(frozen=True)
class NonlinearJacobians:
    """Numeric Jacobian callables for ``f``, ``g`` and ``h`` — see
    :meth:`ComponentDAE.nonlinear_jacobians`. Each takes ``(x, z, u, p)``
    and returns a 2-D numpy array.

    ``Fx``/``Fz``/``Fu`` are what an *implicit* ODE solver needs: with them
    the total derivative ``df/dx`` along the constraint manifold can be
    assembled directly (see
    :meth:`~g2elin_core.timedomain.emt.NonlinearNetworkModel.ode_jacobian`)
    instead of being finite-differenced, which costs one full coupled Newton
    solve per state.
    """

    Gz: object  # d(g)/d(z), shape (n_z, n_z)
    Gx: object  # d(g)/d(x), shape (n_z, n_x)
    Gu: object  # d(g)/d(u), shape (n_z, n_us+n_ug)
    Hz: object  # d(h)/d(z), shape (n_out, n_z)
    Hx: object  # d(h)/d(x), shape (n_out, n_x)
    Hu: object  # d(h)/d(u), shape (n_out, n_us+n_ug)
    Fx: object  # d(f)/d(x), shape (n_x, n_x)
    Fz: object  # d(f)/d(z), shape (n_x, n_z)
    Fu: object  # d(f)/d(u), shape (n_x, n_us+n_ug)


@dataclass(frozen=True)
class ComponentDAE:
    """Symbolic Jacobians for one component, not yet substituted with numbers."""

    Fx: sp.Matrix
    Fz: sp.Matrix
    Fu: sp.Matrix
    Gx: sp.Matrix
    Gz: sp.Matrix
    Gu: sp.Matrix
    Hx: sp.Matrix
    Hz: sp.Matrix
    Hu: sp.Matrix
    n_us: int
    n_ug: int
    n_out_s: int
    n_out_g: int
    state_names: list[str] = field(default_factory=list)
    input_names: list[str] = field(default_factory=list)
    output_names: list[str] = field(default_factory=list)
    # Raw nonlinear equations, kept for EMT simulation (see NonlinearFuncs).
    state_syms: list[sp.Symbol] = field(default_factory=list)
    alg_syms: list[sp.Symbol] = field(default_factory=list)
    input_syms: list[sp.Symbol] = field(default_factory=list)
    diffeq_exprs: list[sp.Expr] = field(default_factory=list)
    algeq_exprs: list[sp.Expr] = field(default_factory=list)
    output_exprs: list[sp.Expr] = field(default_factory=list)
    # Lazily-filled store for the lambdified Jacobians (see cached_jacobians).
    # A dict rather than a plain field because this dataclass is frozen:
    # mutating one is allowed, rebinding it is not.
    _cache: dict = field(default_factory=dict, compare=False, repr=False)

    def param_syms(self) -> list[sp.Symbol]:
        """Every free symbol in the equations that isn't a state/alg/input —
        i.e. the physical parameters (``wb``, ``Ra``, ``Lf``, ...).
        Sorted by name for a stable, reproducible ordering.
        """
        known = set(self.state_syms) | set(self.alg_syms) | set(self.input_syms)
        free: set[sp.Symbol] = set()
        for e in (*self.diffeq_exprs, *self.algeq_exprs, *self.output_exprs):
            free |= e.free_symbols
        return sorted(free - known, key=lambda s: s.name)

    def point_from_subs(self, subs: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Pulls ``(x0, z0, u0, p0)`` numeric arrays — in the exact symbol
        order :meth:`nonlinear_funcs`'s callables expect — out of the same
        substitution dict :meth:`linearize` takes. Lets callers reuse one
        dict to get both the linearized *and* the nonlinear operating point.
        """

        def pull(syms: list[sp.Symbol]) -> np.ndarray:
            return np.array([float(subs[s]) for s in syms], dtype=np.float64)

        return pull(self.state_syms), pull(self.alg_syms), pull(self.input_syms), pull(self.param_syms())

    def nonlinear_funcs(self) -> NonlinearFuncs:
        """Lambdify ``diffeqVec``/``algeqVec``/``outputeqVec`` into numeric
        callables. Not cached here (sympy ``Matrix`` isn't reliably
        hashable) — component modules wrap this in their own
        ``@lru_cache``d accessor, the same way ``sm_dae()`` etc. are cached.
        """
        p_syms = self.param_syms()
        args = (self.state_syms, self.alg_syms, self.input_syms, p_syms)

        def _lambdify(exprs: list[sp.Expr]):
            if not exprs:
                empty = np.zeros(0)
                return lambda x, z, u, p: empty
            fn = sp.lambdify(args, exprs, modules="numpy")
            return lambda x, z, u, p: np.asarray(fn(x, z, u, p), dtype=np.float64)

        return NonlinearFuncs(
            f=_lambdify(self.diffeq_exprs),
            g=_lambdify(self.algeq_exprs),
            h=_lambdify(self.output_exprs),
        )

    def nonlinear_jacobians(self) -> NonlinearJacobians:
        """Lambdify ``Gz``/``Gx``/``Gu``/``Hz``/``Hx``/``Hu`` — the same
        Jacobians :meth:`linearize` uses, but as numeric ``(x, z, u, p) ->
        matrix`` callables instead of being substituted with one fixed
        operating point. This is what a Newton solver needs for the
        *nonlinear* coupled algebraic system (see ``timedomain/emt.py``):
        unlike :meth:`linearize`, the Newton solve has to evaluate the
        Jacobian at whatever point the solver is currently at, not just once
        at an operating point.
        """
        p_syms = self.param_syms()
        args = (self.state_syms, self.alg_syms, self.input_syms, p_syms)

        def _lambdify_matrix(mat: sp.Matrix, shape: tuple[int, int]):
            if mat.rows == 0 or mat.cols == 0:
                empty = np.zeros(shape)
                return lambda x, z, u, p: empty
            fn = sp.lambdify(args, mat, modules="numpy")
            return lambda x, z, u, p: np.asarray(fn(x, z, u, p), dtype=np.float64).reshape(shape)

        n_x, n_z, n_out = len(self.state_syms), len(self.alg_syms), self.n_out_s + self.n_out_g
        n_u = self.n_us + self.n_ug
        return NonlinearJacobians(
            Gz=_lambdify_matrix(self.Gz, (n_z, n_z)),
            Gx=_lambdify_matrix(self.Gx, (n_z, n_x)),
            Gu=_lambdify_matrix(self.Gu, (n_z, n_u)),
            Hz=_lambdify_matrix(self.Hz, (n_out, n_z)),
            Hx=_lambdify_matrix(self.Hx, (n_out, n_x)),
            Hu=_lambdify_matrix(self.Hu, (n_out, n_u)),
            Fx=_lambdify_matrix(self.Fx, (n_x, n_x)),
            Fz=_lambdify_matrix(self.Fz, (n_x, n_z)),
            Fu=_lambdify_matrix(self.Fu, (n_x, n_u)),
        )

    def cached_jacobians(self) -> NonlinearJacobians:
        """:meth:`nonlinear_jacobians`, built once per DAE.

        The component modules already cache each *DAE*, so caching the
        lambdified Jacobians on it means the (slow) lambdify happens once per
        model variant for the life of the process.
        """
        if "jacobians" not in self._cache:
            self._cache["jacobians"] = self.nonlinear_jacobians()
        return self._cache["jacobians"]

    def linearize(self, subs: dict) -> LinearComponent:
        """Substitute numeric parameter + operating-point values and eliminate z.

        Evaluated through the lambdified Jacobians rather than by symbolic
        substitution. The two agree to ~1e-22 -- they are the same
        expressions -- but ``Matrix.subs`` on expressions this size is
        expensive enough to dominate everything built on it: measured on one
        synchronous machine, 952 ms against 0.74 ms, and a network's
        linearisation calls this once per component.
        """
        try:
            x0, z0, u0, p0 = self.point_from_subs(subs)
        except KeyError as e:
            raise ValueError(f"no value given for symbol {e.args[0]}") from e

        j = self.cached_jacobians()
        Fx, Fz, Fu = j.Fx(x0, z0, u0, p0), j.Fz(x0, z0, u0, p0), j.Fu(x0, z0, u0, p0)
        Gx, Gz, Gu = j.Gx(x0, z0, u0, p0), j.Gz(x0, z0, u0, p0), j.Gu(x0, z0, u0, p0)
        Hx, Hz, Hu = j.Hx(x0, z0, u0, p0), j.Hz(x0, z0, u0, p0), j.Hu(x0, z0, u0, p0)

        if Gz.shape[0] or Gz.shape[1]:
            try:
                inv_gz = np.linalg.inv(Gz)
            except np.linalg.LinAlgError as e:
                # A singular algebraic Jacobian means this component's own
                # operating point is degenerate (e.g. from an unusual
                # dispatch/parameter combination in a hand-built network) --
                # every preset in this codebase avoids this by construction,
                # but nothing prevents it structurally. Re-raised as a
                # ValueError (-> HTTP 422 at the API layer, matching
                # operating_point.py's own network-content-problem checks)
                # instead of letting a bare LinAlgError surface as an
                # unhandled 500 with no indication of what's wrong.
                raise ValueError(
                    f"a component's algebraic Jacobian is singular at its operating point ({e}) -- "
                    "this usually means a degenerate dispatch or parameter combination (e.g. a DER "
                    "dispatched at ~0 pu, or two units fighting over the same voltage setpoint); check "
                    "this network's DER dispatch (p_set_mw/q_set_mvar/v_set_pu) and try again"
                ) from e
        else:
            inv_gz = Gz  # 0x0, matmul with it is a shape-consistent no-op

        A = Fx - Fz @ inv_gz @ Gx
        B = Fu - Fz @ inv_gz @ Gu
        C = Hx - Hz @ inv_gz @ Gx
        D = Hu - Hz @ inv_gz @ Gu

        return LinearComponent(
            A=A,
            B=B,
            C=C,
            D=D,
            n_us=self.n_us,
            n_ug=self.n_ug,
            n_out_s=self.n_out_s,
            n_out_g=self.n_out_g,
            state_names=list(self.state_names),
            input_names=list(self.input_names),
            output_names=list(self.output_names),
        )


def build_dae(
    *,
    state_vec: list[sp.Symbol],
    alg_vec: list[sp.Symbol],
    input_vec: list[sp.Symbol],
    diffeq_vec: list[sp.Expr],
    algeq_vec: list[sp.Expr],
    output_vec: list[sp.Expr],
    n_us: int,
    n_ug: int,
    n_out_s: int,
    n_out_g: int,
    state_names: list[str],
    input_names: list[str],
    output_names: list[str],
) -> ComponentDAE:
    """Jacobian bookkeeping shared by every ``sym*`` port — mirrors the
    ``Fx = jacobian(diffeqVec, stateVec)`` block repeated in every MATLAB file.
    """
    f = sp.Matrix(diffeq_vec) if diffeq_vec else sp.Matrix(0, 1, [])
    g = sp.Matrix(algeq_vec) if algeq_vec else sp.Matrix(0, 1, [])
    h = sp.Matrix(output_vec) if output_vec else sp.Matrix(0, 1, [])

    x = sp.Matrix(state_vec) if state_vec else sp.Matrix(0, 1, [])
    z = sp.Matrix(alg_vec) if alg_vec else sp.Matrix(0, 1, [])
    u = sp.Matrix(input_vec) if input_vec else sp.Matrix(0, 1, [])

    def jac(expr_mat: sp.Matrix, wrt: sp.Matrix) -> sp.Matrix:
        if expr_mat.rows == 0 or wrt.rows == 0:
            return sp.zeros(expr_mat.rows, wrt.rows)
        return expr_mat.jacobian(wrt)

    return ComponentDAE(
        Fx=jac(f, x),
        Fz=jac(f, z),
        Fu=jac(f, u),
        Gx=jac(g, x),
        Gz=jac(g, z),
        Gu=jac(g, u),
        Hx=jac(h, x),
        Hz=jac(h, z),
        Hu=jac(h, u),
        n_us=n_us,
        n_ug=n_ug,
        n_out_s=n_out_s,
        n_out_g=n_out_g,
        state_names=state_names,
        input_names=input_names,
        output_names=output_names,
        state_syms=list(state_vec),
        alg_syms=list(alg_vec),
        input_syms=list(input_vec),
        diffeq_exprs=list(diffeq_vec),
        algeq_exprs=list(algeq_vec),
        output_exprs=list(output_vec),
    )


def _to_numpy(mat: sp.Matrix, subs: dict) -> np.ndarray:
    if mat.rows == 0 or mat.cols == 0:
        return np.zeros((mat.rows, mat.cols))
    substituted = mat.subs(subs)
    free = substituted.free_symbols
    if free:
        raise ValueError(f"unsubstituted symbols remain: {sorted(str(s) for s in free)}")
    return np.array(substituted.evalf(), dtype=np.float64).reshape(mat.rows, mat.cols)


def eq_symbol(sym: sp.Symbol) -> sp.Symbol:
    """The ``<name>_0`` equilibrium-point symbol for ``sym``, matching the
    MATLAB toolbox's ``sym(strcat(char(s),'_0'))`` convention.
    """
    return sp.Symbol(f"{sym.name}_0")


def mode_key(modes: Mapping[str, str] | None) -> tuple[tuple[str, str], ...]:
    """A hashable, canonical form of a mode map, for the ``lru_cache`` on
    each component's ``*_dae()`` builder.

    Dynamic entries are dropped, so a fully-dynamic model hashes to ``()``
    and shares its cache entry with a call that passed no modes at all --
    the full-order DAE is built once however it was asked for.
    """
    return tuple(sorted((k, v) for k, v in (modes or {}).items() if v != DYNAMIC))


def apply_reduction(
    *,
    states: list[tuple[sp.Symbol, sp.Expr, str]],
    modes: Mapping[str, str],
    alg_vec: list[sp.Symbol],
    algeq_vec: list[sp.Expr],
    output_vec: list[sp.Expr],
) -> tuple[list[sp.Symbol], list[sp.Expr], list[str], list[sp.Symbol], list[sp.Expr], list[sp.Expr]]:
    """Lower a component's model order by reclassifying individual states.

    ``states`` is the component's full-order ``(symbol, d/dt expression,
    display name)`` list in its canonical order; ``modes`` maps a state
    symbol's *name* to one of :data:`DYNAMIC`, :data:`ALGEBRAIC` or
    :data:`FROZEN` (anything unnamed stays dynamic). Returns the six lists
    :func:`build_dae` takes, ready to pass straight through.

    The two reductions are mathematically different and are *not*
    interchangeable, which is why they're separate modes rather than one
    "drop this state" flag:

    - :data:`ALGEBRAIC` (residualization / singular perturbation) replaces
      ``dx/dt = expr`` with ``0 = expr`` and moves ``x`` into the algebraic
      vector. The state still takes whatever value the rest of the system
      demands of it, instantaneously. Valid for states much *faster* than
      the phenomena being studied -- stator flux, filter currents, network
      branches, inner control loops.
    - :data:`FROZEN` (truncation) substitutes ``x`` by its equilibrium
      value ``x_0`` everywhere and drops its equation entirely. The state
      simply doesn't move. Valid for states much *slower* than the window
      of interest -- the field flux in the classical machine model, where
      ``E'`` is held constant.

    Applying the fast-state treatment to a slow state (or the reverse) is
    the most common way a hand-rolled reduced model goes quietly wrong: a
    field flux with ``dpsi_fd/dt = 0`` is an *instantaneous* field winding,
    the exact opposite of the classical model's constant one.

    A frozen state's ``x_0`` becomes an ordinary free symbol, so
    :meth:`ComponentDAE.param_syms` picks it up as a parameter and
    :meth:`ComponentDAE.point_from_subs` reads its value from the same
    operating-point substitution dict everything else comes from -- which
    is why reduction needs no new initialization code anywhere.
    """
    frozen_subs: dict[sp.Symbol, sp.Symbol] = {
        sym: eq_symbol(sym) for sym, _, _ in states if modes.get(sym.name, DYNAMIC) == FROZEN
    }

    def sub(expr: sp.Expr) -> sp.Expr:
        return expr.subs(frozen_subs) if frozen_subs else expr

    kept_syms: list[sp.Symbol] = []
    kept_eqs: list[sp.Expr] = []
    kept_names: list[str] = []
    new_alg_syms: list[sp.Symbol] = []
    new_alg_eqs: list[sp.Expr] = []

    for sym, eq, name in states:
        mode = modes.get(sym.name, DYNAMIC)
        if mode == FROZEN:
            continue  # equation dropped; the symbol is now the x_0 parameter
        if mode == ALGEBRAIC:
            new_alg_syms.append(sym)
            new_alg_eqs.append(sub(eq))  # dx/dt = eq  ->  0 = eq
            continue
        if mode != DYNAMIC:
            raise ValueError(f"unknown state mode {mode!r} for {sym.name!r}")
        kept_syms.append(sym)
        kept_eqs.append(sub(eq))
        kept_names.append(name)

    # Reclassified states go after the component's own algebraic variables,
    # in their original state order, so the z layout is deterministic.
    return (
        kept_syms,
        kept_eqs,
        kept_names,
        list(alg_vec) + new_alg_syms,
        [sub(e) for e in algeq_vec] + new_alg_eqs,
        [sub(e) for e in output_vec],
    )


def equilibrium_subs(subs: dict) -> dict:
    """``subs`` plus an ``x_0`` entry for every symbol in it -- the values
    :data:`FROZEN` states are substituted by.

    Components call this on their own operating-point dict so that freezing
    any state works without the dict having to know which ones were frozen.
    The extra entries are harmless when nothing is frozen: ``param_syms``
    only returns symbols that actually appear in the equations.
    """
    return {**subs, **{eq_symbol(sym): value for sym, value in subs.items()}}
