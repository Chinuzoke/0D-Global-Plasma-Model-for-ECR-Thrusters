"""
solver.py
---------
Time-integration and steady-state solving utilities for the 0D global
models. Mirrors the approach in Wachs & Jorns (AIAA-2021-3382), who
integrate the transient ODE system forward in time to a quasi-steady-state
("time resolved global model", their Fig. 1) rather than root-finding the
algebraic steady state directly -- this is important because the ODE
system can have multiple algebraic roots (including the trivial
unignited n_e=0 state), and time-marching from a sensible initial guess
naturally selects the physically relevant ignited branch.

We also provide a root-polishing step (scipy.optimize.fsolve) to snap the
late-time transient result onto the exact steady state to floating-point
tolerance, and a convenience `is_steady` check.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve


def integrate_transient(rhs, y0, p, t_span=(0.0, 5e-3), max_step=None, rtol=1e-8, atol=None):
    """
    Integrate dy/dt = rhs(t, y, p) from y0 over t_span using an implicit
    stiff solver (Radau) -- the ion/neutral/energy balance ODEs are
    numerically stiff (electron energy content equilibrates on a much
    faster timescale than neutral gas depletion), exactly as noted in the
    global-model literature (Lieberman & Lichtenberg Ch. 10).

    Returns the scipy `OdeResult` object (has .t, .y, .success, ...).
    """
    if atol is None:
        # absolute tolerances scaled to typical magnitudes: ne,ng ~ 1e16-1e19 m^-3, Te ~ O(1-10) eV
        atol = [1e6, 1e6, 1e-6]

    def _rhs(t, y):
        return rhs(t, y, p)

    sol = solve_ivp(
        _rhs, t_span, y0, method="Radau",
        max_step=max_step if max_step is not None else (t_span[1] - t_span[0]) / 200.0,
        rtol=rtol, atol=atol, dense_output=True,
    )
    return sol


def polish_steady_state(rhs, y_guess, p, xtol=1e-12):
    """
    Root-find the exact steady state (dy/dt = 0) near y_guess using
    scipy.optimize.fsolve, evaluating rhs at a large fixed time t (so that
    any time-dependent power_profile(t) has settled to its steady/CW value).
    """
    t_ss = 1e9

    def _f(y):
        return rhs(t_ss, y, p)

    sol, info, ier, msg = fsolve(_f, y_guess, xtol=xtol, full_output=True)
    converged = (ier == 1)
    return sol, converged, msg


def run_to_steady_state(rhs, y0, p, t_span=(0.0, 5e-3), verbose=True):
    """
    Convenience wrapper: integrate the transient, then polish the late-time
    result to the exact steady state. Returns (sol, y_ss, converged).
    """
    sol = integrate_transient(rhs, y0, p, t_span=t_span)
    if not sol.success:
        raise RuntimeError(f"Transient integration failed: {sol.message}")

    y_late = sol.y[:, -1]
    y_ss, converged, msg = polish_steady_state(rhs, y_late, p)

    if verbose:
        print(f"  transient reached: n_e={y_late[0]:.3e} m^-3, "
              f"n_g={y_late[1]:.3e} m^-3, Te={y_late[2]:.3f} eV")
        print(f"  steady-state polish {'converged' if converged else 'DID NOT converge'} "
              f"({msg.strip()})")
        print(f"  steady state:      n_e={y_ss[0]:.3e} m^-3, "
              f"n_g={y_ss[1]:.3e} m^-3, Te={y_ss[2]:.3f} eV")

    return sol, y_ss, converged
