"""Numerical (Cowell) propagation under point-mass + J2 gravity.

    a = -mu r/|r|^3 - 3/2 J2 mu Re^2 / |r|^5 * [x(1-5z^2/r^2), y(1-5z^2/r^2), z(3-5z^2/r^2)]

integrated with DOP853 (8th-order Dormand-Prince, dense output) at rtol 1e-12
in the TEOE frame (Z = true pole of epoch). Includes the short-period J2
terms that the analytic secular model omits; used for sensitivity analysis
and for the nonlinear validation of the Part II relative-motion design.

``fit_osculating_to_mean`` finds the osculating initial state whose orbit-
averaged ("mean") elements equal prescribed values -- the numerically
consistent way to honour "550 km circular" as MEAN elements.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from .elements import cart2kep, kep2cart


def j2_accel(r, mu, re, j2):
    x, y, z = r[..., 0], r[..., 1], r[..., 2]
    r2 = x * x + y * y + z * z
    rn = np.sqrt(r2)
    k = 1.5 * j2 * mu * re * re / (r2 * r2 * rn)
    z2 = 5.0 * z * z / r2
    a = -mu * r / (rn ** 3)[..., None]
    a = a + np.stack([-k * x * (1 - z2), -k * y * (1 - z2), -k * z * (3 - z2)], -1)
    return a


def propagate(r0, v0, t_span, mu, re, j2, t_eval=None, rtol=1e-12, atol=1e-10, dense=True):
    def rhs(_t, y):
        return np.concatenate([y[3:], j2_accel(y[:3], mu, re, j2)])

    return solve_ivp(rhs, t_span, np.concatenate([r0, v0]), method="DOP853", rtol=rtol, atol=atol,
                     dense_output=dense, t_eval=t_eval)


class CowellJ2Orbit:
    def __init__(self, r0, v0, mu, re, j2, t_end, t_start=0.0, name="Cowell J2 (numerical)"):
        self.mu, self.re, self.j2, self.name = mu, re, j2, name
        self._sol = propagate(np.asarray(r0, float), np.asarray(v0, float), (t_start, t_end), mu, re, j2)
        if not self._sol.success:
            raise RuntimeError(self._sol.message)

    def state_teoe(self, t):
        y = self._sol.sol(np.asarray(t, dtype=float))
        y = np.moveaxis(np.atleast_2d(y), 0, -1)
        return y[..., :3], y[..., 3:]


def mean_elements_numerical(r0, v0, mu, re, j2, n_orbits=3, samples=2000):
    """Orbit-averaged elements over ``n_orbits`` (secular trend of the angles
    removed by linear regression; intercept = mean value at t=0)."""
    a0 = 1.0 / (2.0 / np.linalg.norm(r0) - np.dot(v0, v0) / mu)
    T = 2 * np.pi * np.sqrt(a0 ** 3 / mu) * n_orbits
    t = np.linspace(0.0, T, samples * n_orbits)
    sol = propagate(r0, v0, (0.0, T), mu, re, j2, t_eval=t, dense=False)
    r, v = sol.y[:3].T, sol.y[3:].T
    a, e, i, raan, argp, nu, u = cart2kep(r, v, mu)
    ex = e * np.cos(argp)
    ey = e * np.sin(argp)
    uu = np.unwrap(u)
    rr = np.unwrap(raan)
    pu = np.polyfit(t, uu, 1)
    pr = np.polyfit(t, rr, 1)
    return {"a": a.mean(), "ex": ex.mean(), "ey": ey.mean(), "i": i.mean(),
            "raan0": pr[1], "raan_rate": pr[0], "u0": pu[1], "u_rate": pu[0]}


def fit_osculating_to_mean(a_m, i_m, raan_m, u_m, mu, re, j2, iters=12, tol=None):
    """Osculating (r0, v0) whose numerical mean elements equal
    (a_m, e_m = 0, i_m, raan_m, u_m).

    Fixed-point (Picard) iteration p <- p + (target - mean(p)) in the
    non-singular space (a, ex, ey, i, raan, u). The mean->osculating map is the
    identity plus O(J2) short-period terms, so the iteration contracts with a
    ratio of order J2 and converges in a few steps (robust where Newton with a
    finite-difference Jacobian is ill-conditioned near e = 0)."""
    def elems_to_state(p):
        a, ex, ey, i, raan, u = p
        e = np.hypot(ex, ey)
        w = np.arctan2(ey, ex) if e > 0 else 0.0
        nu = u - w
        r, v = kep2cart(a, e, i, raan, w, nu, mu)
        return r, v

    def mean_of(p):
        r, v = elems_to_state(p)
        m = mean_elements_numerical(r, v, mu, re, j2, n_orbits=2, samples=1500)
        return np.array([m["a"], m["ex"], m["ey"], m["i"], m["raan0"], m["u0"]])

    target = np.array([a_m, 0.0, 0.0, i_m, raan_m, u_m])
    tol = tol if tol is not None else np.array([1e-6, 1e-9, 1e-9, 1e-10, 1e-10, 1e-10])
    p = target.copy()
    for _ in range(iters):
        res = target - mean_of(p)
        res[4:] = (res[4:] + np.pi) % (2 * np.pi) - np.pi
        p = p + res
        if np.all(np.abs(res) < tol):
            break
    return elems_to_state(p)
