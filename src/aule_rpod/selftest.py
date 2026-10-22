"""Fast built-in self-test (< 30 s): every check returns a result code, so a run
can prove its own health before producing products (exit code != 0 on failure).
"""
from __future__ import annotations

import time
from enum import Enum

import numpy as np


class Result(Enum):
    PASSED = "PASSED"
    FAILED = "FAILED"


def _quat_convention():
    from .attitude.quaternion import attitude_matrix, body_to_reference
    q = np.array([0, 0, np.sin(np.pi / 4), np.cos(np.pi / 4)])
    ok = np.allclose(body_to_reference(q, [1, 0, 0]), [0, 1, 0]) and np.allclose(attitude_matrix(q) @ [0, 1, 0], [1, 0, 0])
    return ok, "q=(0,0,s45,c45): body X -> ref +Y; A(q) maps ref +Y -> body X"


def _attitude_file(cfg):
    from .attitude import read_stk_attitude
    p = read_stk_attitude(cfg.attitude_path)
    ok = p.meta["parsed_points"] == p.meta["declared_points"] == 14401 and p.meta["max_norm_error"] < 1e-12 \
        and p.coordinate_axes.upper().startswith("LVLH") and np.allclose(np.diff(p.times_s), 60.0)
    return ok, f"{p.meta['parsed_points']} pts, axes '{p.coordinate_axes}', sha256 {p.meta['sha256'][:12]}..."


def _attitude_law(cfg):
    from .attitude import read_stk_attitude
    from .attitude.identification import identify_precession_law
    law = identify_precession_law(read_stk_attitude(cfg.attitude_path))
    ok = law.max_residual_rad < 1e-8 and abs(abs(np.degrees(law.rate_p_rad_s)) - 1.2) < 1e-9
    return ok, f"precession {np.degrees(law.rate_p_rad_s):.9f} deg/s, nutation {law.nutation_deg:.6f} deg, max resid {law.max_residual_rad:.1e} rad"


def _sso_rate(cfg):
    from .orbit.sso import sso_inclination_rad
    from .orbit.secular_j2 import SecularJ2Orbit
    e = cfg.earth
    a = e.equatorial_radius_km + cfg.target_orbit.altitude_km
    i = sso_inclination_rad(a, 0.0, e.mu_km3_s2, e.equatorial_radius_km, e.j2, e.tropical_year_days)
    o = SecularJ2Orbit(a, 0.0, i, 0.0, 0.0, 0.0, e.mu_km3_s2, e.equatorial_radius_km, e.j2)
    rate = np.degrees(o.rates[1]) * 86400
    ok = abs(rate - 360 / e.tropical_year_days) < 1e-9
    return ok, f"i = {np.degrees(i):.5f} deg, dOmega/dt = {rate:.8f} deg/day"


def _cw_stm():
    from scipy.integrate import solve_ivp
    from .rpod import cw
    n = 1.0948e-3
    X0 = np.array([-80, 120, 20, 0.1, -0.2, 0.05])
    sol = solve_ivp(lambda t, X: np.r_[X[3:], cw.natural_accel(n, X[:3], X[3:])], (0, 3000), X0, rtol=1e-12, atol=1e-12)
    err = np.abs(cw.propagate(n, X0, 3000) - sol.y[:, -1]).max()
    f = cw.forced_line(n, [0, 10, 0], [0, 0, 0], 0.06)
    ok = err < 1e-6 and abs(f["dv_continuous"] - 2 * n * 10) < 1e-9
    return ok, f"STM vs ODE {err:.1e} m; forced V-bar dV = 2nL to {abs(f['dv_continuous'] - 2 * n * 10):.0e}"


def _two_impulse():
    from .rpod import cw
    n = 1.0948e-3
    dv1, dv2, v0p, vfm = cw.two_impulse(n, [-80, 120, 20], [0, 0, 0], [0, 25, 0], 235.0)
    X = cw.propagate(n, np.r_[-80, 120, 20, v0p], 235.0)
    ok = np.allclose(X[:3], [0, 25, 0], atol=1e-9) and np.allclose(X[3:], vfm, atol=1e-12)
    return ok, f"arrival error {np.abs(X[:3] - [0, 25, 0]).max():.1e} m"


def _event_finder():
    from .events.finder import find_positive_intervals
    g = lambda t: np.sin(2 * np.pi * t / 100.0) - 0.5      # noqa: E731  roots at 100(1/12 + k), 100(5/12 + k)
    iv = find_positive_intervals(g, 0.0, 1000.0, 7.0, tol=1e-6)
    exact = np.array([[100 * (1 / 12 + k), 100 * (5 / 12 + k)] for k in range(10)])
    ok = iv.shape == exact.shape and np.abs(iv - exact).max() < 1e-5
    # hidden double crossing inside one step
    h = lambda t: 0.001 - (t - 50.3) ** 2                   # noqa: E731  positive only on [50.268, 50.332]
    iv2 = find_positive_intervals(h, 0.0, 100.0, 2.0, tol=1e-7)
    ok = ok and len(iv2) == 1 and abs(iv2[0, 1] - iv2[0, 0] - 2 * np.sqrt(0.001)) < 1e-5
    return ok, f"10 intervals to {np.abs(iv - exact).max():.1e} s; 63-ms blip inside a 2-s step found"


def _intervals():
    from .events import intervals as I
    a = [[0, 10], [20, 30]]
    b = [[5, 25]]
    ok = np.allclose(I.intersect(a, b), [[5, 10], [20, 25]]) and np.allclose(I.complement(a, 0, 40), [[10, 20], [30, 40]]) \
        and np.allclose(I.union(a, b), [[0, 30]])
    return ok, "intersect / union / complement"


def _eclipse():
    from .geometry.eclipse import g_full_sun, g_umbra
    s = np.array([[1.496e8, 0, 0]])
    behind = np.array([[-7000.0, 0, 0]])
    lit = np.array([[7000.0, 0, 0]])
    ok = g_umbra(behind, s, 6378.137)[0] > 0 and g_full_sun(lit, s, 6378.137)[0] > 0 and g_full_sun(behind, s, 6378.137)[0] < 0
    return ok, "anti-Sun point in umbra, sub-Sun point lit"


CHECKS = [("QUATERNION_CONVENTION", lambda c: _quat_convention()), ("ATTITUDE_FILE_INTEGRITY", _attitude_file),
          ("ATTITUDE_LAW_IDENTIFIED", _attitude_law), ("SSO_NODAL_RATE", _sso_rate),
          ("CW_STM_AND_FORCED_MOTION", lambda c: _cw_stm()), ("CW_TWO_IMPULSE_TARGETING", lambda c: _two_impulse()),
          ("EVENT_FINDER", lambda c: _event_finder()), ("INTERVAL_ALGEBRA", lambda c: _intervals()),
          ("ECLIPSE_GEOMETRY", lambda c: _eclipse())]


def run(cfg, verbose=True):
    t0 = time.time()
    rows = []
    for code, fn in CHECKS:
        try:
            ok, detail = fn(cfg)
        except Exception as exc:  # pragma: no cover
            ok, detail = False, repr(exc)
        rows.append((code, Result.PASSED if ok else Result.FAILED, detail))
    overall = Result.PASSED if all(r[1] is Result.PASSED for r in rows) else Result.FAILED
    if verbose:
        w = max(len(r[0]) for r in rows)
        for code, res, detail in rows:
            print(f"  {code:<{w}}  {res.value:<6}  {detail}")
        print(f"SELFTEST {overall.value} ({len(rows)} checks, {time.time() - t0:.1f} s)")
    return overall, rows
