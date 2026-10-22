"""Nonlinear validation of the CW design: both spacecraft propagated numerically
under point-mass + J2 gravity (Cowell, DOP853, rtol 1e-12) in the TEOE frame,
with the planned impulses applied in the target's instantaneous LVLH frame and
the final translation flown with the planned continuous control acceleration.

LVLH <-> inertial mapping (rotating frame, omega = (r x v)/|r|^2):
    r_c = r_t + C^T rho
    v_c = v_t + C^T rho_dot + omega x (C^T rho)
where C has rows (x, y, z) LVLH unit vectors in inertial coordinates.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from ..frames import lvlh_matrix
from ..orbit.cowell import j2_accel


def _rel(rt, vt, rc, vc):
    C = lvlh_matrix(rt, vt)[0]
    w = np.cross(rt, vt) / np.dot(rt, rt)
    dr = rc - rt
    dv = vc - vt - np.cross(w, dr)
    return C @ dr, C @ dv


def _to_inertial(rt, vt, rho, rhod):
    C = lvlh_matrix(rt, vt)[0]
    w = np.cross(rt, vt) / np.dot(rt, rt)
    d = C.T @ rho
    return rt + d, vt + C.T @ rhod + np.cross(w, d)


def simulate(plan, rt0_km, vt0_kms, mu, re, j2, X0_rel):
    """Fly the plan in the nonlinear J2 model. Returns relative-state history and
    the deviation from the CW-planned positions at each waypoint (metres)."""
    km = 1e-3
    rt, vt = np.asarray(rt0_km, float), np.asarray(vt0_kms, float)
    rc, vc = _to_inertial(rt, vt, np.asarray(X0_rel[:3]) * km, np.asarray(X0_rel[3:]) * km)

    def rhs(_t, y, u_lvlh=None):
        r_t, v_t, r_c, v_c = y[:3], y[3:6], y[6:9], y[9:]
        a_t = j2_accel(r_t, mu, re, j2)
        a_c = j2_accel(r_c, mu, re, j2)
        if u_lvlh is not None:
            C = lvlh_matrix(r_t, v_t)[0]
            a_c = a_c + C.T @ (u_lvlh(_t, r_t, v_t, r_c, v_c) * km)
        return np.concatenate([v_t, a_t, v_c, a_c])

    y = np.concatenate([rt, vt, rc, vc])
    t_now = 0.0
    log = []
    burns = {b.id: b for b in plan.burns}
    p = plan.params
    n_ref = np.linalg.norm(np.cross(rt, vt)) / np.dot(rt, rt)

    def apply(bid):
        nonlocal y
        b = burns[bid]
        r_t, v_t = y[:3], y[3:6]
        C = lvlh_matrix(r_t, v_t)[0]
        y[9:] = y[9:] + C.T @ (b.dv * km)

    def coast(T, u=None):
        nonlocal y, t_now
        sol = solve_ivp(lambda t, yy: rhs(t, yy, u), (t_now, t_now + T), y, method="DOP853", rtol=1e-12, atol=1e-12,
                        dense_output=True)
        ts = np.linspace(t_now, t_now + T, max(int(T), 2))
        for tk in ts:
            yk = sol.sol(tk)
            rho, rhod = _rel(yk[:3], yk[3:6], yk[6:9], yk[9:])
            log.append((tk, *(rho / km), *(rhod / km)))
        y = sol.y[:, -1]
        t_now += T

    def rel_now():
        rho, rhod = _rel(y[:3], y[3:6], y[6:9], y[9:])
        return rho / km, rhod / km

    checks = {}
    apply("B1")
    coast(p["T1"])
    checks["WP1"] = np.linalg.norm(rel_now()[0] - p["WP1"])
    apply("B2")
    coast(p["T2"])
    checks["HP"] = np.linalg.norm(rel_now()[0] - p["HP"])
    apply("B3")
    checks["HP_residual_velocity_mm_s"] = 1000 * np.linalg.norm(rel_now()[1])
    coast(p["T_hold"])
    checks["hold_drift_m"] = np.linalg.norm(rel_now()[0] - p["HP"])
    apply("B4")
    # final translation: planned control = - natural CW acceleration along the line (feed-forward)
    a_ax = plan.strategy.axis
    w = -a_ax * p["v_final"]

    def u_ff(_t, r_t, v_t, r_c, v_c):
        rho, _ = _rel(r_t, v_t, r_c, v_c)
        rho = rho / km
        return -np.array([3 * n_ref ** 2 * rho[0] + 2 * n_ref * w[1], -2 * n_ref * w[0], -n_ref ** 2 * rho[2]])

    T_fin = p["rho_hold"] / p["v_final"] if p["d_brake"] == 0 else (p["rho_hold"] - p["d_brake"]) / p["v_final"]
    coast(T_fin, u_ff)
    rho, rhod = rel_now()
    checks["contact_miss_m"] = float(np.linalg.norm(rho - (p["d_brake"] * a_ax if p["d_brake"] else 0.0)))
    checks["contact_velocity_error_mm_s"] = float(1000 * np.linalg.norm(rhod - w))
    return np.array(log), {k: float(v) for k, v in checks.items()}
