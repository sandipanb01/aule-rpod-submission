"""Clohessy-Wiltshire (Hill) relative dynamics about a circular reference orbit.

Equations of motion (x radial, y along-track, z normal; Clohessy & Wiltshire 1960):

    xdd = 3 n^2 x + 2 n yd + ux
    ydd =        - 2 n xd + uy
    zdd = - n^2 z          + uz

State-transition matrix Phi(t) = [[Phi_rr, Phi_rv], [Phi_vr, Phi_vv]], c = cos nt, s = sin nt:

    Phi_rr = [[4-3c, 0, 0], [6(s-nt), 1, 0], [0, 0, c]]
    Phi_rv = [[s/n, 2(1-c)/n, 0], [-2(1-c)/n, (4s-3nt)/n, 0], [0, 0, s/n]]
    Phi_vr = [[3ns, 0, 0], [-6n(1-c), 0, 0], [0, 0, -ns]]
    Phi_vv = [[c, 2s, 0], [-2s, 4c-3, 0], [0, 0, c]]

Two-impulse rendezvous (fixed time of flight T; Vallado 2013 Sec. 6.8):
    v0+ = Phi_rv(T)^-1 [rf - Phi_rr(T) r0]          dv1 = v0+ - v0-
    vf- = Phi_vr(T) r0 + Phi_vv(T) v0+               dv2 = vf+ - vf-
Forced straight-line motion r(t) = r0 + w t needs u(t) = -f(r, w) with
f = [3n^2 x + 2n wy, -2n wx, -n^2 z]; for a V-bar translation of length L this
integrates to the classic dV = 2 n L (Fehse 2003, Sec. 3.3).
"""
from __future__ import annotations

import numpy as np


def stm_blocks(n: float, t):
    t = np.asarray(t, dtype=float)
    c, s = np.cos(n * t), np.sin(n * t)
    z, o = np.zeros_like(t), np.ones_like(t)
    rr = np.stack([np.stack([4 - 3 * c, z, z], -1), np.stack([6 * (s - n * t), o, z], -1), np.stack([z, z, c], -1)], -2)
    rv = np.stack([np.stack([s / n, 2 * (1 - c) / n, z], -1), np.stack([-2 * (1 - c) / n, (4 * s - 3 * n * t) / n, z], -1),
                   np.stack([z, z, s / n], -1)], -2)
    vr = np.stack([np.stack([3 * n * s, z, z], -1), np.stack([-6 * n * (1 - c), z, z], -1), np.stack([z, z, -n * s], -1)], -2)
    vv = np.stack([np.stack([c, 2 * s, z], -1), np.stack([-2 * s, 4 * c - 3, z], -1), np.stack([z, z, c], -1)], -2)
    return rr, rv, vr, vv


def stm(n: float, t):
    rr, rv, vr, vv = stm_blocks(n, t)
    top = np.concatenate([rr, rv], -1)
    bot = np.concatenate([vr, vv], -1)
    return np.concatenate([top, bot], -2)


def propagate(n: float, X0, t):
    """Free motion: X(t) for an array of times (t measured from X0's epoch)."""
    return np.einsum("...ij,j->...i", stm(n, t), np.asarray(X0, dtype=float))


def natural_accel(n: float, r, v):
    r, v = np.asarray(r, float), np.asarray(v, float)
    return np.stack([3 * n * n * r[..., 0] + 2 * n * v[..., 1], -2 * n * v[..., 0], -n * n * r[..., 2]], -1)


def two_impulse(n: float, r0, v0_minus, rf, T: float, vf_plus=(0.0, 0.0, 0.0)):
    rr, rv, vr, vv = stm_blocks(n, np.array(T))
    r0, rf = np.asarray(r0, float), np.asarray(rf, float)
    v0p = np.linalg.solve(rv, rf - rr @ r0)
    vfm = vr @ r0 + vv @ v0p
    dv1 = v0p - np.asarray(v0_minus, float)
    dv2 = np.asarray(vf_plus, float) - vfm
    return dv1, dv2, v0p, vfm


def forced_line(n: float, r0, r1, speed: float, dt: float = 0.5):
    """Constant-speed straight-line translation r0 -> r1. Returns times, positions,
    the required control acceleration history and its dV integral."""
    r0, r1 = np.asarray(r0, float), np.asarray(r1, float)
    L = np.linalg.norm(r1 - r0)
    w = (r1 - r0) / L * speed
    T = L / speed
    t = np.linspace(0.0, T, max(int(T / dt), 2) + 1)
    r = r0 + np.outer(t, w)
    u = -natural_accel(n, r, np.broadcast_to(w, r.shape))
    un = np.linalg.norm(u, axis=1)
    dv = float(np.trapezoid(un, t))
    return {"t": t, "r": r, "v": np.broadcast_to(w, r.shape).copy(), "u": u, "dv_continuous": dv, "T": T, "L": L, "w": w}


def min_range_free_drift(n: float, X0, horizon_s: float, dt: float = 1.0):
    t = np.arange(0.0, horizon_s + dt, dt)
    X = propagate(n, X0, t)
    rng = np.linalg.norm(X[:, :3], axis=1)
    k = int(np.argmin(rng))
    return float(rng[k]), float(t[k]), X
