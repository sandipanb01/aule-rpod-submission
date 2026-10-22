"""Classical orbital elements <-> Cartesian state (two-body relations)."""
from __future__ import annotations

import numpy as np


def kepler_E(M, e, tol=1e-14, maxit=30):
    """Solve Kepler's equation M = E - e sin E (vectorised Newton)."""
    M = np.asarray(M, dtype=float)
    E = M + e * np.sin(M) if e > 0 else M.copy()
    if e == 0:
        return E
    for _ in range(maxit):
        dE = (E - e * np.sin(E) - M) / (1.0 - e * np.cos(E))
        E = E - dE
        if np.max(np.abs(dE)) < tol:
            break
    return E


def true_from_mean(M, e):
    E = kepler_E(M, e)
    return 2.0 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))


def mean_from_true(nu, e):
    E = 2.0 * np.arctan2(np.sqrt(1 - e) * np.sin(nu / 2), np.sqrt(1 + e) * np.cos(nu / 2))
    return E - e * np.sin(E)


def kep2cart(a, e, i, raan, argp, nu, mu):
    """Vectorised elements -> (r, v) [km, km/s] in the frame the angles refer to."""
    a, e, i, raan, argp, nu = np.broadcast_arrays(*(np.asarray(x, dtype=float) for x in (a, e, i, raan, argp, nu)))
    p = a * (1.0 - e ** 2)
    r_mag = p / (1.0 + e * np.cos(nu))
    u = argp + nu
    cO, sO = np.cos(raan), np.sin(raan)
    ci, si = np.cos(i), np.sin(i)
    cu, su = np.cos(u), np.sin(u)
    # unit vectors: radial and transverse
    rhat = np.stack([cO * cu - sO * su * ci, sO * cu + cO * su * ci, su * si], -1)
    that = np.stack([-cO * su - sO * cu * ci, -sO * su + cO * cu * ci, cu * si], -1)
    sqrt_mu_p = np.sqrt(mu / p)
    r = r_mag[..., None] * rhat
    v = (sqrt_mu_p * e * np.sin(nu))[..., None] * rhat + (sqrt_mu_p * (1 + e * np.cos(nu)))[..., None] * that
    return r, v


def cart2kep(r, v, mu):
    """(r, v) -> (a, e, i, raan, argp, nu, u) ; robust for near-circular orbits via u."""
    r = np.atleast_2d(r)
    v = np.atleast_2d(v)
    rn = np.linalg.norm(r, axis=1)
    h = np.cross(r, v)
    hn = np.linalg.norm(h, axis=1)
    i = np.arccos(np.clip(h[:, 2] / hn, -1, 1))
    nvec = np.cross(np.array([0.0, 0.0, 1.0]), h)
    nn = np.linalg.norm(nvec, axis=1)
    raan = np.mod(np.arctan2(nvec[:, 1], nvec[:, 0]), 2 * np.pi)
    evec = (np.cross(v, h) / mu) - r / rn[:, None]
    e = np.linalg.norm(evec, axis=1)
    energy = 0.5 * np.sum(v * v, axis=1) - mu / rn
    a = -mu / (2 * energy)
    # argument of latitude
    u = np.arctan2(np.sum(r * np.cross(h / hn[:, None], nvec / nn[:, None]), axis=1),
                   np.sum(r * nvec / nn[:, None], axis=1))
    argp = np.arctan2(np.sum(evec * np.cross(h / hn[:, None], nvec / nn[:, None]), axis=1),
                      np.sum(evec * nvec / nn[:, None], axis=1))
    nu = np.mod(u - argp, 2 * np.pi)
    return a, e, i, raan, np.mod(argp, 2 * np.pi), nu, np.mod(u, 2 * np.pi)
