"""Reference frames.

GCRS  -- geocentric celestial reference system (ICRS-aligned axes); the
         inertial frame in which the Sun and both spacecraft are expressed.
TEOE  -- "true equator & CIO of orbit epoch" (CIRS frozen at the orbit epoch);
         the quasi-inertial frame in which the target orbit is propagated, so
         that the J2 symmetry axis is the true (CIP) pole. Equivalent to STK's
         'TrueOfOrbitEpoch' propagation frame. TEOE->GCRS is a constant matrix.
ITRS  -- terrestrial frame of the ground stations.
LVLH  -- STK definition (as used by the attitude file):
         X = r/|r| (radial, zenith), Z = (r x v)/|r x v| (orbit normal),
         Y = Z x X (along-track, toward velocity).

GCRS->ITRS follows the IAU 2006/2000A CIO-based chain implemented by ERFA:
    [ITRS] = W(xp, yp, s') . R3(ERA) . C2I(X, Y, s) . [GCRS]
C2I varies by < 1 mas/min, so it is evaluated exactly on 10-min nodes and
interpolated element-wise (error < 1e-11 rad, verified in tests); ERA is
exact at every time.
"""
from __future__ import annotations

import erfa
import numpy as np

from .timesys import TimeSystem


def _rz_frame(angle):
    """Frame rotation about z by ``angle`` (ERFA convention, vectorised)."""
    c, s = np.cos(angle), np.sin(angle)
    R = np.zeros(np.shape(angle) + (3, 3))
    R[..., 0, 0] = c
    R[..., 0, 1] = s
    R[..., 1, 0] = -s
    R[..., 1, 1] = c
    R[..., 2, 2] = 1.0
    return R


class EarthOrientation:
    def __init__(self, ts: TimeSystem, t_start: float, t_end: float, node_step_s: float = 600.0):
        self.ts = ts
        self._step = node_step_s
        self._nodes = np.arange(t_start - 2 * node_step_s, t_end + 3 * node_step_s, node_step_s)
        tt1, tt2 = ts.jd_tt(self._nodes)
        self._c2i_nodes = erfa.c2i06a(tt1, tt2)

    def c2i(self, t):
        t = np.asarray(t, dtype=float)
        if np.any(t < self._nodes[0]) or np.any(t > self._nodes[-1]):
            tt1, tt2 = self.ts.jd_tt(t)
            return erfa.c2i06a(tt1, tt2)
        x = (t - self._nodes[0]) / self._step
        i = np.clip(np.floor(x).astype(int), 0, len(self._nodes) - 2)
        f = (x - i)[..., None, None]
        return (1.0 - f) * self._c2i_nodes[i] + f * self._c2i_nodes[i + 1]

    def era(self, t):
        u1, u2 = self.ts.jd_ut1(t)
        return erfa.era00(u1, u2)

    def pom(self, t):
        xp, yp = self.ts.polar_motion_rad(t)
        tt1, tt2 = self.ts.jd_tt(t)
        return erfa.pom00(xp, yp, erfa.sp00(tt1, tt2))

    def gcrs_to_itrs(self, t):
        """Rotation matrices M (..., 3, 3) with r_itrs = M @ r_gcrs."""
        return self.pom(t) @ _rz_frame(self.era(t)) @ self.c2i(t)

    def teoe_to_gcrs(self, t_epoch: float):
        """Constant matrix mapping TEOE (CIRS frozen at t_epoch) to GCRS."""
        return np.asarray(self.c2i(t_epoch)).T


# --------------------------------------------------------------------- WGS-84
def geodetic_to_itrs(lat_deg, lon_deg, alt_km, re_km, f):
    """Geodetic (WGS-84) -> ITRS position [km] and local geodetic 'up' unit vector."""
    lat, lon = np.radians(lat_deg), np.radians(lon_deg)
    e2 = f * (2.0 - f)
    N = re_km / np.sqrt(1.0 - e2 * np.sin(lat) ** 2)
    r = np.array([(N + alt_km) * np.cos(lat) * np.cos(lon),
                  (N + alt_km) * np.cos(lat) * np.sin(lon),
                  (N * (1.0 - e2) + alt_km) * np.sin(lat)])
    up = np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])
    east = np.array([-np.sin(lon), np.cos(lon), 0.0])
    north = np.cross(up, east)
    return r, up, east, north


def itrs_to_geodetic(r, re_km, f):
    """ITRS [km] -> geodetic latitude/longitude [deg] and altitude [km] (Bowring, 3 its)."""
    r = np.atleast_2d(r)
    x, y, z = r[:, 0], r[:, 1], r[:, 2]
    e2 = f * (2.0 - f)
    p = np.hypot(x, y)
    lon = np.arctan2(y, x)
    lat = np.arctan2(z, p * (1.0 - e2))
    for _ in range(5):
        N = re_km / np.sqrt(1.0 - e2 * np.sin(lat) ** 2)
        alt = p / np.cos(lat) - N
        lat = np.arctan2(z, p * (1.0 - e2 * N / (N + alt)))
    N = re_km / np.sqrt(1.0 - e2 * np.sin(lat) ** 2)
    alt = p / np.cos(lat) - N
    return np.degrees(lat), np.degrees(lon), alt


# ----------------------------------------------------------------------- LVLH
def lvlh_matrix(r, v):
    """STK LVLH: rows are (X, Y, Z) unit vectors expressed in the inertial frame,
    so that  v_lvlh = M @ v_inertial  and  v_inertial = M.T @ v_lvlh."""
    r = np.atleast_2d(r)
    v = np.atleast_2d(v)
    x = r / np.linalg.norm(r, axis=-1, keepdims=True)
    h = np.cross(r, v)
    z = h / np.linalg.norm(h, axis=-1, keepdims=True)
    y = np.cross(z, x)
    return np.stack([x, y, z], axis=-2)
