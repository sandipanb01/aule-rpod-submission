"""Sun ephemerides (geocentric, GCRS axes, km).

Primary  : ERFA ``epv00`` (IAU SOFA; VSOP2000-fitted series). Earth-Sun vector
           accurate to a few km (~1e-8 rad) for 1900-2100.
V&V only : JPL DE421 Chebyshev ephemeris (Folkner et al. 2009) read from the
           ``de421`` package arrays, and the Astronomical Almanac low-precision
           formula (~0.01 deg) as an analytic sanity check.

Geometric (not light-time/aberration corrected) positions are used; the
apparent direction differs by <= 20.5 arcsec (0.006 deg), far below every
tolerance in this analysis (quantified in the validation report).
"""
from __future__ import annotations

from pathlib import Path

import erfa
import numpy as np

from .constants import AU_KM
from .timesys import TimeSystem


class SunERFA:
    """ERFA epv00, evaluated exactly on 10-min nodes and interpolated with cubic
    Hermite splines using epv00's own velocities (interpolation error < 1 m,
    i.e. < 1e-11 rad; verified against direct evaluation in the tests)."""
    name = "ERFA epv00"

    def __init__(self, ts: TimeSystem, t_start: float | None = None, t_end: float | None = None,
                 node_step_s: float = 600.0):
        self.ts = ts
        self._nodes = None
        if t_start is not None and t_end is not None:
            self._h = node_step_s
            self._nodes = np.arange(t_start - 2 * node_step_s, t_end + 3 * node_step_s, node_step_s)
            self._p, self._v = self._direct(self._nodes)

    def _direct(self, t):
        tt1, tt2 = self.ts.jd_tt(t)
        pvh, _ = erfa.epv00(tt1, tt2)                       # heliocentric Earth [au, au/d]
        p = -np.asarray(pvh["p"]) * AU_KM
        v = -np.asarray(pvh["v"]) * AU_KM / 86400.0
        return p, v

    def position_gcrs_km(self, t):
        t = np.asarray(t, dtype=float)
        if self._nodes is None or np.any(t < self._nodes[0]) or np.any(t > self._nodes[-1]):
            return self._direct(t)[0]
        x = (t - self._nodes[0]) / self._h
        i = np.clip(np.floor(x).astype(int), 0, len(self._nodes) - 2)
        s = (x - i)[..., None]
        h = self._h
        p0, p1, v0, v1 = self._p[i], self._p[i + 1], self._v[i] * h, self._v[i + 1] * h
        h00 = 2 * s ** 3 - 3 * s ** 2 + 1
        h10 = s ** 3 - 2 * s ** 2 + s
        h01 = -2 * s ** 3 + 3 * s ** 2
        h11 = s ** 3 - s ** 2
        return h00 * p0 + h10 * v0 + h01 * p1 + h11 * v1


class SunAlmanac:
    """Astronomical Almanac low-precision Sun (ecliptic of date -> mean equator).
    Accuracy ~0.01 deg 1950-2050; used only as an independent cross-check."""
    name = "Astronomical Almanac (analytic)"

    def __init__(self, ts: TimeSystem):
        self.ts = ts

    def position_gcrs_km(self, t):
        tt1, tt2 = self.ts.jd_tt(t)
        n = (tt1 - 2451545.0) + tt2
        L = np.radians((280.460 + 0.9856474 * n) % 360.0)
        g = np.radians((357.528 + 0.9856003 * n) % 360.0)
        lam = L + np.radians(1.915) * np.sin(g) + np.radians(0.020) * np.sin(2 * g)
        eps = np.radians(23.439 - 0.0000004 * n)
        R = (1.00014 - 0.01671 * np.cos(g) - 0.00014 * np.cos(2 * g)) * AU_KM
        mod = np.stack([R * np.cos(lam), R * np.cos(eps) * np.sin(lam), R * np.sin(eps) * np.sin(lam)], -1)
        # rotate mean-of-date equator/equinox -> GCRS (precession, IAU 2006)
        P = erfa.pmat06(tt1, tt2)                      # GCRS -> mean of date
        return np.einsum("...ji,...j->...i", P, mod)


class SunDE421:
    """JPL DE421 from the PyPI 'de421' package (Chebyshev arrays, ICRF, km, TDB)."""
    name = "JPL DE421"

    def __init__(self, ts: TimeSystem, data_dir: str | Path):
        self.ts = ts
        d = Path(data_dir)
        const = {k.decode(): v for k, v in np.load(d / "constants.npy", allow_pickle=True)}
        self.jalpha, self.jomega = float(const["jalpha"]), float(const["jomega"])
        self.earth_share = 1.0 / (1.0 + float(const["EMRAT"]))
        self._sets = {k: np.load(d / f"jpl-{k}.npy") for k in ("sun", "earthmoon", "moon")}

    def _cheb(self, name, jd):
        c = self._sets[name]
        nsets, _, ncoef = c.shape
        days = (self.jomega - self.jalpha) / nsets
        idx, off = np.divmod(jd - self.jalpha, days)
        idx = idx.astype(int)
        T = np.empty((ncoef,) + jd.shape)
        T[0] = 1.0
        T[1] = x = 2.0 * off / days - 1.0
        for k in range(2, ncoef):
            T[k] = 2.0 * x * T[k - 1] - T[k - 2]
        coeff = c[idx]                                  # (N, 3, ncoef)
        return np.einsum("nak,kn->na", coeff, T)

    def position_gcrs_km(self, t):
        tt1, tt2 = self.ts.jd_tt(t)
        jd = np.atleast_1d(tt1 + tt2)
        earth = self._cheb("earthmoon", jd) - self._cheb("moon", jd) * self.earth_share
        return self._cheb("sun", jd) - earth
