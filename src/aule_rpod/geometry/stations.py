"""Ground-station topocentric geometry (WGS-84 geodetic up, no refraction)."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..frames import geodetic_to_itrs


@dataclass
class Station:
    name: str
    lat_deg: float
    lon_deg: float
    alt_m: float
    re_km: float
    f: float

    def __post_init__(self):
        self.r_itrs, self.up, self.east, self.north = geodetic_to_itrs(
            self.lat_deg, self.lon_deg, self.alt_m / 1000.0, self.re_km, self.f)

    def look_angles(self, r_sat_itrs):
        """Elevation [deg], azimuth [deg, from N through E], range [km]."""
        d = np.atleast_2d(r_sat_itrs) - self.r_itrs
        rng = np.linalg.norm(d, axis=-1)
        el = np.degrees(np.arcsin(np.clip(d @ self.up / rng, -1, 1)))
        az = np.degrees(np.mod(np.arctan2(d @ self.east, d @ self.north), 2 * np.pi))
        return el, az, rng

    def elevation_deg(self, r_sat_itrs):
        d = np.atleast_2d(r_sat_itrs) - self.r_itrs
        return np.degrees(np.arcsin(np.clip((d @ self.up) / np.linalg.norm(d, axis=-1), -1, 1)))
