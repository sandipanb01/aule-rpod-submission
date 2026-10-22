"""Scenario: wires configuration -> time, frames, orbit, attitude, Sun, stations,
and exposes every constraint as a continuous, vectorised event function g(t)
(g >= 0  <=>  constraint satisfied) for the event finder.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .attitude import make_attitude_provider, read_stk_attitude
from .attitude.quaternion import angle_between_deg
from .config import MissionConfig
from .ephemeris import SunERFA
from .frames import EarthOrientation, lvlh_matrix
from .geometry.eclipse import g_full_sun, g_umbra, sun_fraction
from .geometry.stations import Station
from .orbit.cowell import CowellJ2Orbit, fit_osculating_to_mean
from .orbit.elements import kep2cart
from .orbit.secular_j2 import SecularJ2Orbit
from .orbit.sso import raan_from_ltan_teoe, sso_inclination_rad


@dataclass
class OrbitDesign:
    a_km: float
    e: float
    i_rad: float
    raan_teoe_rad: float
    u0_rad: float
    era_epoch_rad: float
    ut1_hours_epoch: float

    def summary(self):
        return {"a_km": self.a_km, "altitude_km_equatorial": None, "e": self.e,
                "i_deg": float(np.degrees(self.i_rad)), "raan_teoe_deg": float(np.degrees(self.raan_teoe_rad)),
                "arg_latitude_epoch_deg": float(np.degrees(self.u0_rad)),
                "era_epoch_deg": float(np.degrees(self.era_epoch_rad)), "ut1_hours_epoch": self.ut1_hours_epoch}


class Scenario:
    def __init__(self, cfg: MissionConfig, attitude_kind: str | None = None,
                 propagator: str | None = None, eop: bool = True):
        from .timesys import TimeSystem
        self.cfg = cfg
        sc, ea = cfg.scenario, cfg.earth
        self.ts = TimeSystem(sc.orbit_epoch_utc, sc.mission_end_utc, eop=eop)
        self.t_start = self.ts.t_from_iso(sc.mission_start_utc)
        self.t_end = self.ts.t_from_iso(sc.mission_end_utc)
        self.eo = EarthOrientation(self.ts, -3600.0, self.t_end + 3600.0)
        self.R_teoe_gcrs = self.eo.teoe_to_gcrs(0.0)
        self.re, self.mu, self.j2, self.f = ea.equatorial_radius_km, ea.mu_km3_s2, ea.j2, ea.flattening
        # ---------------------------------------------------------------- orbit
        to = cfg.target_orbit
        a = self.re + to.altitude_km
        inc = sso_inclination_rad(a, to.eccentricity, self.mu, self.re, self.j2, ea.tropical_year_days)
        era0 = float(self.eo.era(0.0))
        ut1h = self.ts.ut1_hours_of_day(0.0)
        raan = raan_from_ltan_teoe(to.ltan_hours, ut1h, era0)
        u0 = np.radians(to.arg_latitude_at_epoch_deg)
        self.design = OrbitDesign(a, to.eccentricity, inc, raan, u0, era0, ut1h)
        kind = propagator or to.propagator
        self.propagator_kind = kind
        if kind == "secular_j2":
            self.orbit = SecularJ2Orbit(a, to.eccentricity, inc, raan, 0.0, u0, self.mu, self.re, self.j2)
        elif kind == "cowell_osculating":
            r0, v0 = kep2cart(a, to.eccentricity, inc, raan, 0.0, u0, self.mu)
            self.orbit = CowellJ2Orbit(r0, v0, self.mu, self.re, self.j2, self.t_end + 7200.0,
                                       name="Cowell J2, elements as osculating")
        elif kind == "cowell_mean":
            r0, v0 = fit_osculating_to_mean(a, inc, raan, u0, self.mu, self.re, self.j2)
            self.orbit = CowellJ2Orbit(r0, v0, self.mu, self.re, self.j2, self.t_end + 7200.0,
                                       name="Cowell J2, elements as mean (numerically fitted)")
        else:
            raise ValueError(kind)
        # ------------------------------------------------------------- attitude
        self.profile = read_stk_attitude(cfg.attitude_path)
        if self.profile.scenario_epoch[:19] != sc.orbit_epoch_utc[:19]:
            raise ValueError("Attitude ScenarioEpoch differs from the orbit epoch: "
                             f"{self.profile.scenario_epoch} vs {sc.orbit_epoch_utc}")
        if not self.profile.coordinate_axes.upper().startswith("LVLH"):
            raise ValueError(f"Expected LVLH reference axes, got {self.profile.coordinate_axes}")
        self.attitude_kind = attitude_kind or cfg.attitude.interpolation
        self.att = make_attitude_provider(self.attitude_kind, self.profile)
        p = np.asarray(cfg.attitude.docking_port_body, dtype=float)
        self.port_body_raw_norm = float(np.linalg.norm(p))
        self.port_body = p / self.port_body_raw_norm
        # ---------------------------------------------------------- Sun/stations
        self.sun = SunERFA(self.ts, -3600.0, self.t_end + 3600.0)
        self.stations = [Station(s.name, s.lat_deg, s.lon_deg, s.alt_m, self.re, self.f)
                         for s in cfg.ground_stations]
        self.shadow_f = self.f            # WGS-84 ellipsoid for Earth shadow (0 -> sphere)
        self._cache_key = None
        self._cache = None

    # ======================================================== state evaluation
    def state_gcrs(self, t):
        r, v = self.orbit.state_teoe(t)
        R = self.R_teoe_gcrs
        return np.atleast_2d(r) @ R.T, np.atleast_2d(v) @ R.T

    def _states(self, t):
        """Memoised bundle of everything the constraint functions need at times t."""
        t = np.atleast_1d(np.asarray(t, dtype=float))
        key = (t.shape, float(t[0]), float(t[-1]), float(t.sum()))
        if self._cache_key == key:
            return self._cache
        r, v = self.state_gcrs(t)
        s = self.sun.position_gcrs_km(t)
        M = self.eo.gcrs_to_itrs(t)
        r_itrs = np.einsum("nij,nj->ni", M, r)
        L = lvlh_matrix(r, v)
        out = {"t": t, "r": r, "v": v, "s": s, "r_itrs": r_itrs, "L": L}
        self._cache_key, self._cache = key, out
        return out

    def sun_dir_lvlh(self, t):
        st = self._states(t)
        d = st["s"] - st["r"]
        d = d / np.linalg.norm(d, axis=1, keepdims=True)
        return np.einsum("nij,nj->ni", st["L"], d)

    def port_axis_lvlh(self, t):
        return self.att.rotation(np.atleast_1d(t)).apply(self.port_body)

    def sun_angle_deg(self, t):
        return angle_between_deg(self.port_axis_lvlh(t), self.sun_dir_lvlh(t))

    def beta_deg(self, t):
        st = self._states(t)
        h = np.cross(st["r"], st["v"])
        h /= np.linalg.norm(h, axis=1, keepdims=True)
        sh = st["s"] / np.linalg.norm(st["s"], axis=1, keepdims=True)
        return np.degrees(np.arcsin(np.sum(h * sh, 1)))

    def elevation_deg(self, k, t):
        return self.stations[k].elevation_deg(self._states(t)["r_itrs"])

    def look_angles(self, k, t):
        return self.stations[k].look_angles(self._states(t)["r_itrs"])

    def _shadow_inputs(self, t):
        """Spacecraft and Sun in TEOE (z = true pole) for the ellipsoidal shadow model."""
        st = self._states(t)
        R = self.R_teoe_gcrs
        return st["r"] @ R, st["s"] @ R

    def g_full_sun(self, t):
        r, s = self._shadow_inputs(t)
        return g_full_sun(r, s, self.re, self.shadow_f)

    def g_umbra(self, t):
        r, s = self._shadow_inputs(t)
        return g_umbra(r, s, self.re, self.shadow_f)

    def sun_fraction(self, t):
        r, s = self._shadow_inputs(t)
        return sun_fraction(r, s, self.re, self.shadow_f)

    # ================================================== constraint functions
    def constraint_functions(self, sun_min=None, sun_max=None, min_el=None):
        c = self.cfg.constraints
        lo = c.sun_angle_min_deg if sun_min is None else sun_min
        hi = c.sun_angle_max_deg if sun_max is None else sun_max
        el0 = c.min_elevation_deg if min_el is None else min_el
        g = {
            "sun_angle_ge_min": lambda t: self.sun_angle_deg(t) - lo,
            "sun_angle_le_max": lambda t: hi - self.sun_angle_deg(t),
            "full_sun": self.g_full_sun,
            "not_umbra": lambda t: -self.g_umbra(t),
        }
        for k, s in enumerate(self.stations):
            g[f"vis::{s.name}"] = (lambda kk: (lambda t: self.elevation_deg(kk, t) - el0))(k)
        return g

    def lighting_function(self):
        return self.g_full_sun if self.cfg.constraints.eclipse_excludes_penumbra else (lambda t: -self.g_umbra(t))
