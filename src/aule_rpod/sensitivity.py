"""Sensitivity, robustness and alternative-interpretation studies (Part I).

1. Feasibility map  - number of strict windows vs (LTAN, initial argument of
   latitude). Uses a fast vectorised counter: everything that does not depend
   on the orbit (Sun, attitude in LVLH, Earth rotation) is precomputed on a
   10-s grid; per case only the closed-form secular-J2 orbit is re-evaluated.
   All constraints are folded into one continuous margin per station,
       m_s(t) = min(el_s - el_min, theta - lo, hi - theta, g_full_sun)   [deg],
   whose zero crossings are linearly interpolated (window lengths to ~0.1 s).
2. Relaxation curves - exact (event-based) window counts vs the illumination
   lower bound, a uniform band widening, the elevation mask and the required
   duration.
3. Alternative interpretations - the full exact pipeline re-run under each
   plausible alternative reading of the specification.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .events import intervals as iv
from .frames import lvlh_matrix
from .geometry.eclipse import g_full_sun
from .orbit.secular_j2 import SecularJ2Orbit
from .windows import compute_constraints, docking_windows


class FastWindowCounter:
    def __init__(self, sc, step: float = 10.0, port_axis_lvlh=None):
        self.sc = sc
        self.step = step
        self.t = np.arange(sc.t_start, sc.t_end + 1e-9, step)
        self.sun = sc.sun.position_gcrs_km(self.t)
        self.pL = sc.port_axis_lvlh(self.t) if port_axis_lvlh is None else port_axis_lvlh
        self.M = sc.eo.gcrs_to_itrs(self.t)
        self.R = sc.R_teoe_gcrs
        d = sc.design
        self.base = d
        self.ltan0 = sc.cfg.target_orbit.ltan_hours

    def margins(self, ltan_h, u0_deg, lo=30.0, hi=60.0, min_el=5.0):
        d, sc = self.base, self.sc
        raan = d.raan_teoe_rad + np.radians(15.0 * (ltan_h - self.ltan0))
        orb = SecularJ2Orbit(d.a_km, d.e, d.i_rad, raan, 0.0, np.radians(u0_deg), sc.mu, sc.re, sc.j2)
        r, v = orb.state_teoe(self.t)
        r, v = r @ self.R.T, v @ self.R.T
        L = lvlh_matrix(r, v)
        s = self.sun - r
        s /= np.linalg.norm(s, axis=1, keepdims=True)
        sL = np.einsum("nij,nj->ni", L, s)
        ang = np.degrees(np.arccos(np.clip(np.sum(sL * self.pL, 1), -1, 1)))
        common = np.minimum.reduce([ang - lo, hi - ang, g_full_sun(r @ self.R, self.sun @ self.R, sc.re, sc.shadow_f)])
        r_itrs = np.einsum("nij,nj->ni", self.M, r)
        out = {}
        for st in sc.stations:
            dvec = r_itrs - st.r_itrs
            el = np.degrees(np.arcsin(np.clip(dvec @ st.up / np.linalg.norm(dvec, axis=1), -1, 1)))
            out[st.name] = np.minimum(el - min_el, common)
        return out

    def _runs(self, m):
        pos = m >= 0
        x = np.r_[False, pos, False].astype(np.int8)
        dd = np.diff(x)
        starts, ends = np.where(dd == 1)[0], np.where(dd == -1)[0] - 1
        t, h = self.t, self.step
        res = []
        for i, j in zip(starts, ends):
            a = t[i] if i == 0 else t[i - 1] + h * (0 - m[i - 1]) / (m[i] - m[i - 1])
            b = t[j] if j == len(t) - 1 else t[j] + h * (0 - m[j]) / (m[j + 1] - m[j])
            res.append((a, b))
        return res

    def windows(self, ltan_h, u0_deg, D=540.0, **kw):
        res = []
        for name, m in self.margins(ltan_h, u0_deg, **kw).items():
            for a, b in self._runs(m):
                if b - a >= D:
                    res.append((name, a, b))
        return sorted(res, key=lambda x: x[1])

    def best_span(self, ltan_h, u0_deg, **kw):
        best = 0.0
        for name, m in self.margins(ltan_h, u0_deg, **kw).items():
            for a, b in self._runs(m):
                best = max(best, b - a)
        return best


def feasibility_map(sc, ltans, u0s, D=540.0, step=10.0):
    fc = FastWindowCounter(sc, step)
    n = np.zeros((len(u0s), len(ltans)), int)
    best = np.zeros_like(n, dtype=float)
    for j, L in enumerate(ltans):
        for i, u in enumerate(u0s):
            m = fc.margins(L, u)
            spans = [b - a for mm in m.values() for a, b in fc._runs(mm)]
            spans = np.array(spans) if spans else np.zeros(1)
            n[i, j] = int(np.sum(spans >= D))
            best[i, j] = spans.max()
    return n, best


def relaxation_curves(sc, D=540.0):
    rows = []
    base = compute_constraints(sc)
    for lo in np.arange(20.0, 30.01, 1.0):
        cs = compute_constraints(sc, sun_min=lo)
        W = docking_windows(cs, D)
        rows.append({"study": "sun_angle_lower_bound_deg", "value": lo, "windows": len(W),
                     "first_window_utc": sc.ts.iso(W[0].feasible_start_s)[0] if W else ""})
    for dlt in np.arange(0.0, 8.01, 1.0):
        cs = compute_constraints(sc, sun_min=30 - dlt, sun_max=60 + dlt)
        W = docking_windows(cs, D)
        rows.append({"study": "uniform_band_widening_deg", "value": dlt, "windows": len(W),
                     "first_window_utc": sc.ts.iso(W[0].feasible_start_s)[0] if W else ""})
    # the exact minimum widening (Part I logic), so the curve shows where the first window appears
    from .part1 import minimum_relaxation_windows
    from .relaxation import pass_relaxations
    d_star, W_star, _ = minimum_relaxation_windows(sc, base, pass_relaxations(sc, base, D), D)
    if W_star:
        rows.append({"study": "uniform_band_widening_deg", "value": float(d_star), "windows": len(W_star),
                     "first_window_utc": sc.ts.iso(W_star[0].feasible_start_s)[0], "kind": "minimum"})
    for el in np.arange(0.0, 5.01, 1.0):
        cs = compute_constraints(sc, min_el=el)
        W = docking_windows(cs, D)
        rows.append({"study": "min_elevation_deg", "value": el, "windows": len(W),
                     "first_window_utc": sc.ts.iso(W[0].feasible_start_s)[0] if W else ""})
    for dmin in np.arange(5.0, 9.01, 0.5):
        W = docking_windows(base, dmin * 60.0)
        rows.append({"study": "required_duration_min", "value": dmin, "windows": len(W),
                     "first_window_utc": sc.ts.iso(W[0].feasible_start_s)[0] if W else ""})
    df = pd.DataFrame(rows)
    df["kind"] = df["kind"].fillna("grid") if "kind" in df else "grid"
    return df


@dataclass
class AltScenario:
    key: str
    description: str
    cfg_updates: dict = field(default_factory=dict)
    attitude_kind: str | None = None
    propagator: str | None = None
    transpose_quaternion: bool = False
    penumbra_allowed: bool = False
    true_sun_ltan: bool = False


ALTERNATIVES = [
    AltScenario("baseline", "As specified: LTAN 10:30 (mean Sun), STK quaternion convention, SLERP, secular J2"),
    AltScenario("ltdn_1030", "LTDN 10:30 instead of LTAN 10:30 (i.e. LTAN 22:30)", {"target_orbit": {"ltan_hours": 22.5}}),
    AltScenario("quat_transposed", "Quaternion read as body->reference (convention mistake test)", transpose_quaternion=True),
    AltScenario("att_nlerp", "Attitude interpolated as STK 'Lagrange order 1' (NLERP)", attitude_kind="nlerp"),
    AltScenario("att_law", "Attitude from the identified closed-form law", attitude_kind="law"),
    AltScenario("cowell_osc", "Numerical J2, entered elements taken as osculating", propagator="cowell_osculating"),
    AltScenario("cowell_mean", "Numerical J2, entered elements taken as mean", propagator="cowell_mean"),
    AltScenario("penumbra_ok", "Eclipse = umbra only (penumbra allowed)", {"constraints": {"eclipse_excludes_penumbra": False}}),
    AltScenario("true_sun_ltan", "LTAN referenced to the true (apparent) Sun", true_sun_ltan=True),
    AltScenario("u0_90", "Initial argument of latitude 90 deg (phase sensitivity)", {"target_orbit": {"arg_latitude_at_epoch_deg": 90.0}}),
    AltScenario("u0_180", "Initial argument of latitude 180 deg", {"target_orbit": {"arg_latitude_at_epoch_deg": 180.0}}),
    AltScenario("u0_270", "Initial argument of latitude 270 deg", {"target_orbit": {"arg_latitude_at_epoch_deg": 270.0}}),
]


class _Transposed:
    name = "transposed quaternion"

    def __init__(self, inner):
        self.inner = inner

    def rotation(self, t):
        return self.inner.rotation(t).inv()


def build_alt_scenario(cfg, alt: AltScenario):
    from .scenario import Scenario
    c = cfg.with_updates(**alt.cfg_updates) if alt.cfg_updates else cfg
    sc = Scenario(c, attitude_kind=alt.attitude_kind, propagator=alt.propagator)
    if alt.transpose_quaternion:
        sc.att = _Transposed(sc.att)
    if alt.true_sun_ltan:
        # RAAN measured from the apparent Sun instead of the mean Sun (differs by the equation of time)
        s = sc.sun.position_gcrs_km(np.array([0.0]))[0]
        s_teoe = sc.R_teoe_gcrs.T @ s
        ra_sun = np.arctan2(s_teoe[1], s_teoe[0])
        raan_true = np.mod(ra_sun + np.radians(15.0 * (c.target_orbit.ltan_hours - 12.0)), 2 * np.pi)
        d = sc.design
        sc.eot_deg = float(np.degrees((raan_true - d.raan_teoe_rad + np.pi) % (2 * np.pi) - np.pi))
        sc.orbit = SecularJ2Orbit(d.a_km, d.e, d.i_rad, raan_true, 0.0, d.u0_rad, sc.mu, sc.re, sc.j2)
        d.raan_teoe_rad = raan_true
        sc._cache_key = None
    return sc


def run_alternatives(cfg, D=540.0):
    rows, wins = [], []
    for alt in ALTERNATIVES:
        sc = build_alt_scenario(cfg, alt)
        cs = compute_constraints(sc)
        W = docking_windows(cs, D)
        best = 0.0
        for P in cs.passes.values():
            for a, b in iv.filter_min_duration(P, D):
                F = iv.intersect_all([[a, b]], cs.lit, cs.band)
                if len(F):
                    best = max(best, iv.durations(F).max())
        by_station = {}
        for w in W:
            by_station[w.station] = by_station.get(w.station, 0) + 1
        rows.append({"scenario": alt.key, "description": alt.description, "strict_windows": len(W),
                     "windows_by_station": ";".join(f"{k}:{v}" for k, v in by_station.items()),
                     "first_window_open_utc": sc.ts.iso(W[0].feasible_start_s)[0] if W else "",
                     "longest_feasible_span_s": round(best, 2),
                     "extra": f"EoT offset {getattr(sc, 'eot_deg', 0.0):+.3f} deg" if alt.true_sun_ltan else ""})
        for k, w in enumerate(W):
            wins.append({"scenario": alt.key, "window_id": f"W{k + 1:02d}", "station": w.station,
                         "window_open_utc": sc.ts.iso(w.feasible_start_s)[0],
                         "latest_start_utc": sc.ts.iso(w.latest_start_s)[0],
                         "feasible_end_utc": sc.ts.iso(w.feasible_end_s)[0],
                         "feasible_duration_s": round(w.feasible_end_s - w.feasible_start_s, 3)})
    return pd.DataFrame(rows), pd.DataFrame(wins)
