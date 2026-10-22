"""Pure-Python cross-checks against independent implementations.

Each check returns a dict with the measured discrepancy, the tolerance it is
judged against, and PASS/FAIL. Tolerances are set by what the analysis needs
(event times to ~1 s; angles to ~0.01 deg), not by what makes a check pass.
"""
from __future__ import annotations

import warnings
from pathlib import Path

import erfa
import numpy as np

from ..attitude.identification import identify_precession_law
from ..attitude.providers import LawAttitude, NlerpAttitude, SlerpAttitude
from ..ephemeris import SunAlmanac, SunDE421
from ..events import intervals as iv
from ..scenario import Scenario
from ..windows import compute_constraints

ARCSEC = np.degrees(1) * 3600


def _verdict(value, tol):
    return "PASS" if value <= tol else "FAIL"


def check_sun(sc, de421_dir: str | Path | None):
    t = np.arange(sc.t_start, sc.t_end, 3600.0)
    p_int = sc.sun.position_gcrs_km(t)
    p_dir = sc.sun._direct(t)[0]

    def ang(a, b):
        c = np.sum(a * b, 1) / np.linalg.norm(a, axis=1) / np.linalg.norm(b, axis=1)
        return np.arccos(np.clip(c, -1, 1)) * ARCSEC

    out = {"interp_vs_direct_arcsec": float(ang(p_int, p_dir).max())}
    out["almanac_vs_erfa_arcsec"] = float(ang(SunAlmanac(sc.ts).position_gcrs_km(t), p_dir).max())
    if de421_dir and Path(de421_dir).exists():
        de = SunDE421(sc.ts, de421_dir).position_gcrs_km(t)
        out["de421_vs_erfa_arcsec"] = float(ang(de, p_dir).max())
        out["de421_vs_erfa_range_km"] = float(np.abs(np.linalg.norm(de, axis=1) - np.linalg.norm(p_dir, axis=1)).max())
    # light-time/aberration effect we neglect: annual aberration bound
    out["neglected_aberration_arcsec"] = 20.5
    worst = max(out.get("de421_vs_erfa_arcsec", 0), out["interp_vs_direct_arcsec"])
    out["verdict"] = _verdict(worst, 36.0)     # 0.01 deg
    return out


def check_frames(sc, n=400):
    rng = np.random.default_rng(1)
    t = np.sort(rng.uniform(sc.t_start, sc.t_end, n))
    M_fast = sc.eo.gcrs_to_itrs(t)
    tt1, tt2 = sc.ts.jd_tt(t)
    u1, u2 = sc.ts.jd_ut1(t)
    xp, yp = sc.ts.polar_motion_rad(t)
    M_ref = erfa.c2t06a(tt1, tt2, u1, u2, xp, yp)
    dR = np.einsum("nij,nkj->nik", M_fast, M_ref)
    from scipy.spatial.transform import Rotation
    ang = Rotation.from_matrix(dR).magnitude()
    out = {"c2t_interp_vs_erfa_direct_rad": float(ang.max())}
    # astropy: independent coordinate-transformation graph and IERS handling
    try:
        from astropy import units as u
        from astropy.coordinates import GCRS, ITRS, CartesianRepresentation
        from astropy.time import Time
        from astropy.utils import iers
        iers.conf.auto_download = False
        r, _ = sc.state_gcrs(t)
        tim = Time(sc.ts.iso(t, 6).astype(str).tolist(), scale="utc")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            g = GCRS(CartesianRepresentation(r[:, 0] * u.km, r[:, 1] * u.km, r[:, 2] * u.km), obstime=tim)
            ri = g.transform_to(ITRS(obstime=tim)).cartesian.xyz.to_value(u.km).T
        mine = np.einsum("nij,nj->ni", M_fast, r)
        out["itrs_position_vs_astropy_m"] = float(np.linalg.norm(ri - mine, axis=1).max() * 1000)
        dels = []
        for k, st in enumerate(sc.stations):
            dels.append(np.abs(st.elevation_deg(ri) - st.elevation_deg(mine)).max())
        out["elevation_vs_astropy_deg"] = float(max(dels))
    except Exception as exc:  # pragma: no cover
        out["astropy_error"] = repr(exc)
    out["verdict"] = _verdict(out.get("elevation_vs_astropy_deg", 0.0), 1e-3)
    return out


def check_orbit(cfg):
    """Secular-J2 baseline vs Cowell J2 under both readings of '550 km circular'
    (elements as MEAN -> fitted osculating state; elements as OSCULATING), plus
    the local mean solar time actually flown at every ascending node."""
    from ..orbit.elements import cart2kep
    base = Scenario(cfg)
    t = np.arange(0.0, base.t_end + 1, 30.0)
    rb, vb = base.orbit.state_teoe(t)
    ub = np.unwrap(cart2kep(rb, vb, base.mu)[6])
    out = {}
    for kind in ("cowell_mean", "cowell_osculating"):
        sc = Scenario(cfg, propagator=kind)
        r, v = sc.orbit.state_teoe(t)
        kk = cart2kep(r, v, sc.mu)
        du = np.degrees(np.unwrap(kk[6]) - ub)
        du -= 360.0 * np.round(du[0] / 360.0)
        rate = np.polyfit(t / 86400.0, du, 1)[0]
        out[kind] = {"phase_diff_deg_at_mission_start": float(np.interp(base.t_start, t, du)),
                     "phase_diff_deg_at_mission_end": float(du[-1]),
                     "phase_drift_deg_per_day": float(rate),
                     "radius_diff_km_range": [float((np.linalg.norm(r, axis=1) - np.linalg.norm(rb, axis=1)).min()),
                                              float((np.linalg.norm(r, axis=1) - np.linalg.norm(rb, axis=1)).max())],
                     "mean_sma_km": float(np.mean(kk[0]))}
    # exact ascending-node crossings of the baseline orbit -> local mean solar time flown
    m = np.floor(ub / (2 * np.pi))
    idx = np.where(np.diff(m) > 0)[0]
    lmts = []
    for i in idx:
        target = 2 * np.pi * m[i + 1]
        tn = t[i] + (target - ub[i]) / (ub[i + 1] - ub[i]) * (t[i + 1] - t[i])
        rn, _ = base.state_gcrs(np.array([tn]))
        ri = base.eo.gcrs_to_itrs(np.array([tn]))[0] @ rn[0]
        lon = np.degrees(np.arctan2(ri[1], ri[0]))
        lmts.append((base.ts.ut1_hours_of_day(tn) + lon / 15.0) % 24.0)
    lmts = np.array(lmts)
    out["ltan_flown_h_range"] = [float(lmts.min()), float(lmts.max())]
    out["ltan_drift_s_over_19_days"] = float((lmts.max() - lmts.min()) * 3600.0)
    out["node_crossings_checked"] = int(len(lmts))
    out["verdict"] = _verdict(max(abs(lmts - cfg.target_orbit.ltan_hours)) * 3600, 5.0)
    return out


def check_attitude(sc):
    law = identify_precession_law(sc.profile)
    t = np.arange(sc.t_start, sc.t_start + 86400, 1.0)
    p = sc.port_body
    pl = LawAttitude(sc.profile, law).rotation(t).apply(p)

    def dev(prov):
        q = prov.rotation(t).apply(p)
        return float(np.degrees(np.arccos(np.clip(np.sum(q * pl, 1), -1, 1))).max())

    out = {"law_fit_max_residual_rad": law.max_residual_rad, "law_fit_rms_residual_rad": law.rms_residual_rad,
           "slerp_vs_law_port_max_deg": dev(SlerpAttitude(sc.profile)),
           "nlerp_vs_law_port_max_deg": dev(NlerpAttitude(sc.profile)),
           "law": law.summary()}
    # convention unit check: q = 90 deg about Z => body X along +Y_ref
    from ..attitude.quaternion import attitude_matrix, body_to_reference
    q = np.array([0, 0, np.sin(np.pi / 4), np.cos(np.pi / 4)])
    bx = body_to_reference(q, [1, 0, 0])
    out["convention_test_bodyX_in_ref"] = np.round(bx, 12).tolist()
    out["convention_test_matches_shuster_A_transpose"] = bool(np.allclose(attitude_matrix(q).T @ [1, 0, 0], bx))
    out["verdict"] = "PASS" if (law.max_residual_rad < 1e-8 and np.allclose(bx, [0, 1, 0])
                                and out["convention_test_matches_shuster_A_transpose"]) else "FAIL"
    return out


def check_event_grid(sc):
    res = {}
    sets = {h: compute_constraints(sc, step=h) for h in (1.0, 2.0, 5.0)}

    def diff(a, b):
        a, b = np.asarray(a), np.asarray(b)
        if a.shape != b.shape:
            return np.inf
        return float(np.abs(a - b).max()) if a.size else 0.0

    worst = 0.0
    for h in (1.0, 5.0):
        d = [diff(sets[h].band, sets[2.0].band), diff(sets[h].lit, sets[2.0].lit), diff(sets[h].umbra, sets[2.0].umbra)]
        d += [diff(sets[h].passes[k], sets[2.0].passes[k]) for k in sets[2.0].passes]
        res[f"max_event_time_diff_{h:g}s_vs_2s_grid_s"] = max(d)
        worst = max(worst, max(d))
    res["n_band_intervals"] = {str(h): int(len(s.band)) for h, s in sets.items()}
    res["verdict"] = _verdict(worst, 0.01)
    return res


def check_end_to_end_astropy(sc, candidates, n=3, dt=0.5):
    """Recompute the longest feasible span of the top-n candidate passes with an
    independent chain: astropy GCRS->ITRS for elevation, astropy (ERFA/VSOP)
    Sun, closed-form attitude law, own LVLH construction; 0.5-s sampling."""
    from astropy import units as u
    from astropy.coordinates import GCRS, ITRS, CartesianRepresentation, get_body_barycentric
    from astropy.time import Time
    from astropy.utils import iers
    iers.conf.auto_download = False
    law = LawAttitude(sc.profile)
    rows = []
    for _, c in candidates.head(n).iterrows():
        a = sc.ts.t_from_iso(c["pass_aos_utc"].rstrip("Z")) - 30
        b = sc.ts.t_from_iso(c["pass_los_utc"].rstrip("Z")) + 30
        t = np.arange(a, b, dt)
        tim = Time(sc.ts.iso(t, 6).astype(str).tolist(), scale="utc")
        r, v = sc.state_gcrs(t)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sun = (get_body_barycentric("sun", tim) - get_body_barycentric("earth", tim)).xyz.to_value(u.km).T
            ri = GCRS(CartesianRepresentation(r[:, 0] * u.km, r[:, 1] * u.km, r[:, 2] * u.km),
                      obstime=tim).transform_to(ITRS(obstime=tim)).cartesian.xyz.to_value(u.km).T
        st = [s for s in sc.stations if s.name == c["station"]][0]
        el = st.elevation_deg(ri)
        x = r / np.linalg.norm(r, axis=1, keepdims=True)
        z = np.cross(r, v)
        z /= np.linalg.norm(z, axis=1, keepdims=True)
        y = np.cross(z, x)
        s = sun - r
        s /= np.linalg.norm(s, axis=1, keepdims=True)
        sL = np.stack([np.sum(s * x, 1), np.sum(s * y, 1), np.sum(s * z, 1)], 1)
        ang = np.degrees(np.arccos(np.clip(np.sum(law.rotation(t).apply(sc.port_body) * sL, 1), -1, 1)))
        from ..geometry.eclipse import g_full_sun
        Rt = sc.R_teoe_gcrs
        ok = (el >= 5) & (ang >= 30) & (ang <= 60) & (g_full_sun(r @ Rt, sun @ Rt, sc.re, sc.shadow_f) >= 0)
        m = np.r_[False, ok, False].astype(int)
        dd = np.diff(m)
        spans = (np.where(dd == -1)[0] - 1 - np.where(dd == 1)[0]) * dt
        best = float(spans.max()) if len(spans) else 0.0
        rows.append({"station": c["station"], "pass_aos_utc": c["pass_aos_utc"],
                     "pipeline_longest_feasible_s": float(c["longest_feasible_s"]),
                     "independent_longest_feasible_s": best,
                     "difference_s": abs(best - float(c["longest_feasible_s"]))})
    worst = max(r["difference_s"] for r in rows)
    return {"rows": rows, "verdict": _verdict(worst, 2.0)}
