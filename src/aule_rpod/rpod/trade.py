"""Approach-strategy trade at the selected docking window.

Every candidate is flown through the same timeline (same leg times, hold and
final speed) so that the comparison isolates the approach axis. Metrics:
  * delta-v (impulsive + continuous) and time to contact vs window length,
  * illumination of the port during the critical phase (hold + final) and Sun
    angle off the chaser camera boresight; camera background (Earth vs space),
  * far-field passive safety (miss distance if B1/B2 are missed),
  * whether the transfer has to pass close to the target (fly-around needed),
  * station-keeping cost of the hold point (V-bar points are equilibria),
  * plume exposure of the braking impulses.
A 'sun-optimal' free axis is found by maximising the minimum margin to the
30-60 deg band over the critical phase (with an Earth-free camera background).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .design import Strategy, build_plan, illumination_metrics, plume_metrics, safety_metrics

AXES = [("vbar_plus", "+V-bar (from ahead)", [0, 1, 0]),
        ("vbar_minus", "-V-bar (from behind)", [0, -1, 0]),
        ("rbar_below", "+R-bar (from below, Earth side)", [-1, 0, 0]),
        ("rbar_above", "-R-bar (from above)", [1, 0, 0]),
        ("hbar_minus", "H-bar (from -orbit-normal side)", [0, 0, -1]),
        ("hbar_plus", "H-bar (from +orbit-normal side)", [0, 0, 1])]


def sun_optimal_axis(sun_fn, t_abs, band=(30.0, 60.0), earth_limit_deg=82.0, step_deg=2.0):
    s = sun_fn(t_abs)
    best, best_axis = -np.inf, None
    for el in np.arange(-90.0, 90.01, step_deg):
        for az in np.arange(0.0, 360.0, step_deg):
            a = np.array([np.cos(np.radians(el)) * np.cos(np.radians(az)),
                          np.cos(np.radians(el)) * np.sin(np.radians(az)), np.sin(np.radians(el))])
            nadir_to_bore = np.degrees(np.arccos(np.clip(-a @ np.array([-1.0, 0, 0]), -1, 1)))
            if nadir_to_bore < earth_limit_deg:
                continue
            th = np.degrees(np.arccos(np.clip(s @ a, -1, 1)))
            m = np.min(np.minimum(th - band[0], band[1] - th))
            if m > best:
                best, best_axis = m, a
    return best_axis, best


def run_trade(n, X0, sun_fn, t0_abs, window_s, timeline, r_ko=5.0, band=(30.0, 60.0), cam_excl=30.0,
              earth_ang_deg=67.0):
    strategies = [Strategy(k, nm, np.asarray(ax, float)) for k, nm, ax in AXES]
    T_crit = np.linspace(timeline["T1"] + timeline["T2"],
                         timeline["T1"] + timeline["T2"] + timeline["T_hold"] + timeline["rho_hold"] / timeline["v_final"], 60)
    ax_opt, margin = sun_optimal_axis(sun_fn, t0_abs + T_crit, band)
    if ax_opt is not None:
        strategies.append(Strategy("sun_optimal", "Sun-optimal inclined axis", ax_opt))
    rows, plans = [], {}
    for st in strategies:
        p = build_plan(n, X0, st, timeline["T1"], timeline["T2"], timeline["T_hold"], timeline["rho_entry"],
                       timeline["rho_hold"], timeline["v_final"], timeline["v_contact"], timeline["d_brake"])
        plans[st.key] = p
        il = illumination_metrics(p, sun_fn, t0_abs, band, cam_excl, earth_ang_deg)
        sf = safety_metrics(n, p, r_ko=r_ko)
        pl = plume_metrics(p)
        braking_in_plume = [r["burn"] for r in pl if r["target_in_plume_cone"] and r["burn"] in ("B2", "B3", "B5")]
        hp_acc = np.linalg.norm(np.array([3 * n * n * p.params["HP"][0], 0.0, -n * n * p.params["HP"][2]]))
        rows.append({
            "strategy": st.key, "name": st.name, "axis_lvlh": np.round(st.axis, 4).tolist(),
            "dv_total_mps": p.dv_total, "dv_impulsive_mps": p.dv_impulsive,
            "dv_continuous_mps": float(sum(p.dv_continuous.values())),
            "t_contact_s": p.t_contact, "window_margin_s": window_s - p.t_contact,
            "sun_on_port_range_deg": [round(il["sun_on_port_min_deg"], 1), round(il["sun_on_port_max_deg"], 1)],
            "sun_on_port_in_band_frac": il["sun_on_port_in_band_fraction"],
            "port_face_sunlit": il["port_face_sunlit"],
            "sun_off_camera_min_deg": il["sun_off_camera_boresight_min_deg"],
            "camera_background": il["camera_background"],
            "miss_if_B1_missed_m": sf["miss_if_B1_fails_m"], "miss_if_B2_missed_m": sf["miss_if_B2_fails_m"],
            "leg1_min_range_m": sf["leg1_min_range_m"],
            "transfer_passes_inside_KOS": bool(sf["leg1_min_range_m"] < timeline["rho_entry"] - 1.0),
            "hold_point_equilibrium": bool(hp_acc < 1e-9),
            "hold_station_keeping_mm_s": 1000 * p.dv_continuous.get("hold station-keeping", 0.0),
            "braking_impulses_in_plume_cone": ",".join(braking_in_plume) or "none",
        })
    df = pd.DataFrame(rows)
    return df, plans, {"sun_optimal_axis": None if ax_opt is None else ax_opt.tolist(), "sun_optimal_margin_deg": margin}


def score(df):
    """Transparent pass/fail screening (no hidden weights): a strategy is
    eligible only if it is safe in the far field, keeps the camera clear of the
    Sun, lights the port inside the band for the whole critical phase, fits in
    the window and needs no fly-around; eligible strategies are ranked by dV."""
    ok = ((df.miss_if_B2_missed_m > 5.0) & (df.sun_off_camera_min_deg >= 30.0) & (df.port_face_sunlit)
          & (df.sun_on_port_in_band_frac >= 0.999) & (df.window_margin_s > 0) & (~df.transfer_passes_inside_KOS)
          & (df.camera_background == "deep space"))
    df = df.copy()
    df["eligible"] = ok
    df["rank"] = df.dv_total_mps.where(ok).rank(method="min")
    return df
