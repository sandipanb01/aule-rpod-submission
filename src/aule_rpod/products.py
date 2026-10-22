"""Product writers: timestamped CSV files and JSON summaries (Part I).

Every CSV carries ISO-8601 UTC timestamps ('...Z') plus seconds since the
mission-window start, so files can be joined without re-parsing dates.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from .events import intervals as iv


def _iso(sc, t, d=3):
    return sc.ts.iso(t, d)


def _tm(sc, t):
    return np.round(np.asarray(t, dtype=float) - sc.t_start, 3)


def passes_table(sc, cs, diags, duration_s):
    rows = []
    dmap = {(d.station, round(d.aos_s, 3)): d for d in diags}
    pid = 0
    for k, st in enumerate(sc.stations):
        for aos, los in cs.passes[st.name]:
            pid += 1
            el_a, az_a, _ = sc.look_angles(k, np.array([aos]))
            el_l, az_l, _ = sc.look_angles(k, np.array([los]))
            d = dmap.get((st.name, round(aos, 3)))
            rows.append({
                "station": st.name,
                "aos_utc": _iso(sc, aos)[0], "los_utc": _iso(sc, los)[0],
                "tca_utc": _iso(sc, d.tca_s)[0] if d else "",
                "duration_s": round(los - aos, 3), "duration_min": round((los - aos) / 60, 3),
                "max_elevation_deg": round(d.max_elevation_deg, 3) if d else np.nan,
                "aos_azimuth_deg": round(float(az_a[0]), 2), "los_azimuth_deg": round(float(az_l[0]), 2),
                "meets_min_duration": bool(los - aos >= duration_s),
                "truncated_by_mission_window": bool(aos <= sc.t_start + 1e-6 or los >= sc.t_end - 1e-6),
                "aos_t_mission_s": _tm(sc, aos), "los_t_mission_s": _tm(sc, los),
            })
    df = pd.DataFrame(rows).sort_values("aos_t_mission_s").reset_index(drop=True)
    df.insert(0, "pass_id", [f"P{i + 1:03d}" for i in range(len(df))])
    return df


def eclipse_table(sc, cs):
    t0, t1 = cs.span
    ecl = iv.complement(cs.lit, t0, t1)
    rows = []
    for k, (a, b) in enumerate(ecl):
        U = iv.intersect([[a, b]], cs.umbra)
        ua, ub = (U[0] if len(U) else (np.nan, np.nan))
        rows.append({
            "eclipse_id": f"E{k + 1:03d}",
            "penumbra_entry_utc": _iso(sc, a)[0], "umbra_entry_utc": _iso(sc, ua)[0] if np.isfinite(ua) else "",
            "umbra_exit_utc": _iso(sc, ub)[0] if np.isfinite(ub) else "", "penumbra_exit_utc": _iso(sc, b)[0],
            "eclipse_duration_s": round(b - a, 3),
            "umbra_duration_s": round(ub - ua, 3) if np.isfinite(ua) else 0.0,
            "penumbra_in_s": round(ua - a, 3) if np.isfinite(ua) else round(b - a, 3),
            "penumbra_out_s": round(b - ub, 3) if np.isfinite(ub) else 0.0,
            "truncated_by_mission_window": bool(a <= t0 + 1e-6 or b >= t1 - 1e-6),
            "entry_t_mission_s": _tm(sc, a), "exit_t_mission_s": _tm(sc, b),
        })
    return pd.DataFrame(rows)


def band_table(sc, cs):
    rows = []
    for k, (a, b) in enumerate(cs.band):
        lit_overlap = iv.total(iv.intersect([[a, b]], cs.lit))
        rows.append({"band_id": f"B{k + 1:04d}", "start_utc": _iso(sc, a)[0], "end_utc": _iso(sc, b)[0],
                     "duration_s": round(b - a, 3), "sunlit_overlap_s": round(lit_overlap, 3),
                     "start_t_mission_s": _tm(sc, a), "end_t_mission_s": _tm(sc, b)})
    return pd.DataFrame(rows)


def sun_angle_timeseries(sc, cs, step_s, full=False):
    t = np.arange(sc.t_start, sc.t_end + 1e-9, step_s)
    ang = sc.sun_angle_deg(t)
    sL = sc.sun_dir_lvlh(t)
    pL = sc.port_axis_lvlh(t)
    frac = sc.sun_fraction(t)
    lo, hi = cs.sun_min_deg, cs.sun_max_deg
    df = pd.DataFrame({
        "utc": _iso(sc, t, 0), "t_mission_s": _tm(sc, t),
        "sun_angle_deg": np.round(ang, 4),
        "in_band_30_60": ((ang >= lo) & (ang <= hi)).astype(int),
        "sun_visible_fraction": np.round(frac, 5),
        "sunlit": (frac >= 1.0 - 1e-12).astype(int),
    })
    if full:
        df["beta_deg"] = np.round(sc.beta_deg(t), 4)
        for j, ax in enumerate("xyz"):
            df[f"sun_dir_lvlh_{ax}"] = np.round(sL[:, j], 6)
        for j, ax in enumerate("xyz"):
            df[f"port_axis_lvlh_{ax}"] = np.round(pL[:, j], 6)
        vis_any = np.zeros(len(t), bool)
        for k, st in enumerate(sc.stations):
            el = sc.elevation_deg(k, t)
            key = st.name.replace(" ", "")
            df[f"elev_{key}_deg"] = np.round(el, 3)
            vis_any |= el >= cs.min_el_deg
        df["any_station_visible"] = vis_any.astype(int)
        df["all_constraints_instant"] = (vis_any & (df["sunlit"].values == 1) & (df["in_band_30_60"].values == 1)).astype(int)
    return df


def event_log(sc, cs):
    ev = []
    t0, t1 = cs.span
    for name, P in cs.passes.items():
        for a, b in P:
            ev += [(a, "AOS", name, ""), (b, "LOS", name, "")]
    ecl = iv.complement(cs.lit, t0, t1)
    for a, b in ecl:
        ev += [(a, "PENUMBRA_ENTRY", "Earth shadow", ""), (b, "PENUMBRA_EXIT", "Earth shadow", "")]
    for a, b in cs.umbra:
        ev += [(a, "UMBRA_ENTRY", "Earth shadow", ""), (b, "UMBRA_EXIT", "Earth shadow", "")]
    for a, b in cs.sun_ge_min:
        ev += [(a, "SUN_ANGLE_RISES_ABOVE", f"{cs.sun_min_deg:g} deg", ""), (b, "SUN_ANGLE_FALLS_BELOW", f"{cs.sun_min_deg:g} deg", "")]
    for a, b in cs.sun_le_max:
        ev += [(a, "SUN_ANGLE_FALLS_BELOW", f"{cs.sun_max_deg:g} deg", ""), (b, "SUN_ANGLE_RISES_ABOVE", f"{cs.sun_max_deg:g} deg", "")]
    ev = [e for e in ev if t0 + 1e-6 < e[0] < t1 - 1e-6]
    ev.sort(key=lambda e: e[0])
    return pd.DataFrame({"utc": _iso(sc, [e[0] for e in ev]), "t_mission_s": _tm(sc, [e[0] for e in ev]),
                         "event": [e[1] for e in ev], "subject": [e[2] for e in ev]})


def windows_table(sc, windows, scenario_label="strict"):
    cols = ["scenario", "window_id", "station", "pass_aos_utc", "pass_los_utc", "window_open_utc",
            "latest_start_utc", "feasible_end_utc", "feasible_duration_s", "start_slack_s"]
    rows = []
    for k, w in enumerate(windows):
        rows.append({"scenario": scenario_label, "window_id": f"W{k + 1:02d}", "station": w.station,
                     "pass_aos_utc": _iso(sc, w.pass_aos_s)[0], "pass_los_utc": _iso(sc, w.pass_los_s)[0],
                     "window_open_utc": _iso(sc, w.feasible_start_s)[0],
                     "latest_start_utc": _iso(sc, w.latest_start_s)[0],
                     "feasible_end_utc": _iso(sc, w.feasible_end_s)[0],
                     "feasible_duration_s": round(w.feasible_end_s - w.feasible_start_s, 3),
                     "start_slack_s": round(w.slack_s, 3)})
    return pd.DataFrame(rows, columns=cols)


def candidates_table(sc, diags, relax, duration_s):
    rmap = {(r.station, round(r.aos_s, 3)): r for r in relax}
    rows = []
    for d in diags:
        if not d.qualifies_duration:
            continue
        r = rmap.get((d.station, round(d.aos_s, 3)))
        status = "VALID" if d.is_window else "NOT VALID"
        rows.append({
            "station": d.station, "pass_aos_utc": _iso(sc, d.aos_s)[0], "pass_los_utc": _iso(sc, d.los_s)[0],
            "pass_duration_s": round(d.duration_s, 3), "max_elevation_deg": round(d.max_elevation_deg, 2),
            "status": status,
            "longest_feasible_start_utc": _iso(sc, d.best_start_s)[0] if np.isfinite(d.best_start_s) else "",
            "longest_feasible_end_utc": _iso(sc, d.best_end_s)[0] if np.isfinite(d.best_end_s) else "",
            "longest_feasible_s": round(d.best_feasible_s, 3),
            "deficit_vs_required_s": round(d.deficit_s, 3),
            "limited_at_start_by": d.limited_start_by, "limited_at_end_by": d.limited_end_by,
            "sun_angle_min_in_pass_deg": round(d.sun_angle_min_deg, 2),
            "sun_angle_max_in_pass_deg": round(d.sun_angle_max_deg, 2),
            "sunlit_fraction_of_pass": round(d.lit_fraction, 4),
            "min_band_widening_for_window_deg": round(r.delta_band_deg, 3) if r and np.isfinite(r.delta_band_deg) else np.nan,
            "min_lower_bound_relaxation_deg": round(r.delta_lo_deg, 3) if r and np.isfinite(r.delta_lo_deg) else np.nan,
            "relaxed_window_start_utc": _iso(sc, r.best_start_band_s)[0] if r and np.isfinite(r.delta_band_deg) else "",
        })
    df = pd.DataFrame(rows)
    if len(df):
        df = df.sort_values(["deficit_vs_required_s", "min_band_widening_for_window_deg"]).reset_index(drop=True)
        df.insert(0, "rank", np.arange(1, len(df) + 1))
    return df


def write_csv(df: pd.DataFrame, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    return path


def write_json(obj, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)

    def default(o):
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, (np.bool_,)):
            return bool(o)
        return str(o)

    path.write_text(json.dumps(obj, indent=2, default=default))
    return path
