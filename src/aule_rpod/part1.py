"""Part I pipeline: docking-window analysis -> CSV products, Gantt, summary."""
from __future__ import annotations

import time
from pathlib import Path

import numpy as np

from . import products as P
from .config import MissionConfig
from .events import intervals as iv
from .relaxation import pass_relaxations
from .scenario import Scenario
from .windows import compute_constraints, docking_windows, pass_diagnostics


def minimum_relaxation_windows(sc, cs, relax, duration_s):
    """Windows obtained with the smallest uniform widening of the illumination band
    that admits at least one window (other constraints unchanged)."""
    finite = [r for r in relax if np.isfinite(r.delta_band_deg)]
    if not finite:
        return None, [], None
    d_star = min(r.delta_band_deg for r in finite)
    d_use = np.ceil(d_star * 100.0) / 100.0              # round up to 0.01 deg
    cs_r = compute_constraints(sc, sun_min=cs.sun_min_deg - d_use, sun_max=cs.sun_max_deg + d_use)
    return d_use, docking_windows(cs_r, duration_s), cs_r


def run(cfg: MissionConfig, outdir: Path, make_figures: bool = True, verbose: bool = True):
    t_wall = time.time()
    outdir = Path(outdir)
    csv_dir, fig_dir = outdir / "csv", outdir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    csv_dir.mkdir(parents=True, exist_ok=True)
    sc = Scenario(cfg)
    D = cfg.constraints.docking_duration_s
    cs = compute_constraints(sc)
    W = docking_windows(cs, D)
    diags = pass_diagnostics(sc, cs, D)
    relax = pass_relaxations(sc, cs, D)
    d_star, W_relaxed, cs_relaxed = minimum_relaxation_windows(sc, cs, relax, D)

    df_pass = P.passes_table(sc, cs, diags, D)
    df_ecl = P.eclipse_table(sc, cs)
    df_band = P.band_table(sc, cs)
    df_cand = P.candidates_table(sc, diags, relax, D)
    df_win = P.windows_table(sc, W, "strict (as specified)")
    df_wrel = P.windows_table(sc, W_relaxed, f"illumination band widened by {d_star:.2f} deg")
    df_ts = P.sun_angle_timeseries(sc, cs, cfg.numerics.timeseries_step_s, full=False)
    df_state = P.sun_angle_timeseries(sc, cs, 60.0, full=True)
    df_ev = P.event_log(sc, cs)

    files = {
        "ground_station_passes": P.write_csv(df_pass, csv_dir / "ground_station_passes.csv"),
        "eclipse_windows": P.write_csv(df_ecl, csv_dir / "eclipse_windows.csv"),
        "sun_angle_band_intervals": P.write_csv(df_band, csv_dir / "sun_angle_band_intervals.csv"),
        "sun_angle_timeseries_10s": P.write_csv(df_ts, csv_dir / "sun_angle_timeseries_10s.csv"),
        "state_timeline_60s": P.write_csv(df_state, csv_dir / "state_timeline_60s.csv"),
        "docking_windows": P.write_csv(df_win, csv_dir / "docking_windows.csv"),
        "docking_window_candidates": P.write_csv(df_cand, csv_dir / "docking_window_candidates.csv"),
        "docking_windows_min_relaxation": P.write_csv(df_wrel, csv_dir / "docking_windows_min_relaxation.csv"),
        "constraint_event_log": P.write_csv(df_ev, csv_dir / "constraint_event_log.csv"),
    }
    if make_figures:
        from .viz import gantt
        relaxed = {f"band widened {d_star:.2f} deg": W_relaxed} if W_relaxed else None
        files["gantt_png"] = gantt.gantt_overview(sc, cs, W, D, fig_dir / "gantt_docking_windows.png",
                                                  relaxed=relaxed, candidates=df_cand)
        files["gantt_html"] = gantt.gantt_interactive(sc, cs, W, D, fig_dir / "gantt_docking_windows_interactive.html",
                                                      relaxed=relaxed, candidates=df_cand)
        files["candidates_png"] = gantt.candidate_panels(sc, cs, df_cand, D, fig_dir / "near_miss_pass_anatomy.png")

    passes_q = {k: int(len(iv.filter_min_duration(v, D))) for k, v in cs.passes.items()}
    first_relaxed = W_relaxed[0] if W_relaxed else None
    summary = {
        "config_fingerprint": cfg.fingerprint(),
        "attitude_file_sha256": sc.profile.meta["sha256"],
        "orbit": {**sc.design.summary(), "propagator": sc.orbit.name,
                  "nodal_period_min": float(sc.orbit.nodal_period_s / 60) if hasattr(sc.orbit, "nodal_period_s") else None},
        "beta_deg_range": [float(x) for x in np.round(np.percentile(sc.beta_deg(np.linspace(sc.t_start, sc.t_end, 2000)), [0, 100]), 3)],
        "passes_total": {k: int(len(v)) for k, v in cs.passes.items()},
        "passes_ge_required_duration": passes_q,
        "eclipses": int(len(df_ecl)),
        "eclipses_cut_by_mission_window": int(df_ecl["truncated_by_mission_window"].sum()),
        "full_eclipse_duration_s_range": [float(df_ecl.loc[~df_ecl.truncated_by_mission_window, "eclipse_duration_s"].min()),
                                          float(df_ecl.loc[~df_ecl.truncated_by_mission_window, "eclipse_duration_s"].max())],
        "full_eclipse_umbra_s_range": [float(df_ecl.loc[~df_ecl.truncated_by_mission_window, "umbra_duration_s"].min()),
                                       float(df_ecl.loc[~df_ecl.truncated_by_mission_window, "umbra_duration_s"].max())],
        "sun_band_fraction_of_time": iv.total(cs.band) / (sc.t_end - sc.t_start),
        "sun_band_intervals_ge_required": int(len(iv.filter_min_duration(cs.band, D))),
        "lit_and_band_intervals_ge_required": int(len(iv.filter_min_duration(iv.intersect(cs.lit, cs.band), D))),
        "strict_docking_windows": len(W),
        "best_candidate": df_cand.iloc[0].to_dict() if len(df_cand) else None,
        "min_uniform_band_widening_deg": d_star,
        "windows_under_min_relaxation": df_wrel.to_dict(orient="records"),
        "first_window_under_min_relaxation_utc": sc.ts.iso(first_relaxed.feasible_start_s)[0] if first_relaxed else None,
        "runtime_s": round(time.time() - t_wall, 2),
        "files": {k: str(v) for k, v in files.items()},
    }
    P.write_json(summary, outdir / "part1_summary.json")
    if verbose:
        print(f"[part1] strict windows: {len(W)} | qualifying passes: {passes_q} | "
              f"min band widening {d_star:.2f} deg -> {len(W_relaxed)} window(s) | {summary['runtime_s']} s")
    return {"scenario": sc, "constraints": cs, "windows": W, "diagnostics": diags, "relaxations": relax,
            "relaxed_windows": W_relaxed, "relaxed_constraints": cs_relaxed, "delta_star": d_star,
            "candidates": df_cand, "summary": summary}
