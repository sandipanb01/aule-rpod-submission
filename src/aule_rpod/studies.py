"""Orchestrates Part I sensitivity/robustness studies and writes their products."""
from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pandas as pd

from . import products as P
from .attitude.identification import identify_precession_law
from .scenario import Scenario
from .sensitivity import feasibility_map, relaxation_curves, run_alternatives
from .viz import figures as F
from .windows import compute_constraints


def run(cfg, outdir: Path, quick: bool = False, verbose: bool = True, reuse_maps: bool = False):
    t_wall = time.time()
    outdir = Path(outdir)
    (outdir / "csv").mkdir(parents=True, exist_ok=True)
    (outdir / "figures").mkdir(parents=True, exist_ok=True)
    sc = Scenario(cfg)
    D = cfg.constraints.docking_duration_s
    res = {}
    # ---- attitude characterisation & geometry explainer
    law = identify_precession_law(sc.profile)
    P.write_json(law.summary(), outdir / "attitude_law.json")
    F.attitude_characterisation(sc, law, outdir / "figures" / "attitude_characterisation.png")
    cs = compute_constraints(sc)
    F.geometry_explainer(sc, cs, outdir / "figures" / "geometry_why_only_wellington.png")
    # ---- feasibility maps
    if quick:
        ltans, u0s = np.arange(6.0, 18.01, 1.0), np.arange(0.0, 360.0, 30.0)
        ltz, u0z = np.arange(9.5, 11.51, 0.25), np.arange(0.0, 360.0, 30.0)
    else:
        ltans, u0s = np.arange(0.0, 24.0, 0.25), np.arange(0.0, 360.0, 10.0)
        ltz, u0z = np.arange(9.0, 12.001, 0.05), np.arange(0.0, 360.0, 5.0)
    cache = outdir / "feasibility_maps.npz"
    if reuse_maps and cache.exists():
        z = np.load(cache)
        ltans, u0s, n_full, b_full = z["ltans"], z["u0s"], z["n_full"], z["b_full"]
        ltz, u0z, n_zoom, b_zoom = z["ltz"], z["u0z"], z["n_zoom"], z["b_zoom"]
    else:
        n_full, b_full = feasibility_map(sc, ltans, u0s, D)
        n_zoom, b_zoom = feasibility_map(sc, ltz, u0z, D)
        np.savez(cache, ltans=ltans, u0s=u0s, n_full=n_full, b_full=b_full,
                 ltz=ltz, u0z=u0z, n_zoom=n_zoom, b_zoom=b_zoom)
    F.feasibility_map(ltans, u0s, n_full, outdir / "figures" / "feasibility_map_ltan_full_day.png",
                      title="Strict-window count over all LTANs (10-day mission, given attitude profile)")
    F.feasibility_map(ltz, u0z, n_zoom, outdir / "figures" / "feasibility_map_ltan_zoom.png", best=None,
                      title="Zoom around the specified LTAN 10:30: windows need LTAN <= ~10:10 (u0-dependent)")
    k0 = int(np.argmin(np.abs(ltz - 10.5)))
    zoom_row = {"ltan_10_30_fraction_of_phases_with_window": float(np.mean(n_zoom[:, k0] > 0)),
                "ltan_10_30_best_span_s_over_phases": [float(b_zoom[:, k0].min()), float(b_zoom[:, k0].max())]}
    df_map = pd.DataFrame([{"ltan_h": L, "u0_deg": u, "strict_windows": int(n_zoom[i, j]),
                            "best_span_s": round(float(b_zoom[i, j]), 1)}
                           for j, L in enumerate(ltz) for i, u in enumerate(u0z)])
    P.write_csv(df_map, outdir / "csv" / "feasibility_map_ltan_phase.csv")
    # ---- relaxation curves
    df_rel = relaxation_curves(sc, D)
    P.write_csv(df_rel, outdir / "csv" / "relaxation_curves.csv")
    F.relaxation_curves(df_rel, outdir / "figures" / "relaxation_curves.png")
    # ---- alternative interpretations
    df_alt, df_altw = run_alternatives(cfg, D)
    P.write_csv(df_alt, outdir / "csv" / "alternative_interpretations_summary.csv")
    P.write_csv(df_altw, outdir / "csv" / "alternative_interpretations_windows.csv")
    res.update({"zoom": zoom_row, "relaxation": df_rel, "alternatives": df_alt, "alt_windows": df_altw,
                "runtime_s": round(time.time() - t_wall, 1)})
    P.write_json({"zoom": zoom_row, "alternatives": df_alt.to_dict(orient="records"),
                  "relaxation": df_rel.to_dict(orient="records"), "runtime_s": res["runtime_s"]},
                 outdir / "sensitivity_summary.json")
    if verbose:
        print(df_alt[["scenario", "strict_windows", "windows_by_station", "first_window_open_utc",
                      "longest_feasible_span_s"]].to_string(index=False))
        print(f"[studies] done in {res['runtime_s']} s")
    return res
