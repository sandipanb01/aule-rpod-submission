"""Runs the full V&V campaign and writes validation_report.{json,md}."""
from __future__ import annotations

import time
from pathlib import Path

import pandas as pd

from .. import products as P
from ..config import REPO_ROOT
from ..scenario import Scenario
from . import crosscheck as X


def run(cfg, outdir: Path, de421_dir=None, orekit_data=None, candidates_csv=None, verbose=True):
    t0 = time.time()
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    sc = Scenario(cfg)
    rep = {}
    rep["V1_attitude"] = X.check_attitude(sc)
    rep["V2_sun_ephemeris"] = X.check_sun(sc, de421_dir)
    rep["V3_frames"] = X.check_frames(sc)
    rep["V4_orbit"] = X.check_orbit(cfg)
    rep["V5_event_grid"] = X.check_event_grid(sc)
    if candidates_csv and Path(candidates_csv).exists():
        rep["V6_end_to_end_independent"] = X.check_end_to_end_astropy(sc, pd.read_csv(candidates_csv).fillna(""))
    od = Path(orekit_data or REPO_ROOT / "data" / "orekit-data")
    try:
        from . import orekit_xval
        rep["V7_orekit"] = orekit_xval.run(sc, str(od))
    except Exception as exc:  # optional dependency
        rep["V7_orekit"] = {"verdict": "SKIPPED", "reason": repr(exc)}
    rep["runtime_s"] = round(time.time() - t0, 1)
    P.write_json(rep, outdir / "validation_report.json")
    lines = ["| ID | Check | Key result | Verdict |", "|---|---|---|---|"]
    a = rep["V1_attitude"]
    lines.append(f"| V1 | Attitude law & quaternion convention | fit residual {a['law_fit_max_residual_rad']:.1e} rad; "
                 f"SLERP vs law {a['slerp_vs_law_port_max_deg']:.3f} deg; NLERP vs law {a['nlerp_vs_law_port_max_deg']:.3f} deg | {a['verdict']} |")
    s = rep["V2_sun_ephemeris"]
    lines.append(f"| V2 | Sun: ERFA vs JPL DE421 / analytic | DE421 {s.get('de421_vs_erfa_arcsec', float('nan')):.2f}\" ; "
                 f"almanac {s['almanac_vs_erfa_arcsec']:.1f}\" ; interpolation {s['interp_vs_direct_arcsec']:.1e}\" | {s['verdict']} |")
    f = rep["V3_frames"]
    lines.append(f"| V3 | GCRS->ITRS vs ERFA direct & astropy | {f['c2t_interp_vs_erfa_direct_rad']:.1e} rad; "
                 f"ITRS pos {f.get('itrs_position_vs_astropy_m', float('nan')):.3f} m; elev {f.get('elevation_vs_astropy_deg', float('nan')):.1e} deg | {f['verdict']} |")
    o = rep["V4_orbit"]
    lines.append(f"| V4 | Secular J2 vs Cowell J2; LTAN flown | mean reading: {o['cowell_mean']['phase_drift_deg_per_day']:+.4f} deg/day phase drift; "
                 f"osculating reading: {o['cowell_osculating']['phase_drift_deg_per_day']:+.2f} deg/day; LTAN flown "
                 f"{o['ltan_flown_h_range'][0]:.5f}-{o['ltan_flown_h_range'][1]:.5f} h over {o['node_crossings_checked']} nodes | {o['verdict']} |")
    e = rep["V5_event_grid"]
    lines.append(f"| V5 | Event-finder grid independence (1/2/5 s) | max diff {max(v for k, v in e.items() if k.startswith('max')):.1e} s | {e['verdict']} |")
    if "V6_end_to_end_independent" in rep:
        r6 = rep["V6_end_to_end_independent"]
        w = max(r["difference_s"] for r in r6["rows"])
        lines.append(f"| V6 | Independent astropy chain, top candidates | longest-span difference <= {w:.2f} s | {r6['verdict']} |")
    k = rep["V7_orekit"]
    if k.get("verdict") != "SKIPPED":
        lines.append(f"| V7 | Orekit 13: propagation, elevation, AOS/LOS, eclipse | {k['cowell_vs_orekit_max_position_diff_m']:.3f} m; "
                     f"{k['elevation_vs_orekit_max_deg']:.1e} deg; AOS/LOS {k['aos_los_vs_orekit_max_s'] * 1e3:.2f} ms; "
                     f"penumbra {k.get('penumbra_events_vs_orekit_max_s', float('nan')) * 1e3:.1f} ms | {k['verdict']} |")
    else:
        lines.append(f"| V7 | Orekit | {k.get('reason', '')[:60]} | SKIPPED |")
    (outdir / "validation_report.md").write_text("\n".join(lines) + "\n")
    if verbose:
        print("\n".join(lines))
        print(f"[validation] {rep['runtime_s']} s")
    return rep
