"""Part II pipeline: approach-strategy selection and maneuver design at the first
docking window identified in Part I, with analytic CW delta-v, nonlinear J2
validation, safety analysis and Monte Carlo dispersions."""
from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pandas as pd

from . import products as P
from .events import intervals as iv
from .part1 import minimum_relaxation_windows
from .relaxation import pass_relaxations
from .rpod import cw, montecarlo, nonlinear, safety
from .rpod.design import Strategy, build_plan, plume_metrics, safety_metrics
from .rpod.trade import run_trade, score
from .scenario import Scenario
from .windows import compute_constraints, docking_windows


def select_window(sc, cfg):
    """First strict window if any; otherwise the first window under the minimum
    uniform widening of the illumination band (all other constraints kept)."""
    D = cfg.constraints.docking_duration_s
    cs = compute_constraints(sc)
    W = docking_windows(cs, D)
    if W:
        w, label, d = W[0], "strict", 0.0
    else:
        relax = pass_relaxations(sc, cs, D)
        d, Wr, _ = minimum_relaxation_windows(sc, cs, relax, D)
        w, label = Wr[0], f"minimum relaxation (illumination band widened by {d:.2f} deg)"
    return {"t0": w.feasible_start_s, "t_end": w.feasible_end_s, "station": w.station, "label": label,
            "band_widening_deg": d, "pass_aos": w.pass_aos_s, "pass_los": w.pass_los_s, "lit": cs.lit}


def analytic_worksheet(n, X0, plan):
    """Every intermediate quantity of the analytic CW solution (for the report)."""
    p = plan.params
    out = {"n_rad_s": n, "orbital_period_s": 2 * np.pi / n}
    r0, v0 = np.asarray(X0[:3], float), np.asarray(X0[3:], float)
    for leg, (ra, va, rb, T) in {"leg1": (r0, v0, p["WP1"], p["T1"]),
                                 "leg2": (p["WP1"], plan.burns[0].v_after * 0, p["HP"], p["T2"])}.items():
        rr, rv, vr, vv = cw.stm_blocks(n, np.array(T))
        if leg == "leg2":
            va = plan.burns[1].v_before
        v_dep = np.linalg.solve(rv, rb - rr @ ra)
        v_arr = vr @ ra + vv @ v_dep
        out[leg] = {"T_s": T, "nT_rad": n * T, "cos_nT": float(np.cos(n * T)), "sin_nT": float(np.sin(n * T)),
                    "Phi_rr": rr.round(9).tolist(), "Phi_rv": rv.round(6).tolist(),
                    "Phi_vr": vr.round(12).tolist(), "Phi_vv": vv.round(9).tolist(),
                    "r_start_m": ra.tolist(), "v_before_mps": np.asarray(va).round(6).tolist(),
                    "r_target_m": np.asarray(rb).tolist(), "v_departure_mps": v_dep.round(6).tolist(),
                    "v_arrival_mps": v_arr.round(6).tolist()}
    out["burns"] = [{"id": b.id, "dv_mps": b.dv.round(6).tolist(), "dv_mag_mps": round(b.mag, 6), "purpose": b.purpose}
                    for b in plan.burns]
    L = p["rho_hold"]
    out["final_translation"] = {"length_m": L, "speed_mps": p["v_final"], "duration_s": L / p["v_final"],
                                "radial_accel_required_mps2": 2 * n * p["v_final"],
                                "continuous_dv_2nL_mps": 2 * n * L,
                                "numerical_integral_mps": plan.dv_continuous["final translation compensation"]}
    out["dv_impulsive_total_mps"] = plan.dv_impulsive
    out["dv_continuous_total_mps"] = float(sum(plan.dv_continuous.values()))
    out["dv_total_mps"] = plan.dv_total
    return out


def run(cfg, outdir: Path, make_figures=True, verbose=True, mc_samples=None):
    t_wall = time.time()
    outdir = Path(outdir)
    (outdir / "csv").mkdir(parents=True, exist_ok=True)
    (outdir / "figures").mkdir(parents=True, exist_ok=True)
    p2 = cfg.part2
    sc = Scenario(cfg)
    win = select_window(sc, cfg)
    t0 = win["t0"]
    window_s = win["t_end"] - t0
    r_t, v_t = sc.state_gcrs(np.array([t0]))
    r_t, v_t = r_t[0], v_t[0]
    n = float(np.linalg.norm(np.cross(r_t, v_t)) / np.dot(r_t, r_t))
    X0 = np.r_[p2.initial_relative_position_lvlh_m, p2.initial_relative_velocity_lvlh_mps].astype(float)
    sun_fn = lambda t: sc.sun_dir_lvlh(np.atleast_1d(t))          # noqa: E731
    utc = lambda t: sc.ts.iso(t)[0]                                # noqa: E731
    vbar = Strategy("vbar_plus", "+V-bar (from ahead)", np.array([0.0, 1.0, 0.0]))
    # --------------------------------------------------- timeline budget & trade
    T_final = p2.hold_point_distance_m / p2.final_approach_speed_mps
    budget = window_s - p2.window_margin_s - p2.hold_duration_s - T_final
    rows = []
    for T1 in np.arange(0.45 * budget, 0.86 * budget, 2.0):
        T2 = budget - T1
        pl = build_plan(n, X0, vbar, T1, T2, p2.hold_duration_s, p2.keep_out_sphere_radius_m, p2.hold_point_distance_m,
                        p2.final_approach_speed_mps, p2.final_approach_speed_mps, 0.0)
        sf = safety_metrics(n, pl, r_ko=p2.collision_radius_m)
        rows.append({"T1": float(T1), "T2": float(T2), "dv_total": pl.dv_total, "B1": pl.burns[0].mag,
                     "B2": pl.burns[1].mag, "B3": pl.burns[2].mag, "miss_B2": sf["miss_if_B2_fails_m"],
                     "t_contact": pl.t_contact})
    df_tl = pd.DataFrame(rows)
    feas = df_tl[(df_tl.B3 <= 0.15) & (df_tl.miss_B2 >= 2 * p2.collision_radius_m)]
    chosen = (feas if len(feas) else df_tl).sort_values("dv_total").iloc[0].to_dict()
    timeline = dict(T1=round(chosen["T1"]), T2=round(budget - round(chosen["T1"])), T_hold=p2.hold_duration_s,
                    rho_entry=p2.keep_out_sphere_radius_m, rho_hold=p2.hold_point_distance_m, v_final=p2.final_approach_speed_mps,
                    v_contact=p2.final_approach_speed_mps, d_brake=0.0)
    # --------------------------------------------------- strategy trade
    df_tr, plans, extra = run_trade(n, X0, sun_fn, t0, window_s, timeline, r_ko=p2.collision_radius_m,
                                    cam_excl=p2.camera_sun_exclusion_deg)
    df_tr = score(df_tr)
    best = df_tr[df_tr.eligible].sort_values("rank").iloc[0]["strategy"]
    plan = plans[best]
    # --------------------------------------------------- safety, validation, MC
    drift = safety.missed_burn_drift(n, plan)
    cam = safety.cam_scan(n, plan, dv_cam=0.15, r_ko=p2.collision_radius_m)
    R = sc.R_teoe_gcrs
    nl_log, nl = nonlinear.simulate(plan, r_t @ R, v_t @ R, sc.mu, sc.re, sc.j2, X0)
    mc = p2.monte_carlo
    mc_sum, mc_raw = montecarlo.run(n, X0, plan, mc_samples or mc.samples, mc.seed, mc.nav_position_sigma_frac_of_range,
                                    mc.nav_position_sigma_floor_m, mc.nav_velocity_sigma_mps,
                                    mc.burn_magnitude_sigma_frac, mc.burn_pointing_sigma_deg,
                                    r_collision=p2.collision_radius_m, r_kos=p2.keep_out_sphere_radius_m,
                                    corridor_deg=p2.corridor_half_angle_deg, dt_ctrl=mc.control_period_s,
                                    tau_lat=mc.lateral_time_constant_s, v_lat_max=mc.lateral_speed_limit_mps,
                                    deadband=mc.thruster_deadband_mps, tau_nav=mc.nav_correlation_time_s)
    ws = analytic_worksheet(n, X0, plan)
    # robustness: same design logic at the other candidate windows
    rob = [design_at(sc, cfg, t0, win["t_end"], X0, f"selected ({win['label']})")]
    for lab, a_iso, b_iso in (("first window at a round 4 deg widening (Jul 16)", "2026-07-16T22:56:29.857", "2026-07-16T23:06:16.458"),
                              ("longest strictly-feasible span, 7.7 min (Jul 12)", "2026-07-12T23:12:07.051", "2026-07-12T23:19:49.186")):
        rob.append(design_at(sc, cfg, sc.ts.t_from_iso(a_iso), sc.ts.t_from_iso(b_iso), X0, lab))
    # initial-velocity sensitivity: chaser on a natural circular orbit 80 m below (drifting ahead)
    yd = -1.5 * n * X0[0]
    alt = build_plan(n, np.r_[X0[:3], 0.0, yd, 0.0], plan.strategy, timeline["T1"], timeline["T2"],
                     timeline["T_hold"], p2.keep_out_sphere_radius_m, timeline["rho_hold"], timeline["v_final"], timeline["v_contact"], 0.0)
    # --------------------------------------------------- products
    burns_rows = []
    for b in plan.burns:
        burns_rows.append({"burn_id": b.id, "epoch_utc": utc(t0 + b.t_s), "t_from_window_start_s": round(b.t_s, 3),
                           "dv_x_radial_mps": round(b.dv[0], 6), "dv_y_alongtrack_mps": round(b.dv[1], 6),
                           "dv_z_crosstrack_mps": round(b.dv[2], 6), "dv_magnitude_mps": round(b.mag, 6),
                           "pos_x_m": round(b.r[0], 3), "pos_y_m": round(b.r[1], 3), "pos_z_m": round(b.r[2], 3),
                           "range_m": round(float(np.linalg.norm(b.r)), 3),
                           "v_before_mps": np.round(b.v_before, 5).tolist(), "v_after_mps": np.round(b.v_after, 5).tolist(),
                           "purpose": b.purpose})
    for tag, t_mcc, leg in (("MCC1", 0.5 * plan.params["T1"], "leg1"),
                            ("MCC2", plan.params["T1"] + 0.5 * plan.params["T2"], "leg2")):
        q99 = mc_sum[f"{leg}_mcc_dv_mm_s_mean_p99"][1]
        burns_rows.append({"burn_id": tag, "epoch_utc": utc(t0 + t_mcc), "t_from_window_start_s": round(t_mcc, 3),
                           "dv_x_radial_mps": 0.0, "dv_y_alongtrack_mps": 0.0, "dv_z_crosstrack_mps": 0.0,
                           "dv_magnitude_mps": 0.0,
                           "purpose": f"mid-course correction slot: nominal 0, re-targets {'WP1' if leg == 'leg1' else 'HP'} "
                                      f"from the navigated state (Monte Carlo p99 {q99:.1f} mm/s)"})
    for k_, v_ in plan.dv_continuous.items():
        if v_ > 0:
            ph = [ph_ for ph_ in plan.phases if ph_[0].startswith("final translation" if "final" in k_ else "hold")][0]
            burns_rows.append({"burn_id": "C-" + k_.split()[0], "epoch_utc": utc(t0 + ph[1]),
                               "t_from_window_start_s": round(ph[1], 3), "dv_magnitude_mps": round(v_, 6),
                               "purpose": f"continuous {k_} over {ph[2] - ph[1]:.0f} s (2 n v radial accel.)"})
    df_burns = pd.DataFrame(burns_rows).sort_values("t_from_window_start_s", kind="stable").reset_index(drop=True)
    sun_traj = sun_fn(t0 + plan.t)
    ang_port = np.degrees(np.arccos(np.clip(sun_traj @ plan.strategy.axis, -1, 1)))
    rng = np.linalg.norm(plan.r, axis=1)
    df_traj = pd.DataFrame({"utc": sc.ts.iso(t0 + plan.t), "t_from_window_start_s": np.round(plan.t, 3),
                            "phase": plan.phase_of_t, "x_radial_m": plan.r[:, 0].round(4), "y_alongtrack_m": plan.r[:, 1].round(4),
                            "z_crosstrack_m": plan.r[:, 2].round(4), "vx_mps": plan.v[:, 0].round(6),
                            "vy_mps": plan.v[:, 1].round(6), "vz_mps": plan.v[:, 2].round(6), "range_m": rng.round(4),
                            "sun_angle_on_port_deg": ang_port.round(3),
                            "sun_angle_off_camera_boresight_deg": (180 - ang_port).round(3)})
    files = {
        "maneuver_plan": P.write_csv(df_burns, outdir / "csv" / "maneuver_plan.csv"),
        "trajectory": P.write_csv(df_traj, outdir / "csv" / "approach_trajectory_lvlh_1s.csv"),
        "strategy_trade": P.write_csv(df_tr, outdir / "csv" / "strategy_trade.csv"),
        "timeline_trade": P.write_csv(df_tl, outdir / "csv" / "timeline_trade.csv"),
        "missed_burn_safety": P.write_csv(pd.DataFrame([{k: v for k, v in d.items() if k != "trajectory"} for d in drift]),
                                          outdir / "csv" / "missed_burn_safety.csv"),
        "cam_scan": P.write_csv(pd.DataFrame(cam), outdir / "csv" / "collision_avoidance_scan.csv"),
        "plume": P.write_csv(pd.DataFrame(plume_metrics(plan)), outdir / "csv" / "plume_exposure.csv"),
        "window_robustness": P.write_csv(pd.DataFrame(rob), outdir / "csv" / "window_robustness.csv"),
        "nonlinear_log": P.write_csv(pd.DataFrame(nl_log, columns=["t_s", "x_m", "y_m", "z_m", "vx_mps", "vy_mps", "vz_mps"]),
                                     outdir / "csv" / "nonlinear_j2_relative_trajectory.csv"),
    }
    if make_figures:
        from .viz import rpod_figs as RF
        t_hp = plan.params["T1"] + plan.params["T2"]
        files["fig_trajectory"] = RF.trajectory(plan, outdir / "figures" / "approach_trajectory.png", p2.keep_out_sphere_radius_m,
                                                p2.corridor_half_angle_deg, sun_fn(t0 + t_hp)[0], drift=drift[1:3], r_col=p2.collision_radius_m)
        files["fig_timeline"] = RF.timeline(plan, outdir / "figures" / "part2_timeline.png", t0, window_s,
                                            lambda t: np.degrees(np.arccos(np.clip(sun_fn(t) @ plan.strategy.axis, -1, 1))),
                                            win["station"], utc_fn=utc)
        files["fig_trade"] = RF.trade_small_multiples(df_tr, plans, sun_fn, t0, outdir / "figures" / "strategy_trade_sun.png")
        files["fig_mc"] = RF.dispersions(mc_raw, plan, outdir / "figures" / "monte_carlo.png", p2.corridor_half_angle_deg)
        pv = plans["vbar_plus"]                    # the flown split (whole seconds), not the sweep's continuous optimum
        flown = {"T1": pv.params["T1"], "dv_total": pv.dv_total, "B3": pv.burns[2].mag,
                 "miss_B2": safety_metrics(n, pv, r_ko=p2.collision_radius_m)["miss_if_B2_fails_m"]}
        files["fig_tl"] = RF.timeline_trade(rows, flown, outdir / "figures" / "timeline_trade.png",
                                            b3_max=0.15, miss_min=2 * p2.collision_radius_m)
        files["fig_cam"] = RF.cam_figure(n, plan, cam, outdir / "figures" / "collision_avoidance.png", 0.15, cw)
    summary = {
        "window": {"label": win["label"], "station": win["station"], "start_utc": utc(t0), "end_utc": utc(win["t_end"]),
                   "duration_s": window_s, "band_widening_deg": win["band_widening_deg"]},
        "n_rad_s": n, "initial_state_lvlh": X0.tolist(),
        "sun_lvlh_at_window_start": sun_fn(t0)[0].round(5).tolist(),
        "sun_lvlh_at_window_end": sun_fn(win["t_end"])[0].round(5).tolist(),
        "timeline": timeline, "timeline_choice": chosen,
        "selected_strategy": best, "sun_optimal": extra,
        "burns": df_burns.to_dict(orient="records"),
        "dv_impulsive_mps": plan.dv_impulsive, "dv_continuous_mps": float(sum(plan.dv_continuous.values())),
        "dv_total_mps": plan.dv_total, "t_contact_s": plan.t_contact, "contact_utc": utc(t0 + plan.t_contact),
        "window_margin_s": window_s - plan.t_contact,
        "sun_on_port_critical_phase_deg": [float(ang_port[plan.t >= t_hp].min()), float(ang_port[plan.t >= t_hp].max())],
        "missed_burns": [{k: v for k, v in d.items() if k != "trajectory"} for d in drift],
        "cam": cam, "nonlinear_validation": nl, "monte_carlo": mc_sum,
        "initial_velocity_sensitivity": {"assumed_v0": X0[3:].tolist(), "alt_v0_natural_circular": [0, yd, 0],
                                         "dv_total_alt_mps": alt.dv_total, "B1_alt_mps": alt.burns[0].mag},
        "analytic_worksheet": ws,
        "window_robustness": rob,
        "runtime_s": round(time.time() - t_wall, 1),
        "files": {k: str(v) for k, v in files.items()},
    }
    P.write_json(summary, outdir / "part2_summary.json")
    P.write_json(ws, outdir / "analytic_dv_worksheet.json")
    if verbose:
        print(f"[part2] window {summary['window']['start_utc']} ({win['station']}, {window_s:.1f} s, {win['label']})")
        print(f"[part2] strategy {best}; T1={timeline['T1']} T2={timeline['T2']}; dV={plan.dv_total:.4f} m/s; "
              f"contact {plan.t_contact:.1f} s (margin {window_s - plan.t_contact:.1f} s); {summary['runtime_s']} s")
    return {"plan": plan, "summary": summary, "trade": df_tr, "plans": plans, "scenario": sc, "window": win}


def design_at(sc, cfg, t0, t_end, X0, label):
    """Re-run timeline sizing + strategy trade for another window (robustness)."""
    p2 = cfg.part2
    window_s = t_end - t0
    r_t, v_t = sc.state_gcrs(np.array([t0]))
    n = float(np.linalg.norm(np.cross(r_t[0], v_t[0])) / np.dot(r_t[0], r_t[0]))
    sun_fn = lambda t: sc.sun_dir_lvlh(np.atleast_1d(t))          # noqa: E731
    T_final = p2.hold_point_distance_m / p2.final_approach_speed_mps
    budget = window_s - p2.window_margin_s - p2.hold_duration_s - T_final
    vbar = Strategy("vbar_plus", "+V-bar (from ahead)", np.array([0.0, 1.0, 0.0]))
    best = None
    for T1 in np.arange(0.45 * budget, 0.86 * budget, 2.0):
        pl = build_plan(n, X0, vbar, T1, budget - T1, p2.hold_duration_s, p2.keep_out_sphere_radius_m,
                        p2.hold_point_distance_m, p2.final_approach_speed_mps, p2.final_approach_speed_mps, 0.0)
        sf = safety_metrics(n, pl, r_ko=p2.collision_radius_m)
        if pl.burns[2].mag <= 0.15 and sf["miss_if_B2_fails_m"] >= 2 * p2.collision_radius_m:
            if best is None or pl.dv_total < best[1]:
                best = (T1, pl.dv_total)
    T1 = round(best[0]) if best else round(0.65 * budget)
    tl = dict(T1=T1, T2=round(budget - T1), T_hold=p2.hold_duration_s, rho_entry=p2.keep_out_sphere_radius_m,
              rho_hold=p2.hold_point_distance_m, v_final=p2.final_approach_speed_mps,
              v_contact=p2.final_approach_speed_mps, d_brake=0.0)
    df, plans, _ = run_trade(n, X0, sun_fn, t0, window_s, tl, r_ko=p2.collision_radius_m,
                             cam_excl=p2.camera_sun_exclusion_deg)
    df = score(df)
    el = df[df.eligible].sort_values("rank")
    sel = el.iloc[0]["strategy"] if len(el) else "none eligible"
    row = df[df.strategy == "vbar_plus"].iloc[0]
    return {"case": label, "window_start_utc": sc.ts.iso(t0)[0], "available_s": round(window_s, 1),
            "T1_s": tl["T1"], "T2_s": tl["T2"], "selected_strategy": sel,
            "vbar_dv_total_mps": round(float(row["dv_total_mps"]), 4),
            "vbar_sun_on_port_deg": row["sun_on_port_range_deg"], "vbar_margin_s": round(float(row["window_margin_s"]), 1),
            "eligible_strategies": ",".join(el["strategy"].tolist())}
