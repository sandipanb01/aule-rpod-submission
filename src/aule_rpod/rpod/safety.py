"""Trajectory safety: missed-burn drift, collision-avoidance manoeuvre (CAM) design.

Passive safety (Fehse 2003, Ch. 4): after any single missed burn the free-drift
trajectory must stay outside the keep-out sphere (KOS) for >= 2 orbits.
Active safety: inside the corridor, where passive safety is impossible by
construction (the vehicle is meant to reach the target), a single pre-planned
CAM impulse must, from ANY point of the corridor/hold/final phases, (i) stop
the closing motion immediately and (ii) never bring the chaser back inside the
KOS for >= 2 orbits. Candidate CAM directions are scanned and the most robust
one is selected.
"""
from __future__ import annotations

import numpy as np

from . import cw

CAM_DIRECTIONS = {
    "+X (radial up)": [1, 0, 0], "-X (radial down)": [-1, 0, 0],
    "+Y (posigrade)": [0, 1, 0], "-Y (retrograde)": [0, -1, 0],
    "+Z (normal)": [0, 0, 1], "-Z (anti-normal)": [0, 0, -1],
    "away from target": None,
}


def cam_scan(n, plan, dv_cam=0.15, r_ko=5.0, horizon_orbits=2.0, dt_fail=2.0, phases=None):
    """Breakout (collision-avoidance) burn of dv_cam fired from every dt_fail-th point of the
    close-range phases, in each candidate direction, then free drift for horizon_orbits.
    Per direction, worst case over all trigger points:
      worst_additional_closing_m          min over the drift of (range - range at trigger); < 0 closes
      closest_return_10min_to_2orbits_m   closest approach between 10 min and the horizon
    A direction is safe when it never closes by more than 5 cm and never returns inside r_ko."""
    phases = phases or ["leg 2 corridor", "hold (go/no-go)", "final translation", "final braking/contact"]
    T = 2 * np.pi / n * horizon_orbits
    tt = np.arange(0.0, T, 2.0)
    Phi = cw.stm(n, tt)
    sel = np.where(np.isin(plan.phase_of_t, phases))[0][:: max(1, int(dt_fail))]
    rows = []
    for name, d in CAM_DIRECTIONS.items():
        worst_close, worst_long = np.inf, np.inf
        for i in sel:
            r, v = plan.r[i], plan.v[i]
            rng0 = np.linalg.norm(r)
            if rng0 < 0.3:
                continue
            u = (r / rng0) if d is None else np.asarray(d, float)
            X = np.einsum("tij,j->ti", Phi, np.r_[r, v + dv_cam * u])
            rng = np.linalg.norm(X[:, :3], axis=1)
            worst_close = min(worst_close, rng.min() - rng0)          # <0 => it closes further
            late = rng[tt > 600.0]
            worst_long = min(worst_long, late.min())
        rows.append({"cam_direction": name, "dv_cam_mps": dv_cam,
                     "worst_additional_closing_m": float(worst_close),
                     "closest_return_10min_to_2orbits_m": float(worst_long),
                     "safe": bool(worst_close > -0.05 and worst_long > r_ko)})
    return rows


def missed_burn_drift(n, plan, horizon_orbits=2.0):
    T = 2 * np.pi / n * horizon_orbits
    out = []
    for b in plan.burns:
        rmin, tmin, X = cw.min_range_free_drift(n, np.r_[b.r, b.v_before], T, dt=1.0)
        out.append({"failed_burn": b.id, "purpose": b.purpose, "t_burn_s": b.t_s,
                    "range_at_burn_m": float(np.linalg.norm(b.r)), "min_range_m": rmin,
                    "time_of_min_range_s": tmin, "trajectory": X})
    return out
