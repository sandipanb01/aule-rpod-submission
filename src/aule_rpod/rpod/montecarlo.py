"""Closed-loop Monte Carlo of the approach (navigation + burn-execution errors).

Guidance re-targets every impulse from the NAVIGATED state, as flight software
would; the TRUE state is propagated with the CW model. Each transfer leg carries
a mid-course correction (MCC) at half time. The 10-m final translation is flown
by a discrete closed-loop controller: every ``dt_ctrl`` seconds the commanded
velocity is ``v_final`` toward the port along the corridor axis plus a lateral
correction ``-lat / tau_lat``; a correction fires when the commanded change
exceeds the thruster deadband. The CW (Coriolis) drift during the translation
is therefore compensated by the controller itself, not by an assumed 2nL term.

Error models (1 sigma, per axis):
  navigation position  max(floor, frac * range); velocity sig_v. Independent
                       draws at each impulsive burn; first-order Gauss-Markov
                       processes (correlation time tau_nav) during the
                       final translation, where the controller samples them
                       every dt_ctrl.
  execution            magnitude scale factor sig_mag; pointing error sig_pt
                       about a random axis.

All samples are propagated together (vectorised), so 5,000 runs take seconds.
Outputs: hold-point dispersions, corridor and keep-out compliance, contact
conditions in the terms of the IDSS initial-contact table (closing rate,
lateral rate, lateral misalignment), and the delta-v distribution.
"""
from __future__ import annotations

import numpy as np

from . import cw

# IDSS IDD Rev. E, Table 3.3.1.1-2 (initial contact conditions, translational subset)
IDSS_CLOSING_RATE_MPS = (0.05, 0.10)
IDSS_LATERAL_RATE_MPS = 0.04
IDSS_LATERAL_MISALIGNMENT_M = 0.10


def _phi(n, t):
    return cw.stm(n, np.asarray(t, dtype=float))


def _propagate(n, X, t):
    """X (N,6) propagated by a single time t -> (N,6)."""
    return X @ _phi(n, t).T


def _trajectory(n, X, tt):
    """X (N,6) on a grid tt (K,) -> (N,K,6)."""
    return np.einsum("kij,nj->nki", _phi(n, tt), X)


def _target(n, r0, rf, T):
    """Batch two-impulse targeting: departure velocity that reaches rf in T."""
    rr, rv, _, _ = cw.stm_blocks(n, np.array(float(T)))
    return np.linalg.solve(rv, (np.asarray(rf, float)[None, :] - r0 @ rr.T).T).T


def _perturb(dv, rng, sig_mag, sig_pt):
    """Apply magnitude and pointing errors to a batch of impulses (N,3)."""
    m = np.linalg.norm(dv, axis=1)
    out = np.zeros_like(dv)
    ok = m > 1e-12
    if not ok.any():
        return out
    u = dv[ok] / m[ok, None]
    ref = np.where(np.abs(u[:, :1]) < 0.9, np.array([[1.0, 0, 0]]), np.array([[0, 1.0, 0]]))
    p = np.cross(u, ref)
    p /= np.linalg.norm(p, axis=1)[:, None]
    q = np.cross(u, p)
    th = rng.normal(0.0, sig_pt, ok.sum())
    ph = rng.uniform(0.0, 2 * np.pi, ok.sum())
    u2 = np.cos(th)[:, None] * u + np.sin(th)[:, None] * (np.cos(ph)[:, None] * p + np.sin(ph)[:, None] * q)
    out[ok] = (m[ok] * (1.0 + rng.normal(0.0, sig_mag, ok.sum())))[:, None] * u2
    return out


def _nav(X, rng, sig_r_frac, sig_r_floor, sig_v):
    N = len(X)
    sr = np.maximum(sig_r_floor, sig_r_frac * np.linalg.norm(X[:, :3], axis=1))
    return X + np.hstack([rng.normal(0.0, 1.0, (N, 3)) * sr[:, None], rng.normal(0.0, sig_v, (N, 3))])


def _off_axis_deg(r, a):
    along = r @ a
    lat = np.linalg.norm(r - along[..., None] * a, axis=-1)
    return np.degrees(np.arctan2(lat, along))


def run(n, X0, plan, samples, seed, sig_r_frac, sig_r_floor, sig_v, sig_mag, sig_pt_deg, r_collision=5.0,
        r_kos=25.0, corridor_deg=10.0, dt_ctrl=2.0, tau_lat=20.0, v_lat_max=0.02, deadband=0.002, tau_nav=30.0):
    rng = np.random.default_rng(seed)
    p = plan.params
    a = np.asarray(plan.strategy.axis, float)
    a = a / np.linalg.norm(a)
    T1, T2, Th = float(p["T1"]), float(p["T2"]), float(p["T_hold"])
    WP1, HP = np.asarray(p["WP1"], float), np.asarray(p["HP"], float)
    v_fin = float(p["v_final"])
    sig_pt = np.radians(sig_pt_deg)
    N = int(samples)
    X = np.tile(np.asarray(X0, float), (N, 1))
    dv_total = np.zeros(N)
    out = {}

    def impulse(X, dv_cmd):
        e = _perturb(dv_cmd, rng, sig_mag, sig_pt)
        X = X.copy()
        X[:, 3:] += e
        return X, np.linalg.norm(e, axis=1)

    def transfer(X, rf, T, tag):
        """Leg with a departure burn and a mid-course correction at T/2, both
        re-targeted from the navigated state. Returns the true trajectory too."""
        nonlocal dv_total
        Xn = _nav(X, rng, sig_r_frac, sig_r_floor, sig_v)
        X, m = impulse(X, _target(n, Xn[:, :3], rf, T) - Xn[:, 3:])
        dv_total += m
        out[f"{tag}_departure_dv"] = m
        tt = np.linspace(0.0, 0.5 * T, 40)
        tr_a = _trajectory(n, X, tt)
        X = tr_a[:, -1]
        Xn = _nav(X, rng, sig_r_frac, sig_r_floor, sig_v)
        X, m = impulse(X, _target(n, Xn[:, :3], rf, 0.5 * T) - Xn[:, 3:])
        dv_total += m
        out[f"{tag}_mcc_dv"] = m
        tr_b = _trajectory(n, X, tt)
        return tr_b[:, -1], np.concatenate([tr_a, tr_b[:, 1:]], axis=1)

    # ------------------------------------------------ legs 1 and 2 (impulsive, re-targeted)
    X, tr1 = transfer(X, WP1, T1, "leg1")
    out["min_range_leg1"] = np.linalg.norm(tr1[..., :3], axis=-1).min(axis=1)
    X, tr2 = transfer(X, HP, T2, "leg2")
    tr12 = np.concatenate([tr1, tr2], axis=1)
    rng12 = np.linalg.norm(tr12[..., :3], axis=-1)
    off12 = _off_axis_deg(tr12[..., :3], a)
    out["min_range_legs"] = rng12.min(axis=1)
    out["kos_outside_corridor"] = ((rng12 < r_kos - 1e-9) & (off12 > corridor_deg)).any(axis=1)
    out["hp_pos_err"] = np.linalg.norm(X[:, :3] - HP, axis=1)
    # ------------------------------------------------ B3: stop at the hold point (null navigated velocity)
    Xn = _nav(X, rng, sig_r_frac, sig_r_floor, sig_v)
    X, m = impulse(X, -Xn[:, 3:])
    dv_total += m
    out["hp_vel_err"] = np.linalg.norm(X[:, 3:], axis=1)
    # ------------------------------------------------ hold: free drift (V-bar equilibrium) with residuals
    X = _propagate(n, X, Th)
    out["hold_end_offset"] = np.linalg.norm(X[:, :3] - HP, axis=1)
    out["corridor_ok_hold_end"] = _off_axis_deg(X[:, :3], a) <= corridor_deg
    # ------------------------------------------------ closed-loop final translation (B4 + corrections)
    phi = np.exp(-dt_ctrl / tau_nav)
    gm_r = rng.normal(0.0, 1.0, (N, 3))
    gm_v = rng.normal(0.0, 1.0, (N, 3))
    P = _phi(n, dt_ctrl)
    active = np.ones(N, bool)
    dv_ft = np.zeros(N)
    n_burns = np.zeros(N, int)
    c_lat = np.full(N, np.nan)
    c_close = np.full(N, np.nan)
    c_latrate = np.full(N, np.nan)
    c_time = np.full(N, np.nan)
    corridor_ft = np.ones(N, bool)
    k_max = int(np.ceil(2.0 * np.linalg.norm(HP) / v_fin / dt_ctrl))
    for k in range(k_max):
        if not active.any():
            break
        r, v = X[:, :3], X[:, 3:]
        sr = np.maximum(sig_r_floor, sig_r_frac * np.linalg.norm(r, axis=1))
        rn = r + gm_r * sr[:, None]
        vn = v + gm_v * sig_v
        lat_n = rn - (rn @ a)[:, None] * a
        lat_cmd = -lat_n / tau_lat
        s_ = np.linalg.norm(lat_cmd, axis=1)
        lat_cmd *= np.minimum(1.0, v_lat_max / np.maximum(s_, 1e-15))[:, None]
        dv_cmd = (-v_fin * a + lat_cmd) - vn
        fire = active & (np.linalg.norm(dv_cmd, axis=1) > deadband)
        if fire.any():
            e = _perturb(dv_cmd[fire], rng, sig_mag, sig_pt)
            X[fire, 3:] += e
            m = np.linalg.norm(e, axis=1)
            dv_ft[fire] += m
            n_burns[fire] += 1
        X1 = X @ P.T
        s0, s1 = X[:, :3] @ a, X1[:, :3] @ a
        lat1 = np.linalg.norm(X1[:, :3] - s1[:, None] * a, axis=1)
        corridor_ft &= ~(active & (s1 > 0) & (lat1 > np.maximum(s1 * np.tan(np.radians(corridor_deg)),
                                                                   IDSS_LATERAL_MISALIGNMENT_M)))
        hit = active & (s1 <= 0.0)
        if hit.any():
            f = (s0[hit] / (s0[hit] - s1[hit]))[:, None]
            rc = X[hit, :3] + f * (X1[hit, :3] - X[hit, :3])
            vc = X[hit, 3:] + f * (X1[hit, 3:] - X[hit, 3:])
            c_lat[hit] = np.linalg.norm(rc - (rc @ a)[:, None] * a, axis=1)
            c_close[hit] = -(vc @ a)
            c_latrate[hit] = np.linalg.norm(vc - (vc @ a)[:, None] * a, axis=1)
            c_time[hit] = (k + f[:, 0]) * dt_ctrl
            active[hit] = False
        X = np.where(active[:, None], X1, X)
        gm_r = phi * gm_r + np.sqrt(1 - phi * phi) * rng.normal(0.0, 1.0, (N, 3))
        gm_v = phi * gm_v + np.sqrt(1 - phi * phi) * rng.normal(0.0, 1.0, (N, 3))
    dv_total += dv_ft
    out.update(dv_final_translation=dv_ft, final_translation_burns=n_burns, lateral_at_contact=c_lat,
               contact_speed=c_close, contact_lateral_rate=c_latrate, final_translation_time=c_time,
               corridor_ok_final=corridor_ft, dv_total=dv_total)
    out = {k: np.asarray(v) for k, v in out.items()}
    contacted = np.isfinite(c_lat)
    idss = (contacted & (c_close >= IDSS_CLOSING_RATE_MPS[0]) & (c_close <= IDSS_CLOSING_RATE_MPS[1])
            & (c_latrate <= IDSS_LATERAL_RATE_MPS) & (c_lat <= IDSS_LATERAL_MISALIGNMENT_M))
    out["idss_ok"] = idss
    pc = lambda x, q: float(np.nanpercentile(x, q))                     # noqa: E731
    summ = {
        "samples": N,
        "model": {"mcc_at_half_leg": True, "control_period_s": dt_ctrl, "lateral_time_constant_s": tau_lat,
                  "lateral_speed_limit_mps": v_lat_max, "deadband_mps": deadband, "nav_correlation_time_s": tau_nav},
        "hp_position_error_m_mean_p99": [float(out["hp_pos_err"].mean()), pc(out["hp_pos_err"], 99)],
        "hp_residual_velocity_mm_s_p99": 1000 * pc(out["hp_vel_err"], 99),
        "hold_end_offset_m_p99": pc(out["hold_end_offset"], 99),
        "corridor_compliance_at_hold_end": float(out["corridor_ok_hold_end"].mean()),
        "min_range_leg1_m_p01": pc(out["min_range_leg1"], 1),
        "collision_radius_entry_probability": float(np.mean(out["min_range_legs"] < r_collision)),
        "keep_out_entry_outside_corridor_probability": float(out["kos_outside_corridor"].mean()),
        "leg1_mcc_dv_mm_s_mean_p99": [1000 * float(out["leg1_mcc_dv"].mean()), 1000 * pc(out["leg1_mcc_dv"], 99)],
        "leg2_mcc_dv_mm_s_mean_p99": [1000 * float(out["leg2_mcc_dv"].mean()), 1000 * pc(out["leg2_mcc_dv"], 99)],
        "final_translation_dv_mps_mean_p99": [float(dv_ft.mean()), pc(dv_ft, 99)],
        "final_translation_burns_median": float(np.median(n_burns)),
        "final_translation_corridor_compliance": float(corridor_ft.mean()),
        "contact_fraction": float(contacted.mean()),
        "contact_closing_speed_mps_p01_p99": [pc(c_close, 1), pc(c_close, 99)],
        "contact_lateral_rate_mm_s_p99": 1000 * pc(c_latrate, 99),
        "contact_lateral_misalignment_m_mean_p99": [float(np.nanmean(c_lat)), pc(c_lat, 99)],
        "final_translation_duration_s_p01_p99": [pc(c_time, 1), pc(c_time, 99)],
        "idss_translational_envelope_fraction": float(idss.mean()),
        "dv_total_mps_mean_p99": [float(dv_total.mean()), pc(dv_total, 99)],
    }
    return summ, out
