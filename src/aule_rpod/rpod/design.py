"""Approach-strategy definitions, maneuver-plan builder and strategy evaluation.

Waypoint architecture (same structure as crewed-vehicle ISS approaches: far
waypoint -> corridor entry -> hold point with go/no-go -> final translation):

    P0  (initial point, t0 = window start)
      | Leg 1: CW two-impulse transfer, time T1           Burn 1 (departure)
    WP1 (corridor entry on the approach axis, rho1)
      | Leg 2: CW two-impulse transfer, time T2           Burn 2 (= leg-1 arrival + leg-2 departure)
    HP  (hold point, 10 m in front of the port)           Burn 3 (stop)
      | Hold T_hold (go/no-go); station keeping if HP is not an equilibrium
      | Final translation at v_final along -axis           Burn 4 (start) + continuous compensation
      | braking to v_contact at d_brake                   Burn 5
    Contact (port at the LVLH origin)

The target re-orients its docking port along the approach axis (assignment
Part II premise) and holds that LVLH-fixed attitude during the approach.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import cw

EARTH_ANGULAR_RADIUS_DEG = None   # set from altitude at run time


@dataclass
class Strategy:
    key: str
    name: str
    axis: np.ndarray            # unit vector: target -> chaser along the final corridor (= port normal)


@dataclass
class Burn:
    id: str
    t_s: float                  # from window start
    dv: np.ndarray              # LVLH impulse [m/s]
    r: np.ndarray               # position at burn [m]
    v_before: np.ndarray
    v_after: np.ndarray
    purpose: str

    @property
    def mag(self):
        return float(np.linalg.norm(self.dv))


@dataclass
class Plan:
    strategy: Strategy
    params: dict
    burns: list
    phases: list                # (name, t_start, t_end)
    t: np.ndarray
    r: np.ndarray
    v: np.ndarray
    phase_of_t: np.ndarray
    dv_continuous: dict = field(default_factory=dict)

    @property
    def dv_impulsive(self):
        return float(sum(b.mag for b in self.burns))

    @property
    def dv_total(self):
        return self.dv_impulsive + float(sum(self.dv_continuous.values()))

    @property
    def t_contact(self):
        return float(self.phases[-1][2])


def build_plan(n, X0, strat: Strategy, T1, T2, T_hold, rho_entry, rho_hold, v_final, v_contact, d_brake,
               dt=1.0):
    a = strat.axis / np.linalg.norm(strat.axis)
    r0, v0 = np.asarray(X0[:3], float), np.asarray(X0[3:], float)
    WP1, HP = rho_entry * a, rho_hold * a
    # Leg 1 and leg 2 (CW two-impulse)
    dv1, _, v0p, v1m = cw.two_impulse(n, r0, v0, WP1, T1)
    _, dv3, v1p, v2m = cw.two_impulse(n, WP1, v1m, HP, T2)
    dv2 = v1p - v1m
    burns = [Burn("B1", 0.0, dv1, r0, v0, v0p, "departure from initial point toward corridor entry"),
             Burn("B2", T1, dv2, WP1, v1m, v1p, "corridor entry: brake + start slow corridor leg"),
             Burn("B3", T1 + T2, dv3, HP, v2m, np.zeros(3), "stop at 10 m hold point")]
    t_list, r_list, v_list, ph = [], [], [], []

    def add_free(t_start, T, r_s, v_s, label):
        tt = np.arange(0.0, T, dt)
        X = cw.propagate(n, np.r_[r_s, v_s], tt)
        t_list.append(t_start + tt)
        r_list.append(X[:, :3])
        v_list.append(X[:, 3:])
        ph.extend([label] * len(tt))

    add_free(0.0, T1, r0, v0p, "leg 1 transfer")
    add_free(T1, T2, WP1, v1p, "leg 2 corridor")
    # hold
    t_h0 = T1 + T2
    sk_acc = -cw.natural_accel(n, HP, np.zeros(3))
    dv_sk = float(np.linalg.norm(sk_acc) * T_hold)
    th = np.arange(0.0, T_hold, dt)
    t_list.append(t_h0 + th)
    r_list.append(np.repeat(HP[None], len(th), 0))
    v_list.append(np.zeros((len(th), 3)))
    ph.extend(["hold (go/no-go)"] * len(th))
    # final translation (forced straight line) -> optional braking -> contact
    t_f0 = t_h0 + T_hold
    braking = (d_brake > 0) and (v_contact < v_final)
    burns.append(Burn("B4", t_f0, -a * v_final, HP, np.zeros(3), -a * v_final, "start final translation"))
    segs = []
    if braking:
        brake_pt = d_brake * a
        f1 = cw.forced_line(n, HP, brake_pt, v_final, dt=0.5)
        f2 = cw.forced_line(n, brake_pt, np.zeros(3), v_contact, dt=0.5)
        burns.append(Burn("B5", t_f0 + f1["T"], a * (v_final - v_contact), brake_pt, -a * v_final, -a * v_contact,
                          "brake to contact speed"))
        segs = [(f1, t_f0, "final translation"), (f2, t_f0 + f1["T"], "final braking/contact")]
    else:
        f1 = cw.forced_line(n, HP, np.zeros(3), v_final, dt=0.5)
        segs = [(f1, t_f0, "final translation")]
    for seg, t_s, lab in segs:
        k = seg["t"] < seg["T"] - 1e-9
        t_list.append(t_s + seg["t"][k])
        r_list.append(seg["r"][k])
        v_list.append(seg["v"][k])
        ph.extend([lab] * int(k.sum()))
    t_contact = t_f0 + sum(s_[0]["T"] for s_ in segs)
    v_c = v_contact if braking else v_final
    t_list.append(np.array([t_contact]))
    r_list.append(np.zeros((1, 3)))
    v_list.append((-a * v_c)[None])
    ph.append("contact")
    phases = [("leg 1 transfer", 0.0, T1), ("leg 2 corridor", T1, T1 + T2), ("hold (go/no-go)", t_h0, t_f0)]
    phases += [(lab, t_s, t_s + seg["T"]) for seg, t_s, lab in segs]
    params = dict(T1=T1, T2=T2, T_hold=T_hold, rho_entry=rho_entry, rho_hold=rho_hold, v_final=v_final,
                  v_contact=v_c, d_brake=d_brake if braking else 0.0, WP1=WP1, HP=HP)
    cont = {"hold station-keeping": dv_sk, "final translation compensation": segs[0][0]["dv_continuous"]}
    if braking:
        cont["braking-leg compensation"] = segs[1][0]["dv_continuous"]
    return Plan(strat, params, burns, phases, np.concatenate(t_list), np.vstack(r_list), np.vstack(v_list),
                np.array(ph), cont)


# ------------------------------------------------------------------- metrics
def illumination_metrics(plan: Plan, sun_lvlh_fn, t0_abs, band=(30.0, 60.0), cam_excl=30.0, earth_ang_deg=67.0,
                         cam_half_fov=15.0):
    """Sun on the port & camera geometry from HP arrival to contact (critical phase)."""
    a = plan.strategy.axis
    t_hp = plan.params["T1"] + plan.params["T2"]
    tt = np.linspace(t_hp, plan.t_contact, 200)
    s = sun_lvlh_fn(t0_abs + tt)
    port = np.degrees(np.arccos(np.clip(s @ a, -1, 1)))
    off_bore = 180.0 - port                                   # camera boresight = -a
    nadir_to_bore = np.degrees(np.arccos(np.clip(-a @ np.array([-1.0, 0, 0]), -1, 1)))
    earth_in_fov = nadir_to_bore < earth_ang_deg + cam_half_fov
    tt_all = np.linspace(0.0, plan.t_contact, 300)
    port_all = np.degrees(np.arccos(np.clip(sun_lvlh_fn(t0_abs + tt_all) @ a, -1, 1)))
    return {"sun_on_port_min_deg": float(port.min()), "sun_on_port_max_deg": float(port.max()),
            "sun_on_port_in_band_fraction": float(np.mean((port >= band[0]) & (port <= band[1]))),
            "sun_on_port_whole_window_range_deg": [float(port_all.min()), float(port_all.max())],
            "sun_off_camera_boresight_min_deg": float(off_bore.min()),
            "camera_sun_exclusion_ok": bool(off_bore.min() >= cam_excl),
            "port_face_sunlit": bool(port.max() < 90.0),
            "camera_background": "Earth" if earth_in_fov else "deep space",
            "boresight_angle_from_nadir_deg": float(nadir_to_bore)}


def safety_metrics(n, plan: Plan, horizon_orbits=2.0, r_ko=5.0):
    """Missed-burn free-drift analysis: minimum range over `horizon_orbits` periods."""
    T = 2 * np.pi / n * horizon_orbits
    out = {}
    for b in plan.burns[:3]:
        rmin, tmin, _ = cw.min_range_free_drift(n, np.r_[b.r, b.v_before], T, dt=1.0)
        out[f"miss_if_{b.id}_fails_m"] = rmin
        out[f"t_closest_if_{b.id}_fails_s"] = tmin
    # path clearance of the planned transfer legs (excluding corridor/final phases)
    k = np.isin(plan.phase_of_t, ["leg 1 transfer"])
    out["leg1_min_range_m"] = float(np.linalg.norm(plan.r[k], axis=1).min())
    out["passively_safe_far_field"] = bool(out["miss_if_B1_fails_m"] > r_ko and out["miss_if_B2_fails_m"] > r_ko)
    return out


def plume_metrics(plan: Plan, half_angle_deg=30.0):
    """For every impulse, the thruster exhaust points along -dv. Report the angle
    between exhaust and the line of sight to the target and flag exposure when the
    target lies inside the plume cone."""
    rows = []
    for b in plan.burns:
        if b.mag < 1e-9:
            continue
        ex = -b.dv / b.mag
        rng = np.linalg.norm(b.r)
        los = -b.r / rng if rng > 1e-9 else np.zeros(3)
        ang = float(np.degrees(np.arccos(np.clip(ex @ los, -1, 1)))) if rng > 1e-9 else 0.0
        rows.append({"burn": b.id, "range_m": rng, "exhaust_to_target_angle_deg": ang,
                     "target_in_plume_cone": bool(ang < half_angle_deg and rng > 0.5), "dv_mps": b.mag})
    return rows
