"""Docking-window logic (Part I).

Definitions (assignment text -> implementation):
  * Pass P_s,k      : maximal interval with elevation >= 5 deg at station s.
  * Lit             : full sunlight (no umbra and no penumbra) -- ASSUMPTION that
                      any Earth occultation counts as 'eclipse'.
  * Band            : 30 deg <= angle(port axis, target->Sun) <= 60 deg.
  * Feasible set    : F = P_s,k  n  Lit  n  Band   (per pass, per station)
  * Docking window  : every maximal [a, b] in F with b - a >= D (D = 9 min),
                      inside a single pass of duration >= D. The window of
                      admissible START epochs is [a, b - D]: an approach that
                      starts at t0 in [a, b - D] keeps all constraints for D.
Near-miss diagnostics report, for every qualifying pass, the longest feasible
sub-interval, the deficit to D and which constraint bounds it at each end.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
from scipy.optimize import minimize_scalar

from .events import intervals as iv
from .events.finder import find_positive_intervals


@dataclass
class ConstraintSet:
    span: tuple
    passes: dict          # station -> (N,2)
    lit: np.ndarray       # full sun (or not-umbra, per config)
    umbra: np.ndarray
    band: np.ndarray      # sun angle within [lo, hi]
    sun_ge_min: np.ndarray
    sun_le_max: np.ndarray
    sun_min_deg: float
    sun_max_deg: float
    min_el_deg: float


def compute_constraints(sc, sun_min=None, sun_max=None, min_el=None, step=None, tol=None, span=None) -> ConstraintSet:
    c = sc.cfg.constraints
    lo = c.sun_angle_min_deg if sun_min is None else sun_min
    hi = c.sun_angle_max_deg if sun_max is None else sun_max
    el0 = c.min_elevation_deg if min_el is None else min_el
    step = step or sc.cfg.numerics.grid_step_s
    tol = tol or sc.cfg.numerics.root_tol_s
    t0, t1 = span or (sc.t_start, sc.t_end)
    g = sc.constraint_functions(lo, hi, el0)
    fi = lambda fn: find_positive_intervals(fn, t0, t1, step, tol)  # noqa: E731
    passes = {s.name: fi(g[f"vis::{s.name}"]) for s in sc.stations}
    ge_min = fi(g["sun_angle_ge_min"])
    le_max = fi(g["sun_angle_le_max"])
    lit = fi(sc.lighting_function())
    umbra = fi(sc.g_umbra)
    return ConstraintSet((t0, t1), passes, lit, umbra, iv.intersect(ge_min, le_max), ge_min, le_max, lo, hi, el0)


@dataclass
class DockingWindow:
    station: str
    pass_aos_s: float
    pass_los_s: float
    feasible_start_s: float
    feasible_end_s: float
    latest_start_s: float
    slack_s: float


def docking_windows(cs: ConstraintSet, duration_s: float):
    out = []
    for name, P in cs.passes.items():
        for aos, los in iv.filter_min_duration(P, duration_s):
            F = iv.intersect_all([[aos, los]], cs.lit, cs.band)
            for a, b in iv.filter_min_duration(F, duration_s):
                out.append(DockingWindow(name, aos, los, a, b, b - duration_s, (b - a) - duration_s))
    out.sort(key=lambda w: w.feasible_start_s)
    return out


@dataclass
class PassDiagnostic:
    station: str
    aos_s: float
    los_s: float
    duration_s: float
    tca_s: float
    max_elevation_deg: float
    qualifies_duration: bool
    lit_fraction: float
    band_fraction: float
    sun_angle_min_deg: float
    sun_angle_max_deg: float
    best_start_s: float
    best_end_s: float
    best_feasible_s: float
    deficit_s: float
    limited_start_by: str
    limited_end_by: str
    is_window: bool


def _bound_reason(t, cs, edge, pass_iv):
    """Which constraint switches at time t (the edge of a feasible sub-interval)?"""
    eps = 1e-3
    if abs(t - pass_iv[edge]) < eps:
        return "pass AOS" if edge == 0 else "pass LOS"
    probe = t - 0.05 if edge == 0 else t + 0.05
    reasons = []
    if iv.containing(cs.lit, probe) is None:
        reasons.append("eclipse")
    if iv.containing(cs.sun_ge_min, probe) is None:
        reasons.append(f"sun angle < {cs.sun_min_deg:g} deg")
    if iv.containing(cs.sun_le_max, probe) is None:
        reasons.append(f"sun angle > {cs.sun_max_deg:g} deg")
    return " & ".join(reasons) or "span edge"


def pass_diagnostics(sc, cs: ConstraintSet, duration_s: float, only_qualifying=False):
    rows = []
    for k, st in enumerate(sc.stations):
        for aos, los in cs.passes[st.name]:
            dur = los - aos
            if only_qualifying and dur < duration_s:
                continue
            res = minimize_scalar(lambda s: -float(sc.elevation_deg(k, np.array([s]))[0]),
                                  bounds=(aos, los), method="bounded", options={"xatol": 1e-3})
            tca, elmax = float(res.x), -float(res.fun)
            F = iv.intersect_all([[aos, los]], cs.lit, cs.band)
            tt = np.linspace(aos, los, max(int(dur / 2.0), 3))
            ang = sc.sun_angle_deg(tt)
            if len(F):
                j = int(np.argmax(iv.durations(F)))
                a, b = F[j]
                ls, le = _bound_reason(a, cs, 0, (aos, los)), _bound_reason(b, cs, 1, (aos, los))
            else:
                a = b = np.nan
                ls = le = "no feasible instant"
            best = 0.0 if np.isnan(a) else b - a
            rows.append(PassDiagnostic(
                st.name, aos, los, dur, tca, elmax, dur >= duration_s,
                iv.total(iv.intersect([[aos, los]], cs.lit)) / dur,
                iv.total(iv.intersect([[aos, los]], cs.band)) / dur,
                float(ang.min()), float(ang.max()), a, b, best, max(0.0, duration_s - best), ls, le,
                (dur >= duration_s) and (best >= duration_s - 1e-9)))
    rows.sort(key=lambda r: r.aos_s)
    return rows


def as_dicts(rows):
    return [asdict(r) for r in rows]
