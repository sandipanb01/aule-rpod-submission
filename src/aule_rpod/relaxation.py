"""Minimum-relaxation analysis: how far is each pass from being a docking window?

For every ground-station pass of duration >= D, compute on a dense (0.5 s) grid
the smallest relaxation that would make it a valid D-minute window, holding the
other constraints fixed:

  delta_band   : uniform widening of the illumination band to [lo-d, hi+d]
  delta_lo     : lowering of the lower bound only (upper bound kept)
  delta_hi     : raising of the upper bound only
  deficit_s    : D minus the longest strictly feasible sub-interval

computed with sliding-window min/max over the sampled sun angle (exact to the
grid; the grid is 0.5 s, i.e. < 0.01 deg of sun-angle change).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import maximum_filter1d, minimum_filter1d

from .events import intervals as iv


@dataclass
class PassRelaxation:
    station: str
    aos_s: float
    los_s: float
    delta_band_deg: float
    best_start_band_s: float
    delta_lo_deg: float
    best_start_lo_s: float
    delta_hi_deg: float
    sun_min_in_best_deg: float
    sun_max_in_best_deg: float


def pass_relaxations(sc, cs, duration_s: float, dt: float = 0.5):
    lo, hi = cs.sun_min_deg, cs.sun_max_deg
    out = []
    for name, P in cs.passes.items():
        for aos, los in iv.filter_min_duration(P, duration_s):
            t = np.arange(aos, los + 1e-9, dt)
            ang = sc.sun_angle_deg(t)
            lit = np.array([iv.containing(cs.lit, x) is not None for x in t])
            w = int(round(duration_s / dt)) + 1
            if len(t) < w:
                continue
            n_starts = len(t) - w + 1
            c = w // 2
            mn = minimum_filter1d(ang, w, mode="nearest")[c:c + n_starts]
            mx = maximum_filter1d(ang, w, mode="nearest")[c:c + n_starts]
            dark = maximum_filter1d((~lit).astype(float), w, mode="nearest")[c:c + n_starts] > 0
            need_band = np.maximum.reduce([lo - mn, mx - hi, np.zeros_like(mn)])
            need_band[dark] = np.inf
            need_lo = np.where(mx <= hi, np.maximum(lo - mn, 0.0), np.inf)
            need_lo[dark] = np.inf
            need_hi = np.where(mn >= lo, np.maximum(mx - hi, 0.0), np.inf)
            need_hi[dark] = np.inf
            j = int(np.argmin(need_band))
            k = int(np.argmin(need_lo))
            out.append(PassRelaxation(name, aos, los, float(need_band[j]), float(t[j]), float(need_lo[k]),
                                      float(t[k]), float(np.min(need_hi)), float(mn[j]), float(mx[j])))
    out.sort(key=lambda r: r.aos_s)
    return out
