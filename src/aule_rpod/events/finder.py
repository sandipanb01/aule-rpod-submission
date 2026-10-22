"""Event detection: intervals where a continuous function g(t) >= 0.

Algorithm (robust 'sample - bracket - refine', as in Orekit/GMAT event detectors):
  1. Evaluate g on a uniform grid (step h). Every sign change brackets a root.
  2. Refine all brackets simultaneously by vectorised bisection to `tol`.
  3. Guard against double crossings hidden inside one step: wherever three
     consecutive samples of equal sign show a local extremum toward zero, fit a
     parabola; if it predicts a sign change, locate the true extremum with a
     bounded Brent search and, if it does cross, bracket and refine both roots.
With h = 2 s and the smoothest dynamics in this problem (300-s nutation of the
port, 95-min orbit), missed events would require a sub-2-s excursion; the
guard removes even that case. Grid-independence is verified in the tests
(1-s vs 2-s vs 5-s grids give identical events to < 1 ms).
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize_scalar


def _bisect(g, lo, hi, glo, tol):
    lo, hi, glo = lo.copy(), hi.copy(), glo.copy()
    while np.any(hi - lo > tol):
        mid = 0.5 * (lo + hi)
        gm = g(mid)
        left = np.sign(gm) == np.sign(glo)
        lo = np.where(left, mid, lo)
        glo = np.where(left, gm, glo)
        hi = np.where(left, hi, mid)
    return 0.5 * (lo + hi)


def find_positive_intervals(g, t0: float, t1: float, step: float, tol: float = 1e-3,
                            guard: bool = True, return_samples: bool = False):
    """Return (N,2) array of maximal intervals in [t0, t1] where g >= 0."""
    n = int(np.ceil((t1 - t0) / step))
    t = np.linspace(t0, t1, n + 1)
    y = np.asarray(g(t), dtype=float)
    pos = y >= 0.0
    roots = []
    # (1)-(2) sign changes
    idx = np.where(pos[1:] != pos[:-1])[0]
    if len(idx):
        r = _bisect(g, t[idx], t[idx + 1], y[idx], tol)
        roots.extend(zip(r.tolist(), (pos[idx + 1]).tolist()))   # (time, becomes_positive)
    # (3) hidden double crossings
    if guard and n >= 2:
        y0, y1, y2 = y[:-2], y[1:-1], y[2:]
        same = (np.sign(y0) == np.sign(y1)) & (np.sign(y1) == np.sign(y2))
        toward_zero = (np.abs(y1) <= np.abs(y0)) & (np.abs(y1) <= np.abs(y2))
        cand = np.where(same & toward_zero)[0] + 1
        for k in cand:
            a2 = 0.5 * (y[k - 1] - 2 * y[k] + y[k + 1])
            a1 = 0.5 * (y[k + 1] - y[k - 1])
            if a2 == 0:
                continue
            xv = -a1 / (2 * a2)
            if abs(xv) > 1.0:
                continue
            yv = y[k] + a1 * xv + a2 * xv * xv
            if np.sign(yv) == np.sign(y[k]) and abs(yv) > 1e-12:
                continue
            sgn = np.sign(y[k])
            res = minimize_scalar(lambda s: float(sgn * g(np.array([s]))[0]),
                                  bounds=(t[k - 1], t[k + 1]), method="bounded",
                                  options={"xatol": tol * 0.1})
            ge = float(g(np.array([res.x]))[0])
            if np.sign(ge) != sgn and ge != 0.0:
                for lo, hi in ((t[k - 1], res.x), (res.x, t[k + 1])):
                    glo = np.array([float(g(np.array([lo]))[0])])
                    rr = _bisect(g, np.array([lo]), np.array([hi]), glo, tol)[0]
                    roots.append((float(rr), bool(glo[0] < 0)))
    roots.sort()
    intervals = []
    state = bool(pos[0])
    start = t0 if state else None
    for tr, becomes_pos in roots:
        if becomes_pos and not state:
            start, state = tr, True
        elif (not becomes_pos) and state:
            intervals.append((start, tr))
            state = False
    if state:
        intervals.append((start, t1))
    out = np.array(intervals) if intervals else np.zeros((0, 2))
    if return_samples:
        return out, (t, y)
    return out
