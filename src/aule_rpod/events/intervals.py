"""Interval algebra on sorted, disjoint closed intervals stored as (N, 2) arrays."""
from __future__ import annotations

import numpy as np

EMPTY = np.zeros((0, 2))


def as_intervals(x) -> np.ndarray:
    a = np.asarray(x, dtype=float).reshape(-1, 2)
    return a


def normalize(a, merge_gap: float = 0.0) -> np.ndarray:
    a = as_intervals(a)
    if len(a) == 0:
        return EMPTY.copy()
    a = a[np.argsort(a[:, 0])]
    out = [a[0].copy()]
    for s, e in a[1:]:
        if s <= out[-1][1] + merge_gap:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append(np.array([s, e]))
    return np.array(out)


def intersect(a, b) -> np.ndarray:
    a, b = as_intervals(a), as_intervals(b)
    i = j = 0
    out = []
    while i < len(a) and j < len(b):
        s = max(a[i, 0], b[j, 0])
        e = min(a[i, 1], b[j, 1])
        if s < e:
            out.append((s, e))
        if a[i, 1] < b[j, 1]:
            i += 1
        else:
            j += 1
    return np.array(out) if out else EMPTY.copy()


def intersect_all(*sets) -> np.ndarray:
    res = sets[0]
    for s in sets[1:]:
        res = intersect(res, s)
    return res


def union(*sets) -> np.ndarray:
    return normalize(np.vstack([as_intervals(s) for s in sets])) if sets else EMPTY.copy()


def complement(a, lo: float, hi: float) -> np.ndarray:
    a = intersect(normalize(a), [[lo, hi]])
    out = []
    cur = lo
    for s, e in a:
        if s > cur:
            out.append((cur, s))
        cur = max(cur, e)
    if cur < hi:
        out.append((cur, hi))
    return np.array(out) if out else EMPTY.copy()


def durations(a) -> np.ndarray:
    a = as_intervals(a)
    return a[:, 1] - a[:, 0]


def filter_min_duration(a, dmin: float) -> np.ndarray:
    a = as_intervals(a)
    return a[durations(a) >= dmin - 1e-9]


def total(a) -> float:
    return float(durations(a).sum())


def containing(a, t: float):
    """Index of the interval containing t, or None."""
    a = as_intervals(a)
    k = np.searchsorted(a[:, 0], t, side="right") - 1
    if k >= 0 and a[k, 0] <= t <= a[k, 1]:
        return int(k)
    return None
