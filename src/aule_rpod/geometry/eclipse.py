"""Earth shadow: conical umbra/penumbra model on the WGS-84 ellipsoid.

From the spacecraft, the Earth disk has apparent radius  a_E = asin(R_E/|r|),
the Sun disk  a_S = asin(R_S/|s - r|), and their centres are separated by
d = angle(-r, s - r)  (Montenbruck & Gill 2000, Sec. 3.4.2). Then

    d >= a_E + a_S            full sunlight
    |a_E - a_S| < d < a_E+a_S penumbra (partial occultation)
    d <= a_E - a_S            umbra (total)

Oblateness: the ellipsoid x^2/a^2 + y^2/a^2 + z^2/b^2 = 1 is mapped onto the
sphere of radius a by the affine scaling z -> z a/b applied to BOTH the
spacecraft and the Sun position. Straight lines (Sun rays, limb tangents) are
preserved by affine maps, so umbra/penumbra boundaries are reproduced to within
the (0.34 %) distortion of the 0.27-deg solar cone, i.e. a few hundredths of a
second; verified against Orekit's ellipsoidal OccultationEngine in the V&V.
Set flattening = 0 for the classical spherical model.

Event functions (zero at the boundaries, positive inside the named state):
    g_full_sun = d - (a_E + a_S)          [deg]
    g_umbra    = (a_E - a_S) - d          [deg]
"""
from __future__ import annotations

import numpy as np

from ..constants import R_SUN_KM


def _scale(v, f):
    if f == 0.0:
        return v
    v = np.array(v, dtype=float, copy=True)
    v[..., 2] *= 1.0 / (1.0 - f)
    return v


def shadow_geometry(r, s, re_km, f=0.0):
    """r: spacecraft (N,3) km, s: Sun (N,3) km, geocentric, Earth-equator-aligned frame."""
    r, s = _scale(r, f), _scale(s, f)
    rs = s - r
    rn = np.linalg.norm(r, axis=-1)
    rsn = np.linalg.norm(rs, axis=-1)
    a_E = np.arcsin(np.clip(re_km / rn, -1, 1))
    a_S = np.arcsin(np.clip(R_SUN_KM / rsn, -1, 1))
    d = np.arccos(np.clip(np.sum(-r * rs, -1) / (rn * rsn), -1, 1))
    return a_E, a_S, d


def g_full_sun(r, s, re_km, f=0.0):
    a_E, a_S, d = shadow_geometry(r, s, re_km, f)
    return np.degrees(d - (a_E + a_S))


def g_umbra(r, s, re_km, f=0.0):
    a_E, a_S, d = shadow_geometry(r, s, re_km, f)
    return np.degrees((a_E - a_S) - d)


def sun_fraction(r, s, re_km, f=0.0):
    """Visible fraction of the solar disk (0 = umbra, 1 = full Sun); circular-disk overlap."""
    a_E, a_S, d = shadow_geometry(r, s, re_km, f)
    frac = np.ones_like(d)
    umb = d <= a_E - a_S
    frac[umb] = 0.0
    pen = (~umb) & (d < a_E + a_S)
    if np.any(pen):
        a, b, c = a_S[pen], a_E[pen], d[pen]
        x = (c * c + a * a - b * b) / (2 * c)
        y = np.sqrt(np.clip(a * a - x * x, 0, None))
        area = a * a * np.arccos(np.clip(x / a, -1, 1)) + b * b * np.arccos(np.clip((c - x) / b, -1, 1)) - c * y
        frac[pen] = 1.0 - area / (np.pi * a * a)
    return frac
