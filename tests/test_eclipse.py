import numpy as np

from aule_rpod.geometry.eclipse import g_full_sun, g_umbra, sun_fraction

S = np.array([[1.496e8, 0.0, 0.0]])


def test_states():
    assert g_umbra(np.array([[-7000.0, 0, 0]]), S, 6378.137)[0] > 0
    assert g_full_sun(np.array([[7000.0, 0, 0]]), S, 6378.137)[0] > 0
    assert sun_fraction(np.array([[-7000.0, 0, 0]]), S, 6378.137)[0] == 0.0


def test_penumbra_fraction_monotonic():
    ys = np.linspace(6200, 6600, 200)
    r = np.stack([np.full_like(ys, -7000.0), ys, np.zeros_like(ys)], 1)
    f = sun_fraction(r, np.repeat(S, len(ys), 0), 6378.137)
    assert np.all(np.diff(f) >= -1e-12) and f[0] == 0.0 and f[-1] == 1.0


def test_oblateness_shrinks_shadow_over_pole():
    # a point above the pole's shadow edge is lit with the ellipsoid but not with the sphere
    r = np.array([[-7000.0, 0.0, 6370.0]])
    assert g_full_sun(r, S, 6378.137, f=1 / 298.257223563)[0] > g_full_sun(r, S, 6378.137)[0]
