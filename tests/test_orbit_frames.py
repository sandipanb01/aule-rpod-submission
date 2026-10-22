import erfa
import numpy as np

from aule_rpod.frames import geodetic_to_itrs, itrs_to_geodetic
from aule_rpod.orbit.elements import cart2kep, kep2cart
from aule_rpod.orbit.sso import sso_inclination_rad

MU, RE, J2, F = 398600.4418, 6378.137, 1.08262668e-3, 1 / 298.257223563


def test_sso_inclination_550km():
    i = np.degrees(sso_inclination_rad(RE + 550, 0.0, MU, RE, J2, 365.2421897))
    assert abs(i - 97.593) < 0.001


def test_kepler_roundtrip():
    r, v = kep2cart(7000.0, 0.01, 1.2, 0.5, 0.3, 2.0, MU)
    a, e, i, raan, argp, nu, u = cart2kep(r, v, MU)
    assert np.allclose([a[0], e[0], i[0], raan[0], argp[0], nu[0]], [7000.0, 0.01, 1.2, 0.5, 0.3, 2.0], atol=1e-9)


def test_geodetic_roundtrip():
    r, up, *_ = geodetic_to_itrs(-41.152, 175.505, 0.0, RE, F)
    lat, lon, alt = itrs_to_geodetic(r, RE, F)
    assert abs(lat[0] + 41.152) < 1e-9 and abs(lon[0] - 175.505) < 1e-9 and abs(alt[0]) < 1e-9


def test_fast_c2t_matches_erfa(sc):
    t = np.linspace(sc.t_start, sc.t_end, 50)
    M = sc.eo.gcrs_to_itrs(t)
    tt1, tt2 = sc.ts.jd_tt(t)
    u1, u2 = sc.ts.jd_ut1(t)
    xp, yp = sc.ts.polar_motion_rad(t)
    assert np.abs(M - erfa.c2t06a(tt1, tt2, u1, u2, xp, yp)).max() < 1e-10


def test_ltan_is_flown(sc):
    """Node at local mean solar time 10:30 at epoch (UT1 + lon/15)."""
    r, _ = sc.state_gcrs(np.array([0.0]))
    ri = sc.eo.gcrs_to_itrs(np.array([0.0]))[0] @ r[0]
    lmt = sc.ts.ut1_hours_of_day(0.0) + np.degrees(np.arctan2(ri[1], ri[0])) / 15.0
    assert abs(lmt - 10.5) < 1e-6
