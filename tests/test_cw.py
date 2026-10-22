import numpy as np
from scipy.integrate import solve_ivp

from aule_rpod.rpod import cw

N = 1.0948e-3


def test_stm_identity_and_semigroup():
    assert np.allclose(cw.stm(N, 0.0), np.eye(6))
    assert np.allclose(cw.stm(N, 700.0), cw.stm(N, 300.0) @ cw.stm(N, 400.0), atol=1e-12)


def test_stm_matches_ode():
    X0 = np.array([-80, 120, 20, 0.1, -0.2, 0.05])
    sol = solve_ivp(lambda t, X: np.r_[X[3:], cw.natural_accel(N, X[:3], X[3:])], (0, 5000), X0, rtol=1e-12, atol=1e-12)
    assert np.abs(cw.propagate(N, X0, 5000.0) - sol.y[:, -1]).max() < 1e-6


def test_two_impulse_hits_target():
    dv1, dv2, v0p, vfm = cw.two_impulse(N, [-80, 120, 20], [0, 0, 0], [0, 25, 0], 235.0)
    X = cw.propagate(N, np.r_[-80, 120, 20, v0p], 235.0)
    assert np.allclose(X[:3], [0, 25, 0], atol=1e-9)
    assert np.allclose(dv2, -vfm)


def test_forced_vbar_translation_dv_is_2nL():
    f = cw.forced_line(N, [0, 10, 0], [0, 0, 0], 0.05)
    assert abs(f["dv_continuous"] - 2 * N * 10) < 1e-12
    assert np.allclose(f["u"][:, 0], 2 * N * 0.05)                 # radial accel 2 n v


def test_vbar_point_is_equilibrium():
    X = cw.propagate(N, [0, 10, 0, 0, 0, 0], np.linspace(0, 6000, 50))
    assert np.allclose(X[:, :3], [0, 10, 0], atol=1e-12)
