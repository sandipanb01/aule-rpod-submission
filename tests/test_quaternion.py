import numpy as np
from scipy.spatial.transform import Rotation

from aule_rpod.attitude.quaternion import attitude_matrix, body_to_reference


def test_stk_convention_canonical_case():
    q = np.array([0, 0, np.sin(np.pi / 4), np.cos(np.pi / 4)])       # 90 deg about Z, scalar last
    assert np.allclose(body_to_reference(q, [1, 0, 0]), [0, 1, 0])  # body X -> ref +Y


def test_attitude_matrix_is_transpose_of_active_rotation():
    rng = np.random.default_rng(0)
    for _ in range(50):
        q = rng.normal(size=4)
        q /= np.linalg.norm(q)
        A = attitude_matrix(q)
        assert np.allclose(A @ A.T, np.eye(3), atol=1e-12)
        assert np.allclose(A.T, Rotation.from_quat(q).as_matrix(), atol=1e-12)


def test_attitude_file_is_reproduced_by_identified_law(sc):
    from aule_rpod.attitude.identification import identify_precession_law
    law = identify_precession_law(sc.profile)
    assert law.max_residual_rad < 1e-9
    assert abs(np.degrees(law.rate_p_rad_s) - 1.2) < 1e-9
    assert abs(law.nutation_deg - 10.0) < 1e-6
    assert np.allclose(np.abs(law.axis_ref), [0, 1, 0], atol=1e-9)    # precession about +Y (V-bar)
    assert np.allclose(np.abs(law.axis_body), [1, 0, 0], atol=1e-9)   # spin about body X


def test_port_axis_stays_within_5_to_15_deg_of_vbar(sc):
    """Exact law: port normal is always 10 +/- 5 deg from +V-bar. SLERP between the
    60-s samples deviates by <= 0.19 deg (documented interpolation error)."""
    from aule_rpod.attitude.providers import LawAttitude
    t = np.linspace(sc.t_start, sc.t_end, 20000)
    law = LawAttitude(sc.profile)
    off = np.degrees(np.arccos(np.clip(law.rotation(t).apply(sc.port_body)[:, 1], -1, 1)))
    tilt = np.degrees(np.arctan2(0.0872, 0.9962))                  # port normal 5.0025 deg off body X
    assert abs(off.min() - (10.0 - tilt)) < 1e-3 and abs(off.max() - (10.0 + tilt)) < 1e-3
    off_s = np.degrees(np.arccos(np.clip(sc.port_axis_lvlh(t)[:, 1], -1, 1)))
    assert off_s.min() > 4.8 and off_s.max() < 15.2
