"""Closed-loop Monte Carlo: with every error source switched off it must fly the
nominal plan; with errors on it must be reproducible from its seed."""
import numpy as np

from aule_rpod.rpod import montecarlo
from aule_rpod.rpod.design import Strategy, build_plan

N_MEAN = 1.0948236928858023e-3
X0 = np.r_[-80.0, 120.0, 20.0, 0.0, 0.0, 0.0]


def _plan():
    return build_plan(N_MEAN, X0, Strategy("vbar_plus", "+V-bar", np.array([0.0, 1.0, 0.0])),
                      235.0, 122.0, 30.0, 25.0, 10.0, 0.06, 0.06, 0.0)


def test_noise_free_run_reproduces_nominal_plan():
    plan = _plan()
    s, o = montecarlo.run(N_MEAN, X0, plan, 50, 1, 0.0, 0.0, 0.0, 0.0, 0.0)
    assert s["hp_position_error_m_mean_p99"][1] < 1e-9
    assert s["leg1_mcc_dv_mm_s_mean_p99"][1] < 1e-6 and s["leg2_mcc_dv_mm_s_mean_p99"][1] < 1e-6
    lo, hi = s["contact_closing_speed_mps_p01_p99"]
    assert abs(lo - 0.06) < 1e-4 and abs(hi - 0.06) < 1e-4
    # the controller replaces the analytic 2nL term; deadband keeps it within a few mm/s of it
    assert abs(s["dv_total_mps_mean_p99"][0] - plan.dv_total) < 5e-3
    assert s["contact_lateral_misalignment_m_mean_p99"][1] < 0.02
    assert s["idss_translational_envelope_fraction"] == 1.0


def test_seeded_runs_are_reproducible():
    plan = _plan()
    args = (N_MEAN, X0, plan, 300, 7, 0.01, 0.01, 0.001, 0.01, 0.5)
    a, _ = montecarlo.run(*args)
    b, _ = montecarlo.run(*args)
    assert a == b
    assert a["collision_radius_entry_probability"] == 0.0
