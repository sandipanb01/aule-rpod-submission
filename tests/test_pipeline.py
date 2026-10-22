"""End-to-end regression of the headline results (fast path, ~20 s)."""
import numpy as np

from aule_rpod.events import intervals as I
from aule_rpod.windows import compute_constraints, docking_windows, pass_diagnostics


def test_part1_headline(sc, cfg):
    cs = compute_constraints(sc)
    W = docking_windows(cs, cfg.constraints.docking_duration_s)
    assert len(W) == 0                                     # no strict window
    qual = {k: len(I.filter_min_duration(v, 540.0)) for k, v in cs.passes.items()}
    assert qual == {"Svalbard": 80, "Bengaluru": 15, "Wellington": 19, "Los Angeles": 16}
    d = pass_diagnostics(sc, cs, 540.0, only_qualifying=True)
    best = max(d, key=lambda r: r.best_feasible_s)
    assert best.station == "Wellington" and abs(best.best_feasible_s - 462.1) < 1.0
    assert best.limited_start_by.startswith("sun angle < 30")


def test_part2_plan_fits_window(cfg, sc):
    from aule_rpod.rpod.design import Strategy, build_plan
    n = 1.0948237e-3
    p = build_plan(n, np.r_[-80, 120, 20, 0, 0, 0.0], Strategy("vbar", "+V-bar", np.array([0, 1.0, 0])),
                   235, 122, 30, 25, 10, 0.06, 0.06, 0.0)
    assert p.t_contact < 583.9 - 25
    assert 1.1 < p.dv_total < 1.3
    assert np.allclose(p.r[-1], 0.0)
