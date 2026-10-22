import numpy as np

from aule_rpod.events import intervals as I
from aule_rpod.events.finder import find_positive_intervals


def test_sine_roots():
    g = lambda t: np.sin(2 * np.pi * t / 100.0) - 0.5          # noqa: E731
    iv = find_positive_intervals(g, 0.0, 1000.0, 7.0, tol=1e-7)
    exact = np.array([[100 * (1 / 12 + k), 100 * (5 / 12 + k)] for k in range(10)])
    assert iv.shape == exact.shape and np.abs(iv - exact).max() < 1e-6


def test_hidden_double_crossing_is_found():
    h = lambda t: 0.001 - (t - 50.3) ** 2                        # noqa: E731
    iv = find_positive_intervals(h, 0.0, 100.0, 2.0, tol=1e-8)
    assert len(iv) == 1 and abs((iv[0, 1] - iv[0, 0]) - 2 * np.sqrt(0.001)) < 1e-6


def test_interval_algebra():
    a, b = [[0, 10], [20, 30]], [[5, 25]]
    assert np.allclose(I.intersect(a, b), [[5, 10], [20, 25]])
    assert np.allclose(I.union(a, b), [[0, 30]])
    assert np.allclose(I.complement(a, 0, 40), [[10, 20], [30, 40]])
    assert np.allclose(I.filter_min_duration(a, 10.0), a)
    assert I.containing(a, 25) == 1 and I.containing(a, 15) is None
