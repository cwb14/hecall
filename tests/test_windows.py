import numpy as np

from hecall import windows


def test_windows_collect_information_from_every_reference():
    cov = np.array([[1000] * 10, [250] * 10])   # ref 2 aligns sparsely
    mm = np.zeros_like(cov)
    s, e = windows.build(cov, mm, cmin=1000, maxspan=100)
    assert list(s) == [0, 4] and list(e) == [4, 8]     # 4 blocks needed for ref 2; tail too short
    sums = windows.window_sums(cov, s, e)
    assert (sums >= 1000).all()


def test_window_span_cap_leaves_sparse_region_unevaluated():
    cov = np.array([[1000, 1000, 10, 10, 10, 10, 10, 1000, 1000]])
    s, e = windows.build(cov, np.zeros_like(cov), cmin=1000, maxspan=2)
    # blocks 2-5 cannot reach cmin within the cap; block 6 can (with block 7)
    assert list(zip(s, e)) == [(0, 1), (1, 2), (6, 8), (8, 9)]


def test_typical_divergence_when_references_are_the_parents():
    """Most 10-kb blocks without a mismatch (references = parents): blocks grow until the median
    is informative instead of returning 0 (v0.4.7)."""
    rng = np.random.default_rng(5)
    nb = 2000
    cov = np.full((2, nb), 1000)
    mm = np.zeros((2, nb), dtype=np.int64)
    mm[1] = rng.poisson(50, nb)                  # other reference: 5%
    mm[0, rng.choice(nb, 40, replace=False)] = 1  # closest reference: ~2e-5
    d = windows.typical_closest_divergence({"c": {"cov": cov, "mm": mm}})
    assert 0 < d < 1e-4
