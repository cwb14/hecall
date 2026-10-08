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
