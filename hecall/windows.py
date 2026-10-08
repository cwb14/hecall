"""Equal-information windows along each polyploid chromosome.

A window grows over consecutive 1-kb blocks until every reference has at least `cmin` compared
bases, where cmin gives ~WINDOW_MISMATCHES expected mismatches to the closest reference at this
genome's typical divergence. Windows are therefore short where alignments are dense and divergence
high (fine resolution) and long where they are sparse; a window that would have to span more than
WINDOW_MAX_FACTOR x the median span is not formed (region not evaluable).
"""
import numpy as np

from . import params as PR


def typical_closest_divergence(data):
    """Median over 10-kb windows of the divergence to the closest reference (all subgenomes)."""
    vals = []
    for d in data.values():
        cov, mm = d["cov"], d["mm"]
        nb = cov.shape[1]
        g = max(1, 10_000 // PR.BIN_BP)  # 10-kb blocks
        n = nb // g * g
        if n == 0:
            continue
        c = cov[:, :n].reshape(cov.shape[0], -1, g).sum(axis=2)
        m = mm[:, :n].reshape(mm.shape[0], -1, g).sum(axis=2)
        ok = c >= PR.WINDOW_MIN_COMPARED
        div = np.where(ok, m / np.maximum(c, 1), np.inf).min(axis=0)
        vals.append(div[np.isfinite(div)])
    v = np.concatenate(vals) if vals else np.zeros(0)
    return float(np.median(v)) if v.size else float("nan")


def _ends(cov, cmin):
    """For each start block i: first end block e (exclusive) with every reference >= cmin."""
    K, nb = cov.shape
    e = np.zeros(nb, dtype=np.int64)
    for k in range(K):
        cum = np.concatenate([[0], np.cumsum(cov[k], dtype=np.int64)])
        ek = np.searchsorted(cum, cum[:-1] + cmin, side="left")  # cum[ek] >= cum[i] + cmin
        e = np.maximum(e, ek)
    return e  # value nb+1 means never reached


def build(cov, mm, cmin, maxspan):
    """Greedy left-to-right tiling. Returns (starts, ends) in blocks; blocks not in any window are
    not evaluable."""
    nb = cov.shape[1]
    e = _ends(cov, cmin)
    starts, ends = [], []
    i = 0
    while i < nb:
        if e[i] <= nb and e[i] - i <= maxspan:
            starts.append(i)
            ends.append(int(e[i]))
            i = int(e[i])
        else:
            i += 1
    return np.asarray(starts, dtype=np.int64), np.asarray(ends, dtype=np.int64)


def window_sums(arr, starts, ends):
    """Sum of per-block counts (K, nb) over each window -> (n, K)."""
    cum = np.concatenate([np.zeros((arr.shape[0], 1), dtype=np.int64), np.cumsum(arr, axis=1, dtype=np.int64)], axis=1)
    return (cum[:, ends] - cum[:, starts]).T
