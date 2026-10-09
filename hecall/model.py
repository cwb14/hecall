"""Origin HMM over equal-information windows.

States = polyploid subgenomes. A window's observation is its vector of SNV divergences to every
reference; the emission density of each state is learned (KDE) from the polyploid's own windows
of that subgenome, so native heterogeneity (incomplete lineage sorting, introgression, divergent
reference haplotypes, hard-to-align regions) is part of each state's expected spread rather than
evidence of exchange. Single-site evidence is used only to place breakpoints (see design 2.3).
"""
import math

import numpy as np
from scipy.ndimage import gaussian_filter

from . import params as PR


# ----------------------------------------------------------------------------- features
def features(cov, mm):
    """(n, K) compared bases / mismatches -> (n, K-1) centred log-divergence contrasts.

    Centring removes the local mutation-rate / constraint level shared by all references; what is
    left is how much closer the window is to one reference than to another."""
    y = np.log((mm + 0.5) / (cov + 1.0))
    K = y.shape[1]
    H = np.zeros((K - 1, K))  # orthonormal Helmert contrasts: rows sum to zero
    for i in range(1, K):
        H[i - 1, :i] = 1.0 / i
        H[i - 1, i] = -1.0
        H[i - 1] /= np.linalg.norm(H[i - 1])
    return y @ H.T


def lda(X, lab, S):
    """Projection onto the min(dim, S-1) directions that best separate the labelled states."""
    d = X.shape[1]
    k = min(d, S - 1)
    if d <= k:
        return np.eye(d)
    mu = X.mean(axis=0)
    Sw = np.zeros((d, d))
    Sb = np.zeros((d, d))
    for s in range(S):
        Xs = X[lab == s]
        if Xs.shape[0] < 2:
            continue
        ms = Xs.mean(axis=0)
        Sw += np.cov(Xs, rowvar=False) * (Xs.shape[0] - 1)
        Sb += Xs.shape[0] * np.outer(ms - mu, ms - mu)
    w, V = np.linalg.eig(np.linalg.pinv(Sw) @ Sb)
    order = np.argsort(-w.real)[:k]
    return V[:, order].real


class KDE:
    """Gridded Gaussian KDE (Scott bandwidth) with a floor of one pseudo-window over the grid."""

    def __init__(self, X, lo, hi):
        d = X.shape[1]
        nb = PR.KDE_GRID[d]
        self.edges = [np.linspace(lo[j], hi[j], nb + 1) for j in range(d)]
        H, _ = np.histogramdd(X, bins=self.edges)
        n = max(X.shape[0], 1)
        sd = X.std(axis=0) if X.shape[0] > 1 else np.ones(d)
        bw = np.maximum(sd, 1e-6) * n ** (-1.0 / (d + 4))
        sig = bw / ((np.asarray(hi) - np.asarray(lo)) / nb)
        G = gaussian_filter(H, sig, mode="constant")
        G = G / G.sum() if G.sum() > 0 else np.full_like(G, 1.0 / G.size)
        self.logG = np.log(G + 1.0 / (n * G.size))

    def logpdf(self, X):
        idx = tuple(np.clip(np.searchsorted(self.edges[j], X[:, j]) - 1, 0, len(self.edges[j]) - 2)
                    for j in range(X.shape[1]))
        return self.logG[idx]


# ----------------------------------------------------------------------------- transitions / HMM
def trans_probs(d, rho, S):
    """Stay / move-to-one-specific-other-state probabilities over d bp (uniform switching)."""
    e = np.exp(-S * rho * d / (S - 1))
    stay = 1.0 / S + (S - 1.0) / S * e
    move = (1.0 - stay) / (S - 1)
    return stay, move


def forward_backward(binned, labels, pi, rho, S):
    """Lockstep scaled forward-backward over all chromosomes (padded to the longest).

    binned: per chromosome dict with E (n, S) log emissions and d (n,) bp from previous unit.
    pi: (S, S) start distribution per chromosome label. Returns (posteriors, logL, expected
    switches, total distance)."""
    C = len(binned)
    N = max((x["E"].shape[0] for x in binned), default=0)
    if N == 0:
        return [np.zeros((0, S)) for _ in binned], 0.0, 0.0, 0.0
    E = np.zeros((C, N, S))
    D = np.zeros((C, N))
    n = np.zeros(C, dtype=np.int64)
    for i, x in enumerate(binned):
        m = x["E"].shape[0]
        n[i] = m
        if m:
            E[i, :m] = x["E"]
            D[i, :m] = x["d"]
    mx = E.max(axis=2, keepdims=True)
    e = np.exp(E - mx)
    stay, move = trans_probs(D, rho, S)
    alpha = np.empty((C, N, S))
    scale = np.empty((C, N))
    a = pi[labels] * e[:, 0]
    s = a.sum(axis=1)
    alpha[:, 0] = a / s[:, None]
    scale[:, 0] = s
    for t in range(1, N):
        prev = alpha[:, t - 1]
        pred = prev * stay[:, t, None] + (1.0 - prev) * move[:, t, None]
        a = pred * e[:, t]
        s = a.sum(axis=1)
        alpha[:, t] = a / s[:, None]
        scale[:, t] = s
    beta = np.ones((C, N, S))
    xi_sw = np.zeros(C)
    for t in range(N - 1, 0, -1):
        eb = e[:, t] * beta[:, t]
        tot = eb.sum(axis=1)
        beta[:, t - 1] = (stay[:, t, None] * eb + move[:, t, None] * (tot[:, None] - eb)) / scale[:, t, None]
        prev = alpha[:, t - 1]
        xi_sw += move[:, t] / scale[:, t] * (tot - (prev * eb).sum(axis=1))
    gamma = alpha * beta
    gamma /= gamma.sum(axis=2, keepdims=True)
    valid = np.arange(N)[None, :] < n[:, None]
    logL = float((np.log(scale) + mx[..., 0])[valid].sum())
    return [gamma[i, :n[i]] for i in range(C)], logL, float(xi_sw.sum()), float(D[valid].sum())


def fit_hmm(binned, labels, S, rho=PR.RHO_INIT, verbose=None):
    """EM for the breakpoint rate and the start distribution (emissions fixed)."""
    pi = np.full((S, S), 1e-3 / max(S - 1, 1))
    np.fill_diagonal(pi, 1 - 1e-3)
    last = None
    for it in range(PR.EM_MAX_ITER):
        post, logL, sw, dist = forward_backward(binned, labels, pi, rho, S)
        if verbose:
            verbose(f"  EM {it}: logL={logL:.1f} rho={rho:.3g} switches={sw:.1f}")
        if dist > 0:
            rho = min(max(sw / dist, PR.RHO_MIN), PR.RHO_MAX)
        newpi = np.ones((S, S))
        for g, lab in zip(post, labels):
            if g.size:
                newpi[lab] += g.sum(axis=0)
        pi = newpi / newpi.sum(axis=1, keepdims=True)
        if last is not None and abs(logL - last) <= PR.EM_TOL * abs(last):
            break
        last = logL
    post, logL, _, _ = forward_backward(binned, labels, pi, rho, S)
    return post, rho, pi, logL


def autocorrelation(series, max_lag=None):
    """Autocorrelation of window evidence at lags 0..max_lag, pooled over chromosomes, truncated
    where Geyer's initial positive sequence ends (first non-positive sum of a lag pair -> 0).

    No fixed cap by default (max lag = longest series): in genomes whose native evidence wanders
    over megabases (stand-in references, introgression), a cap would understate the correlation
    of long blocks and overstate their significance (v0.4.6, Avena agadiriana)."""
    if max_lag is None:
        max_lag = max([x.size - 1 for x in series] + [1])
    num = np.zeros(max_lag + 1)
    den = 0.0
    for x in series:
        if x.size < 3:
            continue
        x = x - x.mean()
        m = 1 << int(np.ceil(np.log2(2 * x.size)))
        f = np.fft.rfft(x, m)
        ac = np.fft.irfft(f * np.conj(f), m)[:min(max_lag, x.size - 1) + 1]
        num[:ac.size] += ac
        den += ac[0]
    rho = num / den if den > 0 else np.zeros(max_lag + 1)
    rho[0] = 1.0
    k = 1
    while k + 1 <= max_lag and rho[k] + rho[k + 1] > 0:
        k += 2
    rho[k:] = 0.0
    return rho


def tau_n(rho, n):
    """Variance inflation of a sum of n consecutive correlated window evidences:
    1 + 2 * sum_{k=1}^{n-1} (1 - k/n) rho_k (composite-likelihood / effective-sample-size factor)."""
    n = int(max(n, 1))
    k = np.arange(1, min(n, rho.size))
    return max(1.0, 1.0 + 2.0 * float(((1.0 - k / n) * rho[k]).sum()))


def run(wins, S, verbose=None):
    """wins: per chromosome dict(label, X (n, d) features, d (n,) distances). Rounds: fit KDEs on
    native windows of each state -> HMM -> reassign -> repeat.

    A state's density is learned only from windows on chromosomes labelled with that state and
    currently assigned to it (start: all of them). Windows assigned to s on other subgenomes'
    chromosomes are the candidate exchanges being tested; letting them shape s's density made it
    absorb ambiguous windows (equidistant from all references) round after round (v0.4.0).
    The HMM proposes segments; their evidence is calibrated for serial correlation afterwards
    (calls.py), using the autocorrelation of window evidence in native windows returned here."""
    labels = np.array([w["label"] for w in wins], dtype=np.int64)
    Xall = np.vstack([w["X"] for w in wins if w["X"].shape[0]])
    lo = np.quantile(Xall, PR.KDE_RANGE_Q, axis=0) - PR.KDE_PAD
    hi = np.quantile(Xall, 1 - PR.KDE_RANGE_Q, axis=0) + PR.KDE_PAD
    assign = [np.full(w["X"].shape[0], w["label"], dtype=np.int64) for w in wins]
    L = np.concatenate(assign)
    rho, prev, kdes = PR.RHO_INIT, None, None
    for r in range(PR.PROFILE_ROUNDS + 1):
        A = np.concatenate(assign)
        kdes = [KDE(Xall[(A == s) & (L == s)], lo, hi) for s in range(S)]
        raw = [np.column_stack([k.logpdf(w["X"]) for k in kdes]) if w["X"].shape[0] else np.zeros((0, S))
               for w in wins]
        series = []
        for E, w, a in zip(raw, wins, assign):
            if E.shape[0]:
                own = E[np.arange(E.shape[0]), a]
                other = np.where(np.eye(S, dtype=bool)[a], -np.inf, E).max(axis=1)
                series.append((other - own)[a == w["label"]])
        rho_ac = autocorrelation(series)
        binned = [{"E": E, "d": w["d"]} for E, w in zip(raw, wins)]
        post, rho, pi, logL = fit_hmm(binned, labels, S, rho, verbose)
        assign = [g.argmax(axis=1) if g.size else np.zeros(0, dtype=np.int64) for g in post]
        cur = np.concatenate(assign)
        changed = None if prev is None else int((cur != prev).sum())
        if verbose:
            verbose(f" round {r}: rho={rho:.3g} logL={logL:.1f} windows changed={changed}")
        if changed == 0:
            break
        prev = cur
    return {"kdes": kdes, "lo": lo, "hi": hi, "binned": binned, "post": post, "rho": rho, "pi": pi,
            "logL": logL, "assign": assign, "rho_ac": rho_ac}


def tempered_posteriors(fit, wins, S, rho, temper):
    """Posteriors with window evidence divided by `temper` and the breakpoint rate fixed at rho.

    Used to propose candidates at coarser scales: with correlated windows the raw evidence of short
    noisy runs inside a long exchange beats the switch cost and shatters it; tempered evidence
    (temper ~ integrated autocorrelation time) lets long blocks stand. Scoring is unchanged."""
    labels = np.array([w["label"] for w in wins], dtype=np.int64)
    binned = [{"E": b["E"] / temper, "d": b["d"]} for b in fit["binned"]]
    post, _, _, _ = forward_backward(binned, labels, fit["pi"], rho, S)
    return post


def integrated_tau(rho_ac):
    """Integrated autocorrelation time of window evidence (truncated sequence from autocorrelation)."""
    return max(1.0, 1.0 + 2.0 * float(rho_ac[1:].sum()))


# ----------------------------------------------------------------------------- site profiles (breakpoints only)
def site_profiles(chroms, S, n_pat):
    """log P(diagnostic-site pattern | state); used only to place cuts between windows."""
    cnt = np.zeros((S, n_pat))
    for c in chroms:
        st = c.get("site_state")
        if st is None:
            st = np.full(c["site_pid"].size, c["label"])
        np.add.at(cnt, (st, c["site_pid"]), 1.0)
    p = cnt + PR.PSEUDOCOUNT
    return np.log(p / p.sum(axis=1, keepdims=True))


# ----------------------------------------------------------------------------- segments
def _runs(states):
    if states.size == 0:
        return []
    br = np.flatnonzero(np.diff(states)) + 1
    starts = np.concatenate([[0], br])
    ends = np.concatenate([br, [states.size]])
    return [(int(a), int(b), int(states[a])) for a, b in zip(starts, ends)]


def refine_breakpoint(pos, lp_left, lp_right):
    """Cut between consecutive sites maximising likelihood; P(cut) ∝ exp(loglik).

    Cut index i: sites[:i] left, sites[i:] right (i in 1..n-1). Returns (i, ml_lo, ml_hi, ci_lo,
    ci_hi): the ML cut lies in [ml_lo, ml_hi] (bp); the credible interval likewise."""
    n = pos.size
    if n < 2:
        return None
    cl = np.cumsum(lp_left)
    cr = np.cumsum(lp_right[::-1])[::-1]
    i = np.arange(1, n)
    ll = cl[i - 1] + cr[i]
    w = np.exp(ll - ll.max())
    w /= w.sum()
    k = int(np.argmax(ll))
    cdf = np.cumsum(w)
    tail = (1 - PR.BREAKPOINT_CI) / 2
    klo = int(np.searchsorted(cdf, tail))
    khi = min(int(np.searchsorted(cdf, 1 - tail)), n - 2)
    return (int(i[k]), int(pos[i[k] - 1]) + 1, int(pos[i[k]]), int(pos[i[klo] - 1]) + 1, int(pos[i[khi]]))


def segments(chrom, win, post, binned, slogp):
    """Runs of max-posterior state over windows; breakpoints refined to diagnostic-site resolution
    inside the two windows flanking each switch (cuts kept in order, never crossing)."""
    if post.size == 0:
        return []
    st = post.argmax(axis=1)
    ws, we = win["start"], win["end"]
    pos, pid = chrom["site_pos"], chrom["site_pid"]
    L = chrom["length"]
    segs = []
    for a, b, s in _runs(st):
        segs.append({"state": s, "win_a": a, "win_b": b, "post_mean": float(post[a:b, s].mean()),
                     "llr": float((binned["E"][a:b, s] - binned["E"][a:b, chrom["label"]]).sum())})
    prev_cut = 0
    for left, right in zip(segs[:-1], segs[1:]):
        lo_bp = max(int(ws[left["win_b"] - 1]), prev_cut)
        hi_bp = int(we[right["win_a"]])
        ra, rz = np.searchsorted(pos, [lo_bp, hi_bp])
        bp = refine_breakpoint(pos[ra:rz], slogp[left["state"]][pid[ra:rz]], slogp[right["state"]][pid[ra:rz]])
        if bp is None:  # fewer than two sites: cut at the window boundary
            x = int(ws[right["win_a"]])
            bp = (0, x, x, int(ws[left["win_b"] - 1]), int(we[right["win_a"]]))
        right["bp_left"] = left["bp_right"] = bp[1:]
        prev_cut = bp[2]
    for i, s in enumerate(segs):
        s["terminal_left"] = i == 0
        s["terminal_right"] = i == len(segs) - 1
        s.setdefault("bp_left", (0, 0, 0, 0))
        s.setdefault("bp_right", (L, L, L, L))
        # chromosome ends without evaluated windows are not painted (they go to nocall.bed)
        s["tile_start"] = int(ws[s["win_a"]]) if i == 0 else (s["bp_left"][0] + s["bp_left"][1]) // 2
        s["tile_end"] = int(we[s["win_b"] - 1]) if i == len(segs) - 1 else (s["bp_right"][0] + s["bp_right"][1]) // 2
        s["site_a"], s["site_b"] = (int(x) for x in np.searchsorted(pos, [s["tile_start"], s["tile_end"]]))
        s["qual"] = qual(s["post_mean"])
    return segs


def qual(post_mean):
    return min(PR.QUAL_CAP, -10.0 * math.log10(max(1.0 - post_mean, 10 ** (-PR.QUAL_CAP / 10))))


# ----------------------------------------------------------------------------- resolution
def resolution(fit, wins, S, rho, rho_ac):
    """Per ordered pair (recipient s, donor t): expected evidence per window of donor sequence (mean
    over t's own windows of log f_t - log f_s), typical window span on s, and the shortest HE whose
    serially-calibrated evidence (n * kl / tau_n) beats two breakpoints at rate rho by QUAL 20."""
    E = np.vstack([b["E"] for b in fit["binned"] if b["E"].shape[0]])
    A = np.concatenate(fit["assign"])
    out = []
    for s in range(S):
        ws = [w for w in wins if w["label"] == s]
        spans = [w["end"] - w["start"] for w in ws if w["end"].size]
        span = float(np.median(np.concatenate(spans))) if spans else float("nan")
        med_len = float(np.median([w["length"] for w in ws])) if ws else 0.0
        nmax = int(med_len / span) if span > 0 else 0
        for t in range(S):
            if t == s:
                continue
            sel = A == t
            kl = float((E[sel, t] - E[sel, s]).mean()) if sel.any() else 0.0
            prior = 2 * math.log(rho * span / (S - 1)) if span > 0 else -float("inf")
            lstar = float("inf")
            if kl > 0:
                n = 1
                while n <= nmax:
                    if n * kl / tau_n(rho_ac, n) + prior >= math.log(99):
                        lstar = n * span
                        break
                    n = n + 1 if n < 50 else int(n * 1.05) + 1
            out.append({"recipient": s, "donor": t, "llr_per_window": kl, "window_bp": span,
                        "tau_at_min": tau_n(rho_ac, int(lstar / span)) if math.isfinite(lstar) else float("nan"),
                        "min_he_bp": lstar, "median_chrom_bp": med_len, "resolvable": math.isfinite(lstar)})
    return out
