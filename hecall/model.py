"""Origin HMM: states = polyploid subgenomes, emissions = learned site-pattern profiles.

Data for one polyploid (all chromosomes):
    chroms: list of dicts with keys
        name, length, label (state index),
        site_pos (int64, sorted), site_pid (int pattern id)
Profiles: logp (S, P) log-probabilities of pattern id per state.
"""
import math
import numpy as np

from . import params as PR


# ----------------------------------------------------------------------------- profiles
def pool_rare(chroms, n_pat):
    """Map pattern ids seen < MIN_PATTERN_COUNT times overall to one shared id. Returns map."""
    tot = np.zeros(n_pat, dtype=np.int64)
    for c in chroms:
        tot += np.bincount(c["site_pid"], minlength=n_pat)
    keep = tot >= PR.MIN_PATTERN_COUNT
    remap = np.empty(n_pat, dtype=np.int64)
    remap[keep] = np.arange(keep.sum())
    remap[~keep] = keep.sum()  # 'rare'
    return remap, int(keep.sum()) + int((~keep).any())


def profile_counts(chroms, S, n_pat, site_state):
    """Counts (S, n_pat) of pattern ids per assigned state. site_state: list of arrays."""
    cnt = np.zeros((S, n_pat), dtype=np.float64)
    for c, st in zip(chroms, site_state):
        np.add.at(cnt, (st, c["site_pid"]), 1.0)
    return cnt


def log_profiles(cnt):
    p = cnt + PR.PSEUDOCOUNT
    p /= p.sum(axis=1, keepdims=True)
    return np.log(p)


# ----------------------------------------------------------------------------- transitions
def trans_probs(d, rho, S):
    """Stay / move-to-one-specific-other-state probabilities over d bp (uniform switching)."""
    e = np.exp(-S * rho * d / (S - 1))
    stay = 1.0 / S + (S - 1.0) / S * e
    move = (1.0 - stay) / (S - 1)
    return stay, move


# ----------------------------------------------------------------------------- binning
def bin_chroms(chroms, logp):
    """Per chromosome: nonempty bin index, distances, and (nb, S) emission log-likelihoods."""
    out = []
    for c in chroms:
        b = c["site_pos"] // PR.BIN_BP
        ub, inv = np.unique(b, return_inverse=True)
        S = logp.shape[0]
        E = np.zeros((ub.size, S))
        for k in range(S):
            E[:, k] = np.bincount(inv, weights=logp[k][c["site_pid"]], minlength=ub.size)
        d = np.diff(ub, prepend=ub[:1]).astype(np.float64) * PR.BIN_BP if ub.size else np.zeros(0)
        out.append({"bins": ub, "inv": inv, "E": E, "d": d})
    return out


# ----------------------------------------------------------------------------- forward-backward
def forward_backward(binned, labels, pi, rho, S):
    """Lockstep scaled forward-backward over all chromosomes (padded to the longest).

    pi: (S, S) start distribution per chromosome label. Returns (posteriors list, logL,
    expected switches, total distance).
    """
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
    e = np.exp(E - mx)                              # (C, N, S), padding rows = 1
    stay, move = trans_probs(D, rho, S)              # padding d=0 -> stay=1, move=0
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
        eb = e[:, t] * beta[:, t]                     # (C, S)
        tot = eb.sum(axis=1)
        beta[:, t - 1] = (stay[:, t, None] * eb + move[:, t, None] * (tot[:, None] - eb)) / scale[:, t, None]
        # expected switches between t-1 and t
        prev = alpha[:, t - 1]
        xi_sw += move[:, t] / scale[:, t] * (tot - (prev * eb).sum(axis=1))
    gamma = alpha * beta
    gamma /= gamma.sum(axis=2, keepdims=True)
    valid = np.arange(N)[None, :] < n[:, None]
    logL = float((np.log(scale) + mx[..., 0])[valid].sum())
    post = [gamma[i, :n[i]] for i in range(C)]
    return post, logL, float(xi_sw.sum()), float(D[valid].sum())


def fit(chroms, logp, S, labels, rho=PR.RHO_INIT, verbose=None):
    """EM for rho and the start distribution with fixed emission profiles."""
    binned = bin_chroms(chroms, logp)
    pi = np.full((S, S), 1e-3 / max(S - 1, 1))
    np.fill_diagonal(pi, 1 - 1e-3)
    last = None
    for it in range(PR.EM_MAX_ITER):
        post, logL, sw, dist = forward_backward(binned, labels, pi, rho, S)
        if verbose:
            verbose(f"  EM {it}: logL={logL:.1f} rho={rho:.3g} switches={sw:.1f}")
        if dist > 0:
            rho = min(max(sw / dist, PR.RHO_MIN), PR.RHO_MAX)
        # start distribution per label = genome-wide posterior state fractions on that label
        newpi = np.full((S, S), 1.0)  # pseudocount 1 per cell
        for g, lab in zip(post, labels):
            if g.size:
                newpi[lab] += g.sum(axis=0)
        pi = newpi / newpi.sum(axis=1, keepdims=True)
        if last is not None and abs(logL - last) <= PR.EM_TOL * abs(last):
            break
        last = logL
    post, logL, sw, dist = forward_backward(binned, labels, pi, rho, S)
    return binned, post, rho, pi, logL


# ----------------------------------------------------------------------------- full fit with profile rounds
def run(chroms, S, n_pat, verbose=None):
    """Profile rounds: estimate profiles -> fit HMM -> reassign sites to called state -> repeat."""
    remap, P = pool_rare(chroms, n_pat)
    for c in chroms:
        c["site_pid"] = remap[c["site_pid"]]
    labels = np.array([c["label"] for c in chroms], dtype=np.int64)
    site_state = [np.full(c["site_pid"].size, c["label"], dtype=np.int64) for c in chroms]
    rho = PR.RHO_INIT
    prev_calls = None
    for r in range(PR.PROFILE_ROUNDS + 1):
        cnt = profile_counts(chroms, S, P, site_state)
        logp = log_profiles(cnt)
        binned, post, rho, pi, logL = fit(chroms, logp, S, labels, rho, verbose)
        calls = [g.argmax(axis=1) if g.size else np.zeros(0, dtype=np.int64) for g in post]
        site_state = [cl[x["inv"]] if cl.size else np.zeros(0, dtype=np.int64)
                      for cl, x in zip(calls, binned)]
        changed = None if prev_calls is None else sum(int((a != b).sum()) for a, b in zip(calls, prev_calls))
        if verbose:
            verbose(f" profile round {r}: rho={rho:.3g} logL={logL:.1f} bins changed={changed}")
        if changed == 0:
            break
        prev_calls = calls
    return {"logp": logp, "counts": cnt, "remap": remap, "n_pat": P, "binned": binned,
            "post": post, "rho": rho, "pi": pi, "logL": logL}


# ----------------------------------------------------------------------------- segments
def _runs(states):
    """Runs of equal values: list of (start_idx, end_idx_exclusive, value)."""
    if states.size == 0:
        return []
    br = np.flatnonzero(np.diff(states)) + 1
    starts = np.concatenate([[0], br])
    ends = np.concatenate([br, [states.size]])
    return [(int(a), int(b), int(states[a])) for a, b in zip(starts, ends)]


def refine_breakpoint(pos, lp_left, lp_right):
    """Cut between consecutive sites maximising likelihood; P(cut) ∝ exp(loglik).

    pos: positions of sites in the boundary region (sorted); lp_left/right: per-site log-probs
    under the left / right state. Cut index i means sites[:i] left, sites[i:] right (i in 1..n-1).
    Returns (ml_lo, ml_hi, ci_lo, ci_hi): ML cut lies in (pos[i-1], pos[i]]; CI likewise.
    """
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
    khi = int(np.searchsorted(cdf, 1 - tail))
    khi = min(khi, n - 2)
    return (int(pos[i[k] - 1]) + 1, int(pos[i[k]]), int(pos[i[klo] - 1]) + 1, int(pos[i[khi]]))


def segments(chrom, binned, post, logp, edge_bins=2):
    """Segments of max-posterior state along one chromosome, with site-level breakpoints."""
    g = post
    if g.size == 0:
        return []
    st = g.argmax(axis=1)
    runs = _runs(st)
    pos, pid, inv = chrom["site_pos"], chrom["site_pid"], binned["inv"]
    bins = binned["bins"]
    segs = []
    # site index range per nonempty bin
    first_site = np.searchsorted(inv, np.arange(bins.size), side="left")
    last_site = np.searchsorted(inv, np.arange(bins.size), side="right")
    for r, (a, b, s) in enumerate(runs):
        seg = {"state": s, "bin_a": a, "bin_b": b,
               "site_a": int(first_site[a]), "site_b": int(last_site[b - 1]),
               "post_mean": float(g[a:b, s].mean())}
        segs.append(seg)
    # refine each internal boundary
    for left, right in zip(segs[:-1], segs[1:]):
        la = max(left["bin_a"], left["bin_b"] - edge_bins)
        rb = min(right["bin_b"], right["bin_a"] + edge_bins)
        sa, sb = int(first_site[la]), int(last_site[rb - 1])
        p = pos[sa:sb]
        bp = refine_breakpoint(p, logp[left["state"]][pid[sa:sb]], logp[right["state"]][pid[sa:sb]])
        if bp is None:
            bp = (int(pos[left["site_b"] - 1]) + 1, int(pos[right["site_a"]]),) * 2
        right["bp_left"] = bp
        left["bp_right"] = bp
        # move the site boundary to the ML cut
        cut = int(np.searchsorted(pos, bp[1]))
        left["site_b"] = cut
        right["site_a"] = cut
    L = chrom["length"]
    for s in segs:
        s["terminal_left"] = "bp_left" not in s
        s["terminal_right"] = "bp_right" not in s
        s.setdefault("bp_left", (0, 0, 0, 0))
        s.setdefault("bp_right", (L, L, L, L))
        if s["site_b"] > s["site_a"]:
            s["start"] = int(pos[s["site_a"]])
            s["end"] = int(pos[s["site_b"] - 1]) + 1
        else:
            s["start"] = s["end"] = None
    return segs


def qual(post_mean):
    return min(PR.QUAL_CAP, -10.0 * math.log10(max(1.0 - post_mean, 10 ** (-PR.QUAL_CAP / 10))))


# ----------------------------------------------------------------------------- resolution
def resolution(logp, chroms, S, rho):
    """For each ordered pair (s, t): KL per site, site density on s, min detectable HE length."""
    p = np.exp(logp)
    out = []
    for s in range(S):
        cs = [c for c in chroms if c["label"] == s]
        nsite = sum(c["site_pos"].size for c in cs)
        blen = sum(c["length"] for c in cs)
        dens = nsite / blen if blen else 0.0
        med_len = float(np.median([c["length"] for c in cs])) if cs else 0.0
        for t in range(S):
            if t == s:
                continue
            kl = float((p[t] * (logp[t] - logp[s])).sum())
            need = math.log(99) + 2 * math.log(1.0 / (rho * PR.BIN_BP))
            per_bp = kl * dens
            lstar = need / per_bp if per_bp > 0 else float("inf")
            out.append({"recipient": s, "donor": t, "kl_per_site": kl, "sites_per_kb": dens * 1000,
                        "llr_per_kb": per_bp * 1000, "llr_needed": need, "min_he_bp": lstar,
                        "median_chrom_bp": med_len, "resolvable": lstar <= med_len})
    return out
