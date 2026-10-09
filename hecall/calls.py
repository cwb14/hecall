"""Segments -> HE calls with quality, independent evidence and artifact flags."""
import os

import numpy as np

import math

from . import extract, model
from .model import tau_n
from . import params as PR


# ----------------------------------------------------------------------------- helpers
def project(aln, a, b):
    """Project query interval [a,b) through alignments, robustly.

    On the target chromosome holding most aligned bases, only the dominant collinear chain is used
    (main strand, diagonal offset within max(CHAIN_TOL_MIN_BP, half the interval) of the weighted
    mode), so repeat alignments elsewhere on the target do not stretch the interval. Returns
    dict(target, lo, hi, t_med, q_med, bp) with lo/hi = 5th/95th and t_med = 50th aligned-base-
    weighted percentile of chain target positions, bp = aligned bases on the chain; or None."""
    if aln is None:
        return None
    sel = np.flatnonzero((aln["qs"] < b) & (aln["qe"] > a))
    if sel.size == 0:
        return None
    x0 = np.maximum(a, aln["qs"][sel])
    x1 = np.minimum(b, aln["qe"][sel])
    fwd = aln["fwd"][sel]
    t0 = np.where(fwd, aln["ts"][sel] + (x0 - aln["qs"][sel]), aln["te"][sel] - (x1 - aln["qs"][sel]))
    t1 = np.where(fwd, aln["ts"][sel] + (x1 - aln["qs"][sel]), aln["te"][sel] - (x0 - aln["qs"][sel]))
    w = (x1 - x0).astype(float)
    tg = aln["t"][sel]
    tot = {}
    for t, ww in zip(tg, w):
        tot[t] = tot.get(t, 0.0) + ww
    best = max(tot, key=tot.get)
    m = tg == best
    # keep the dominant collinear chain: pieces on the main strand whose diagonal offset (target -
    # query, or target + query on the reverse strand) lies near the aligned-base-weighted mode
    tmid, qmid = ((t0 + t1) / 2.0)[m], ((x0 + x1) / 2.0)[m]
    fw_m, w_m = fwd[m], w[m]
    off = np.where(fw_m, tmid - qmid, tmid + qmid)
    strand = fw_m if w_m[fw_m].sum() >= w_m[~fw_m].sum() else ~fw_m
    tol = max(PR.CHAIN_TOL_MIN_BP, 0.5 * (b - a))
    binw = max(PR.CHAIN_TOL_MIN_BP / 10, 0.05 * (b - a))
    o_s = off[strand]
    edges = np.arange(o_s.min(), o_s.max() + 2 * binw, binw)
    hist, _ = np.histogram(o_s, bins=edges, weights=w_m[strand])
    mode = edges[int(np.argmax(hist))] + binw / 2
    chain = strand & (np.abs(off - mode) <= tol)
    m = np.flatnonzero(m)[chain]
    m = np.isin(np.arange(tg.size), m)
    lo_, hi_, wm = np.minimum(t0, t1)[m].astype(float), np.maximum(t0, t1)[m].astype(float), w[m]
    qm = ((x0 + x1) / 2.0)[m]

    def pct(q):
        """Target position below which a fraction q of aligned bases lies (uniform within pieces)."""
        lo, hi = lo_.min(), hi_.max()
        goal = q * wm.sum()
        for _ in range(60):
            x = (lo + hi) / 2
            f = (wm * np.clip((x - lo_) / np.maximum(hi_ - lo_, 1.0), 0, 1)).sum()
            lo, hi = (x, hi) if f < goal else (lo, x)
        return (lo + hi) / 2

    return {"target": best, "lo": int(pct(0.05)), "hi": int(pct(0.95)), "t_med": pct(0.5),
            "q_med": float(np.average(qm, weights=wm)), "bp": int(wm.sum()),
            "target_bp": int(tot[best]), "total_bp": int(w.sum())}


def syntenic(aln, a, b, block_a, block_b):
    """Where does the reference sequence matching [a,b) lie, relative to the sequence matching the
    native flanks of the non-native block [block_a, block_b) that contains it?

    Returns (class, core target interval). Classes:
      homoeologous  collinear with at least one native flank (a rearrangement between assemblies can
                    sit on one side)
      translocation not collinear with any flank, but the matching sequence is one collinear chain
                    on another chromosome or position, spanning at least half the segment and
                    >= SYNTENY_FLANK_MIN_BP: an intergenomic exchange between non-homoeologous
                    positions (a transposed repeat or scattered hits form no such chain)
      nonsyntenic   neither (transposed repeat, scattered hits, misassembly)
      None          no flank aligns (cannot tell)
    Flanks: max(block length, SYNTENY_FLANK_MIN_BP) (capped at SYNTENY_FLANK_MAX_BP) on each side."""
    seg = project(aln, a, b)
    if seg is None:
        return None, None
    F = min(max(block_b - block_a, PR.SYNTENY_FLANK_MIN_BP), PR.SYNTENY_FLANK_MAX_BP)
    verdicts = []
    for fa, fb in ((max(0, block_a - F), block_a), (block_b, block_b + F)):
        if fb <= fa:
            continue
        fl = project(aln, fa, fb)
        if fl is None:
            continue
        dq = abs(seg["q_med"] - fl["q_med"])
        verdicts.append(fl["target"] == seg["target"] and
                        abs(seg["t_med"] - fl["t_med"]) <= PR.SYNTENY_SLACK_BP + PR.SYNTENY_FOLD * dq)
    target = f"{seg['target']}:{seg['lo']}-{seg['hi']}"
    if not verdicts:
        return None, target
    if any(verdicts):
        return "homoeologous", target
    if seg["hi"] - seg["lo"] >= max(PR.SYNTENY_FLANK_MIN_BP, 0.5 * (b - a)):
        return "translocation", target
    return "nonsyntenic", target


def code_stats(work, source, chrom, a, b):
    """(covered SNV-comparable bases, mismatches) of one source over [a,b), from the base codes."""
    f = os.path.join(work, "codes", source, f"{chrom}.npy")
    if not os.path.exists(f):
        return 0, 0
    c = np.load(f, mmap_mode="r")[a:b]
    return int(((c >= 1) & (c <= 5)).sum()), int(((c >= 2) & (c <= 5)).sum())


class IdentityNull:
    """Matched-span reference for homoeologue identity: divergence between each subgenome-s
    chromosome and its subgenome-t homoeologue over sliding blocks of the same span, genome-wide.
    p_identity = fraction of such blocks at least as similar (conservative: includes true
    conversion tracts)."""

    def __init__(self, chroms, data):
        self.chroms, self.data, self.cache = chroms, data, {}

    def p(self, s_i, t_i, span_bp, d_obs, n_obs):
        L = max(1, int(round(span_bp / PR.BIN_BP)))
        L = 1 << int(round(math.log2(L)))            # spans pooled on a doubling grid
        key = (s_i, t_i, L)
        if key not in self.cache:
            vals = []
            for c in self.chroms:
                if c["label_idx"] != s_i:
                    continue
                cv = np.concatenate([[0], np.cumsum(self.data[c["name"]]["ss_cov"][t_i], dtype=np.int64)])
                mm = np.concatenate([[0], np.cumsum(self.data[c["name"]]["ss_mm"][t_i], dtype=np.int64)])
                if cv.size <= L:
                    continue
                st = max(1, L // 2)
                i = np.arange(0, cv.size - L, st)
                c_, m_ = cv[i + L] - cv[i], mm[i + L] - mm[i]
                ok = c_ >= PR.IDENTITY_MIN_COVER * L * PR.BIN_BP
                vals.append(m_[ok] / c_[ok])
            v = np.sort(np.concatenate(vals)) if vals else np.zeros(0)
            self.cache[key] = v
        v = self.cache[key]
        if v.size == 0 or n_obs < PR.IDENTITY_MIN_COVER * span_bp:
            return 1.0
        return (1 + np.searchsorted(v, d_obs, side="right")) / (1 + v.size)


def calibrate(he, rho, rho_ac, span, S):
    """log posterior odds and QUAL of each candidate: serially calibrated origin evidence
    (llr / tau_n) + breakpoint prior at rate rho. Homoeologue identity is reported (p_identity) but
    not added: its genome-wide null contains the conversion tracts it would detect."""
    for h in he:
        n = h["n_windows"]
        tn = tau_n(rho_ac, n)
        nbp = int(h["bp_left_ci"] != "terminal") + int(h["bp_right_ci"] != "terminal")
        prior = nbp * math.log(max(rho * span[h["label_idx"]] / (S - 1), 1e-300))
        lo = h["llr"] / tn + prior
        post = 1.0 / (1.0 + math.exp(-lo)) if lo > -700 else 0.0
        h["tau_n"] = tn
        h["log_odds"] = lo
        h["qual"] = round(min(PR.QUAL_CAP, -10 * math.log10(max(1 - post, 10 ** (-PR.QUAL_CAP / 10)))), 1)


# ----------------------------------------------------------------------------- main
def _segments(mchroms, chroms, wins, posts, binned_list, slogp):
    segs = {}
    for m, c, w, post, binned in zip(mchroms, chroms, wins, posts, binned_list):
        sg = model.segments(m, w, post, binned, slogp)
        if not sg:
            L = c["length"]
            sg = [{"state": c["label_idx"], "tile_start": 0, "tile_end": L, "qual": 0.0, "post_mean": 0.0,
                   "nosites": True, "site_a": 0, "site_b": 0, "llr": 0.0, "win_a": 0, "win_b": 0,
                   "bp_left": (0, 0, 0, 0), "bp_right": (L, L, L, L),
                   "terminal_left": True, "terminal_right": True}]
        segs[c["name"]] = sg
    return segs


def _candidates(segs, homsegs, scale, ctx):
    """Non-native segments of one segmentation scale -> candidate HE records with evidence.
    homsegs: the fine-scale segmentation, used to type the homoeologue's origin."""
    (mchroms, chroms, labels, slogp, gaps, work, ref_sources, ref_aln, ss_aln, inull, base, name2chrom) = ctx
    he = []
    for m, c in zip(mchroms, chroms):
        s_i = c["label_idx"]
        pid = m["site_pid"]
        # contiguous blocks of non-native segments (any donor): their native flanks vouch for synteny
        block, run = {}, []
        for k, seg in enumerate(segs[c["name"]] + [None]):
            if seg is not None and seg["state"] != s_i and not seg.get("nosites"):
                run.append(k)
                continue
            for j in run:
                block[j] = (segs[c["name"]][run[0]]["tile_start"], segs[c["name"]][run[-1]]["tile_end"])
            run = []
        for k, seg in enumerate(segs[c["name"]]):
            t_i = seg["state"]
            if t_i == s_i or seg.get("nosites"):
                continue
            a, b = seg["tile_start"], seg["tile_end"]
            sa, sb = seg["site_a"], seg["site_b"]
            lr = slogp[t_i][pid[sa:sb]] - slogp[s_i][pid[sa:sb]]
            st = [code_stats(work, src, c["name"], a, b) for src in ref_sources]
            div = np.array([mm / cv if cv else np.inf for cv, mm in st])
            best = int(np.argmin(div))
            syn, tgt = syntenic(ref_aln[ref_sources[best]].get(c["name"]), a, b, block[k][0], block[k][1])
            pr = project(ss_aln[labels[t_i]].get(c["name"]), a, b)
            hom_type, hom_loc, frac_s = "unpaired", "", float("nan")
            ssc, ssm = code_stats(work, f"ss_{labels[t_i]}", c["name"], a, b)
            if pr is not None and pr["bp"] >= PR.HOMOEOLOG_MIN_COVER * (b - a):
                y, ylo, yhi = pr["target"], pr["lo"], pr["hi"]
                hom_loc = f"{y}:{ylo}-{yhi}"
                ysegs = homsegs.get(y)
                if ysegs is not None and name2chrom[y]["label_idx"] == t_i:
                    bp_s = bp_t = 0
                    for ys in ysegs:
                        ov = min(yhi, ys["tile_end"]) - max(ylo, ys["tile_start"])
                        if ov > 0:
                            if ys["state"] == s_i:
                                bp_s += ov
                            elif ys["state"] == t_i:
                                bp_t += ov
                    if bp_s + bp_t:
                        frac_s = bp_s / (bp_s + bp_t)
                        hom_type = "reciprocal" if frac_s >= 0.5 else "non-reciprocal"
            flags = []
            for key, term in (("bp_left", "terminal_left"), ("bp_right", "terminal_right")):
                lo, hi = seg[key][2], seg[key][3]
                if not seg[term] and any(g0 < hi and g1 > lo for g0, g1 in gaps.get(c["name"], [])):
                    flags.append("gap_at_breakpoint")
                    break
            if seg["terminal_left"] or seg["terminal_right"]:
                flags.append("terminal")
            ss_d = ssm / ssc if ssc else float("nan")
            # a non-reciprocal exchange copies the donor homoeologue, so the two should be closer than
            # homoeologues are genome-wide; if not, the reference signal has another source (a copy
            # missing from the assembly, a stand-in reference, paralogy)
            if (hom_type == "non-reciprocal" and ssc >= PR.HOMOEOLOG_MIN_COVER * (b - a)
                    and ss_d >= base[(s_i, t_i)]):
                flags.append("homoeolog_not_closer")
            p_id = inull.p(s_i, t_i, b - a, ss_d, ssc) if hom_type == "non-reciprocal" and ssc else 1.0
            he.append({"label_idx": s_i, "donor_idx": t_i, "p_identity": p_id, "scale": scale,
                       "chrom": c["name"], "label": labels[s_i], "donor": labels[t_i],
                       "start": a, "end": b, "length": b - a, "n_windows": seg["win_b"] - seg["win_a"],
                       "bp_left_ci": "terminal" if seg["terminal_left"] else f"{seg['bp_left'][2]}-{seg['bp_left'][3]}",
                       "bp_right_ci": "terminal" if seg["terminal_right"] else f"{seg['bp_right'][2]}-{seg['bp_right'][3]}",
                       "qual": round(seg["qual"], 1), "post_mean": seg["post_mean"], "llr": seg["llr"],
                       "n_sites": int(sb - sa), "n_support": int((lr > 0).sum()), "n_against": int((lr < 0).sum()),
                       "closest_ref": best, "ref_div": [float(x) for x in div],
                       "syntenic": syn, "ref_target": tgt or "",
                       "homoeolog": hom_loc, "type": hom_type, "homoeolog_frac_recipient_origin": frac_s,
                       "ss_div": float(ssm / ssc) if ssc else float("nan"), "ss_cover": float(ssc / max(b - a, 1)),
                       "ss_div_baseline": base[(s_i, t_i)],
                       "flags": ",".join(flags) or "."})
    return he


def label_coherence(chroms, segs, S):
    """Share of each label's evaluable sequence painted with that label's own origin.

    The model assumes every label is a subgenome of one origin, with exchanges a minority. Below one
    half, most of the label is of another origin: the labels are not subgenomes (e.g. phased
    haplotypes numbered independently per chromosome, or subgenomes the input mis-assigned), and
    exchange calls involving that label are not interpretable."""
    own, tot = np.zeros(S), np.zeros(S)
    for c in chroms:
        for s in segs[c["name"]]:
            if not s.get("nosites"):
                L = s["tile_end"] - s["tile_start"]
                tot[c["label_idx"]] += L
                own[c["label_idx"]] += L * (s["state"] == c["label_idx"])
    return [float(own[i] / tot[i]) if tot[i] else float("nan") for i in range(S)]


def _filters(h, resolvable, incoherent=()):
    filt = []
    if h["qual"] < PR.QUAL_PASS:
        filt.append("LowQual")
    if not resolvable[(h["label_idx"], h["donor_idx"])]:
        filt.append("unresolvable_pair")
    if h["label_idx"] in incoherent or h["donor_idx"] in incoherent:
        filt.append("incoherent_label")
    if h["syntenic"] == "nonsyntenic":
        filt.append("nonsyntenic")
    return ";".join(filt) or "PASS"


def _select(cands):
    """Non-overlapping calls from candidates of all scales.

    Candidates passing every filter come first, longest first: when a long block and a fragment of
    it both meet the evidence standard, the block describes the supported extent. The remaining
    candidates follow by calibrated log-odds. Each call keeps, as its core, the best-scoring
    candidate (any scale) of the same donor inside it."""
    order = sorted(cands, key=lambda h: (h["filter"] != "PASS",
                                         -h["length"] if h["filter"] == "PASS" else -h["log_odds"]))
    taken, out = {}, []
    for h in order:
        iv = taken.setdefault(h["chrom"], [])
        if any(h["start"] < e and h["end"] > s for s, e in iv):
            continue
        iv.append((h["start"], h["end"]))
        out.append(h)
    by_chrom = {}
    for c in cands:
        by_chrom.setdefault(c["chrom"], []).append(c)
    for h in out:
        inside = [c for c in by_chrom[h["chrom"]] if c["donor"] == h["donor"]
                  and c["start"] >= h["start"] and c["end"] <= h["end"]]
        core = max(inside, key=lambda c: c["log_odds"]) if inside else h
        h["core"] = f"{core['start']}-{core['end']}"
        h["core_qual"] = core["qual"]
    return sorted(out, key=lambda h: (h["chrom"], h["start"]))


def make_calls(mchroms, chroms, labels, fit, wins, slogp, data, gaps, work, ref_sources, log,
               posts, rho_start=None, tables=None):
    """Segmentations at one or more scales -> scored, non-overlapping HE calls.

    posts: [(scale name, per-chromosome posteriors)]; the first is the fine (untempered) scale, whose
    segments are returned for painting and used to type homoeologues. Returns (segs, he, res,
    rho_cal, tables); `tables` (loaded alignment tables) can be passed back in to avoid reloading."""
    S = len(labels)
    inull = IdentityNull(chroms, data)
    segs_by = [(name, _segments(mchroms, chroms, wins, p, fit["binned"], slogp)) for name, p in posts]
    fine = segs_by[0][1]
    coherence = label_coherence(chroms, fine, S)
    incoherent = {i for i in range(S) if coherence[i] < 0.5}
    if incoherent and tables is None:  # once per run (make_calls runs per scale set)
        log("calls: WARNING: " + ", ".join(f"label {labels[i]} is {100 * coherence[i]:.0f}% its own origin"
                                          for i in sorted(incoherent))
            + "; not a coherent subgenome, so calls involving it are filtered (incoherent_label)")
    log("calls: synteny and homoeologue checks", detail=True)
    if tables is None:
        tables = ({s: extract.load_aln(work, s) for s in ref_sources},
                  {lab: extract.load_aln(work, f"ss_{lab}") for lab in labels})
    ref_aln, ss_aln = tables
    name2chrom = {c["name"]: c for c in chroms}
    base = {}
    for s_i in range(S):
        for t_i in range(S):
            if s_i != t_i:
                num = sum(int(data[c["name"]]["ss_mm"][t_i].sum()) for c in chroms if c["label_idx"] == s_i)
                den = sum(int(data[c["name"]]["ss_cov"][t_i].sum()) for c in chroms if c["label_idx"] == s_i)
                base[(s_i, t_i)] = num / den if den else float("nan")
    ctx = (mchroms, chroms, labels, slogp, gaps, work, ref_sources, ref_aln, ss_aln, inull, base, name2chrom)
    cands = []
    for name, sg in segs_by:
        cands += _candidates(sg, fine, name, ctx)

    # empirical-Bayes breakpoint rate: rate implied by the PASS calls, iterated to a fixed point.
    # Several self-consistent rates can exist (a lenient prior lets more calls pass, which justify the
    # lenient prior); starting from the lower bound reaches the smallest one, i.e. the most
    # conservative rate the calls themselves support.
    span = {}
    for s_i in range(S):
        sp = [w["end"] - w["start"] for w in wins if w["label"] == s_i and w["end"].size]
        span[s_i] = float(np.median(np.concatenate(sp))) if sp else 1.0
    evaluated = float(sum((w["end"] - w["start"]).sum() for w in wins))
    rho = rho_start if rho_start is not None else PR.RHO_MIN
    for _ in range(200):
        res = model.resolution(fit, wins, S, rho, fit["rho_ac"])
        if fit.get("no_contrast"):
            for r in res:
                r["resolvable"] = False
        resolvable = {(r["recipient"], r["donor"]): r["resolvable"] for r in res}
        calibrate(cands, rho, fit["rho_ac"], span, S)
        for h in cands:
            h["filter"] = _filters(h, resolvable, incoherent)
        he = _select(cands)
        nbp = sum(int(h["bp_left_ci"] != "terminal") + int(h["bp_right_ci"] != "terminal")
                  for h in he if h["filter"] == "PASS")
        new = min(max(nbp / evaluated, PR.RHO_MIN), PR.RHO_MAX)
        if abs(math.log(new / rho)) < 0.01:
            rho = new
            break
        rho = new
    for r in res:  # genome-wide homoeologue divergence, also when there are no calls
        r["ss_div_baseline"] = base[(r["recipient"], r["donor"])]
        r["recipient_own_origin"] = coherence[r["recipient"]]
    log(f"calls: calibrated breakpoint rate {rho:.3g}/bp; serial correlation tau(5 windows)="
        f"{tau_n(fit['rho_ac'], 5):.1f}, tau(100)={tau_n(fit['rho_ac'], 100):.1f}; "
        + ", ".join(f"scale {n}: {sum(h['scale'] == n and h['filter'] == 'PASS' for h in he)} PASS"
                    for n, _ in segs_by))
    return fine, he, res, rho, tables


def nocall_regions(chroms, wins, data, refs, min_len):
    """Stretches covered by no window (too sparsely aligned) at least one typical window long."""
    out = {}
    B = PR.BIN_BP
    for c, w in zip(chroms, wins):
        L = c["length"]
        edges_s = np.concatenate([w["end"], [L]]) if w["start"].size else np.array([L])
        edges_a = np.concatenate([[0], w["end"]]) if w["start"].size else np.array([0])
        starts = np.concatenate([w["start"], [L]]) if w["start"].size else np.array([L])
        regs = []
        d = data[c["name"]]
        for a, b in zip(edges_a, starts):
            a, b = int(a), int(b)
            if b - a < min_len:
                continue
            b0, b1 = a // B, max(a // B + 1, (b + B - 1) // B)
            span = (b1 - b0) * B
            frac = d["cov"][:, b0:b1].sum(axis=1) / span
            if frac.max() < 0.05:
                reason = "no_alignment"
            else:
                weak = [refs[k][0] for k in np.flatnonzero(frac < 0.5 * frac.max())]
                reason = "sparse_alignment" + (f"({','.join(weak)})" if weak else "")
            regs.append((a, b, reason))
        out[c["name"]] = regs
    return out


def chromosome_summary(chroms, labels, segs, he, nocall, data):
    """Per chromosome: origin composition, no-call share, PASS calls, and how much of it aligns
    (collinear chains) to each other subgenome. Low homoeologue cover marks sequence with no
    counterpart in the assembly: lost homoeologues (LOH in trio-binned or collapsed assemblies),
    aneuploidy, rearranged or unassembled arms. No call can be made there."""
    rows = []
    for c in chroms:
        L = c["length"]
        ssc = data[c["name"]]["ss_cov"]
        hcov = {f"pct_homoeolog_{lab}": ("." if i == c["label_idx"] else round(100 * float(ssc[i].sum()) / L, 2))
                for i, lab in enumerate(labels)}
        bp = np.zeros(len(labels))
        for s in segs[c["name"]]:
            if not s.get("nosites"):
                bp[s["state"]] += s["tile_end"] - s["tile_start"]
        nc = sum(b - a for a, b, _ in nocall[c["name"]])
        maj = labels[int(np.argmax(bp))] if bp.sum() else "."
        hp = [h for h in he if h["chrom"] == c["name"] and h["filter"] == "PASS"]
        rows.append({"chrom": c["name"], "label": c["label"], "length": L, "majority_origin": maj,
                     **{f"pct_{lab}": round(100 * bp[i] / L, 2) for i, lab in enumerate(labels)},
                     "pct_nocall": round(100 * nc / L, 2), **hcov,
                     "n_he_pass": len(hp), "he_bp_pass": sum(h["length"] for h in hp),
                     "label_conflict": maj not in (".", c["label"])})
    return rows
