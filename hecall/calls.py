"""Segments -> HE calls with quality, independent evidence and artifact flags."""
import os

import numpy as np

from . import extract, model
from . import params as PR


# ----------------------------------------------------------------------------- helpers
def project(aln, a, b):
    """Project query interval [a,b) through alignments, robustly.

    Returns dict(target, lo, hi, t_med, q_med, bp) for the target chromosome holding most aligned
    bases, with lo/hi = 5th/95th and t_med = 50th aligned-base-weighted percentile of target
    positions (stray repeat alignments elsewhere on the target do not stretch the interval), or None."""
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
            "q_med": float(np.average(qm, weights=wm)), "bp": int(wm.sum())}


def syntenic(aln, a, b, block_a, block_b):
    """Is the reference sequence matching [a,b) collinear with the sequence matching the native
    flanks of the non-native block [block_a, block_b) that contains it?

    Flanks: max(block length, SYNTENY_FLANK_MIN_BP) (capped at SYNTENY_FLANK_MAX_BP) on each side of
    the block. True if collinear with at least one flank (a rearrangement between assemblies can
    sit on one side); False if no flank agrees (transposed repeat, translocation, misassembly);
    None if no flank aligns. Also returns the segment's core target interval."""
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
    return (any(verdicts) if verdicts else None), target


def code_stats(work, source, chrom, a, b):
    """(covered SNV-comparable bases, mismatches) of one source over [a,b), from the base codes."""
    f = os.path.join(work, "codes", source, f"{chrom}.npy")
    if not os.path.exists(f):
        return 0, 0
    c = np.load(f, mmap_mode="r")[a:b]
    return int(((c >= 1) & (c <= 5)).sum()), int(((c >= 2) & (c <= 5)).sum())


# ----------------------------------------------------------------------------- main
def make_calls(mchroms, chroms, labels, fit, wins, res, slogp, data, gaps, work, ref_sources, log):
    S = len(labels)
    resolvable = {(r["recipient"], r["donor"]): r["resolvable"] for r in res}
    segs = {}
    for m, c, w, post, binned in zip(mchroms, chroms, wins, fit["post"], fit["binned"]):
        sg = model.segments(m, w, post, binned, slogp)
        if not sg:
            L = c["length"]
            sg = [{"state": c["label_idx"], "tile_start": 0, "tile_end": L, "qual": 0.0, "post_mean": 0.0,
                   "nosites": True, "site_a": 0, "site_b": 0, "llr": 0.0, "win_a": 0, "win_b": 0,
                   "bp_left": (0, 0, 0, 0), "bp_right": (L, L, L, L),
                   "terminal_left": True, "terminal_right": True}]
        segs[c["name"]] = sg

    log("calls: synteny and homoeologue checks", detail=True)
    ref_aln = {s: extract.load_aln(work, s) for s in ref_sources}
    ss_aln = {lab: extract.load_aln(work, f"ss_{lab}") for lab in labels}
    name2chrom = {c["name"]: c for c in chroms}

    # genome-wide homoeologue divergence per (label, other label)
    base = {}
    for s_i in range(S):
        for t_i in range(S):
            if s_i != t_i:
                num = sum(int(data[c["name"]]["ss_mm"][t_i].sum()) for c in chroms if c["label_idx"] == s_i)
                den = sum(int(data[c["name"]]["ss_cov"][t_i].sum()) for c in chroms if c["label_idx"] == s_i)
                base[(s_i, t_i)] = num / den if den else float("nan")

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
                ysegs = segs.get(y)
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
            filt = []
            if seg["qual"] < PR.QUAL_PASS:
                filt.append("LowQual")
            if not resolvable[(s_i, t_i)]:
                filt.append("unresolvable_pair")
            if syn is False:
                filt.append("nonsyntenic")
            he.append({
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
                "flags": ",".join(flags) or ".", "filter": ";".join(filt) or "PASS"})
    return segs, he


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


def chromosome_summary(chroms, labels, segs, he, nocall):
    rows = []
    for c in chroms:
        L = c["length"]
        bp = np.zeros(len(labels))
        for s in segs[c["name"]]:
            if not s.get("nosites"):
                bp[s["state"]] += s["tile_end"] - s["tile_start"]
        nc = sum(b - a for a, b, _ in nocall[c["name"]])
        maj = labels[int(np.argmax(bp))] if bp.sum() else "."
        hp = [h for h in he if h["chrom"] == c["name"] and h["filter"] == "PASS"]
        rows.append({"chrom": c["name"], "label": c["label"], "length": L, "majority_origin": maj,
                     **{f"pct_{lab}": round(100 * bp[i] / L, 2) for i, lab in enumerate(labels)},
                     "pct_nocall": round(100 * nc / L, 2),
                     "n_he_pass": len(hp), "he_bp_pass": sum(h["length"] for h in hp),
                     "label_conflict": maj not in (".", c["label"])})
    return rows
