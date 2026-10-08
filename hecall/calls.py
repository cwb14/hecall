"""Segments -> HE calls with quality, independent evidence and artifact flags."""
import os

import numpy as np

from . import extract, model, sites
from . import params as PR


# ----------------------------------------------------------------------------- helpers
def project(aln, a, b):
    """Project query interval [a,b) through alignments; returns (target, t_lo, t_hi, overlap_bp)
    for the target with the largest overlap, or None."""
    if aln is None:
        return None
    sel = np.flatnonzero((aln["qs"] < b) & (aln["qe"] > a))
    if sel.size == 0:
        return None
    best = {}
    for i in sel:
        x0, x1 = max(a, aln["qs"][i]), min(b, aln["qe"][i])
        if aln["fwd"][i]:
            t0, t1 = aln["ts"][i] + (x0 - aln["qs"][i]), aln["ts"][i] + (x1 - aln["qs"][i])
        else:
            t0, t1 = aln["te"][i] - (x1 - aln["qs"][i]), aln["te"][i] - (x0 - aln["qs"][i])
        t = aln["t"][i]
        o = best.setdefault(t, [0, np.inf, -np.inf])
        o[0] += x1 - x0
        o[1] = min(o[1], t0)
        o[2] = max(o[2], t1)
    t, (ov, lo, hi) = max(best.items(), key=lambda kv: kv[1][0])
    return t, int(lo), int(hi), int(ov)


def flank_point(aln, a, b, side):
    """Target (chrom, pos, query pos) just outside [a,b): through the alignment reaching closest to
    a from the left (side='L') or to b from the right ('R'). An alignment spanning the boundary
    counts (it then vouches for collinearity on its own)."""
    if aln is None:
        return None
    if side == "L":
        sel = np.flatnonzero(aln["qs"] < a)
        if sel.size == 0:
            return None
        reach = np.minimum(aln["qe"][sel], a)
        i = sel[np.argmax(reach)]
        q = int(min(aln["qe"][i], a)) - 1
    else:
        sel = np.flatnonzero(aln["qe"] > b)
        if sel.size == 0:
            return None
        reach = np.maximum(aln["qs"][sel], b)
        i = sel[np.argmin(reach)]
        q = int(max(aln["qs"][i], b))
    t = aln["ts"][i] + (q - aln["qs"][i]) if aln["fwd"][i] else aln["te"][i] - (q - aln["qs"][i]) - 1
    return aln["t"][i], int(t), int(q)


def syntenic(aln, a, b):
    """Is the reference alignment covering [a,b) collinear with its flanks? True/False/None."""
    pr = project(aln, a, b)
    if pr is None:
        return None, None
    t, lo, hi, _ = pr
    mid_q, mid_t = (a + b) // 2, (lo + hi) // 2
    verdicts = []
    for side in ("L", "R"):
        fp = flank_point(aln, a, b, side)
        if fp is None or abs(fp[2] - mid_q) > PR.SYNTENY_FLANK_MAX_BP:
            continue
        ft, fpos, fq = fp
        ok = ft == t and abs(fpos - mid_t) <= PR.SYNTENY_SLACK_BP + PR.SYNTENY_FOLD * abs(fq - mid_q)
        verdicts.append(ok)
    target = f"{t}:{lo}-{hi}"
    if not verdicts:
        return None, target
    return all(verdicts), target


def window_rate(num, den, a, b):
    W = PR.WINDOW_BP
    w0, w1 = a // W, max(a // W + 1, (b + W - 1) // W)
    n, d = num[..., w0:w1].sum(axis=-1), den[..., w0:w1].sum(axis=-1)
    return n, d


def code_stats(work, source, chrom, a, b):
    """(covered SNV-comparable bases, mismatches) of one source over [a,b), from the base codes."""
    f = os.path.join(work, "codes", source, f"{chrom}.npy")
    if not os.path.exists(f):
        return 0, 0
    c = np.load(f, mmap_mode="r")[a:b]
    return int(((c >= 1) & (c <= 5)).sum()), int(((c >= 2) & (c <= 5)).sum())


def pattern_names(allkeys, remap, n_pat, K):
    names = [[] for _ in range(n_pat)]
    for k, r in zip(allkeys, remap):
        names[r].append(sites.key_to_pattern(int(k), K))
    return ["|".join(n) if len(n) == 1 else f"rare({len(n)})" for n in names]


# ----------------------------------------------------------------------------- main
def make_calls(mchroms, chroms, labels, fit, res, data, gaps, work, ref_sources, log):
    S = len(labels)
    logp = fit["logp"]
    resolvable = {(r["recipient"], r["donor"]): r["resolvable"] for r in res}
    segs = {}
    for m, c, binned, post in zip(mchroms, chroms, fit["binned"], fit["post"]):
        sg = model.segments(m, binned, post, logp)
        L = c["length"]
        for i, s in enumerate(sg):
            s["tile_start"] = 0 if i == 0 else (s["bp_left"][0] + s["bp_left"][1]) // 2
            s["tile_end"] = L if i == len(sg) - 1 else (s["bp_right"][0] + s["bp_right"][1]) // 2
            s["qual"] = model.qual(s["post_mean"])
        if not sg:
            sg = [{"state": c["label_idx"], "tile_start": 0, "tile_end": L, "qual": 0.0, "post_mean": 0.0,
                   "nosites": True, "site_a": 0, "site_b": 0, "start": None, "end": None,
                   "bp_left": (0, 0, 0, 0), "bp_right": (L, L, L, L),
                   "terminal_left": True, "terminal_right": True}]
        segs[c["name"]] = sg

    log("calls: synteny and homoeologue checks", detail=True)
    ref_aln = {s: extract.load_aln(work, s) for s in ref_sources}
    ss_aln = {lab: extract.load_aln(work, f"ss_{lab}") for lab in labels}
    name2chrom = {c["name"]: c for c in chroms}

    # genome-wide homoeologue divergence baseline per (label, other label): native windows only
    base = {}
    for s_i in range(S):
        for t_i in range(S):
            if s_i == t_i:
                continue
            num = den = 0
            for c in chroms:
                if c["label_idx"] != s_i:
                    continue
                d = data[c["name"]]
                num += int(d["ss_mm"][t_i].sum())
                den += int(d["ss_cov"][t_i].sum())
            base[(s_i, t_i)] = num / den if den else float("nan")

    he = []
    for m, c in zip(mchroms, chroms):
        s_i = c["label_idx"]
        d = data[c["name"]]
        pos, pid = m["site_pos"], m["site_pid"]
        for seg in segs[c["name"]]:
            t_i = seg["state"]
            if t_i == s_i or seg.get("nosites"):
                continue
            a, b = seg["tile_start"], seg["tile_end"]
            sa, sb = seg["site_a"], seg["site_b"]
            lr = logp[t_i][pid[sa:sb]] - logp[s_i][pid[sa:sb]]
            # closest reference inside the segment (lowest SNV divergence)
            st = [code_stats(work, src, c["name"], a, b) for src in ref_sources]
            div = np.array([mm / cv if cv else np.inf for cv, mm in st])
            best = int(np.argmin(div))
            syn, tgt = syntenic(ref_aln[ref_sources[best]].get(c["name"]), a, b)
            # homoeologue
            pr = project(ss_aln[labels[t_i]].get(c["name"]), a, b)
            hom_type, hom_loc, frac_s = "unpaired", "", float("nan")
            ssc, ssm = code_stats(work, f"ss_{labels[t_i]}", c["name"], a, b)
            ss_cov_frac = ssc / max(b - a, 1)
            if pr is not None and pr[3] >= PR.HOMOEOLOG_MIN_COVER * (b - a):
                y, ylo, yhi, _ = pr
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
                    tot = bp_s + bp_t
                    if tot:
                        frac_s = bp_s / tot
                        hom_type = "reciprocal" if frac_s >= 0.5 else "non-reciprocal"
            # flags
            flags = []
            for key in ("bp_left", "bp_right"):
                lo, hi = seg[key][2], seg[key][3]
                if any(g0 < hi and g1 > lo for g0, g1 in gaps.get(c["name"], [])) and \
                        not seg["terminal_" + key.split("_")[1]]:
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
                "start": a, "end": b, "length": b - a,
                "core_start": seg["start"], "core_end": seg["end"],
                "bp_left_ci": "terminal" if seg["terminal_left"] else f"{seg['bp_left'][2]}-{seg['bp_left'][3]}",
                "bp_right_ci": "terminal" if seg["terminal_right"] else f"{seg['bp_right'][2]}-{seg['bp_right'][3]}",
                "qual": round(seg["qual"], 1), "post_mean": seg["post_mean"],
                "llr": float(lr.sum()), "n_sites": int(sb - sa), "n_support": int((lr > 0).sum()),
                "n_against": int((lr < 0).sum()),
                "closest_ref": best, "ref_div": [float(x) for x in div],
                "syntenic": syn, "ref_target": tgt or "",
                "homoeolog": hom_loc, "type": hom_type, "homoeolog_frac_recipient_origin": frac_s,
                "ss_div": float(ssm / ssc) if ssc else float("nan"), "ss_cover": float(ss_cov_frac),
                "ss_div_baseline": base[(s_i, t_i)],
                "flags": ",".join(flags) or ".", "filter": ";".join(filt) or "PASS"})
    return segs, he


def nocall_regions(mchroms, chroms, data, res, refs):
    """Stretches without informative sites longer than the shortest detectable HE."""
    lmin = min((r["min_he_bp"] for r in res), default=float("inf"))
    lmin = max(lmin, PR.BIN_BP)
    out = {}
    W = PR.WINDOW_BP
    for m, c in zip(mchroms, chroms):
        L = c["length"]
        p = m["site_pos"]
        edges = np.concatenate([[0], p, [L]])
        gl = np.diff(edges)
        idx = np.flatnonzero(gl > lmin)
        d = data[c["name"]]
        regs = []
        for i in idx:
            a = int(edges[i]) + (1 if i > 0 else 0)
            b = int(edges[i + 1])
            w0, w1 = a // W, max(a // W + 1, (b + W - 1) // W)
            span = (w1 - w0) * W
            f2 = d["cov2"][w0:w1].sum() / span
            f1 = d["cov"][:, w0:w1].sum(axis=1).max() / span
            reason = "no_differences" if f2 >= 0.5 else ("one_reference" if f1 >= 0.5 else "no_alignment")
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
        ev = bp.sum()
        maj = labels[int(np.argmax(bp))] if ev else "."
        hp = [h for h in he if h["chrom"] == c["name"] and h["filter"] == "PASS"]
        rows.append({"chrom": c["name"], "label": c["label"], "length": L,
                     "majority_origin": maj,
                     **{f"pct_{lab}": round(100 * bp[i] / L, 2) for i, lab in enumerate(labels)},
                     "pct_nocall": round(100 * nc / L, 2),
                     "n_he_pass": len(hp), "he_bp_pass": sum(h["length"] for h in hp),
                     "label_conflict": maj not in (".", c["label"])})
    return rows
