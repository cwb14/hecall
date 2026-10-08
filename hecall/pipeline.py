"""Stage orchestration with checkpoints in OUT/work/.

prepare -> align (FastGA + ALNchain) -> extract (codes per source) -> sites (per chromosome)
-> model (profiles + HMM) -> calls (segments, evidence, flags) -> outputs + figure
"""
import json
import math
import os
import shutil
import time
from multiprocessing import Pool

import numpy as np

from . import align, extract, inputs, model, sites, windows, calls, output, plot
from . import params as PR


class Log:
    def __init__(self, path, verbose):
        self.fh = open(path, "a")
        self.verbose = verbose
        self.t0 = time.time()

    def __call__(self, msg, detail=False):
        line = f"[hecall {time.time() - self.t0:7.0f}s] {msg}"
        self.fh.write(line + "\n")
        self.fh.flush()
        if self.verbose or not detail:
            print(line, flush=True)


# ----------------------------------------------------------------------------- sites stage (worker)
def _sites_worker(job):
    workdir, chrom, ref_sources, ss_sources = job
    out = os.path.join(workdir, "sites", f"{chrom['name']}.npz")
    if os.path.exists(out):
        return out
    L = chrom["length"]
    mat = np.vstack([extract.load_codes(workdir, s, chrom["name"], L) for s in ref_sources])
    pos, key = sites.informative_sites(mat, PR.MIN_REFS_COVERED)
    B = PR.BIN_BP
    nb = (L + B - 1) // B

    def block_sum(mask):
        pad = np.zeros(nb * B, dtype=bool)
        pad[:L] = mask
        return pad.reshape(nb, B).sum(axis=1).astype(np.int32)

    cov = np.vstack([block_sum((m >= 1) & (m <= 5)) for m in mat])
    mm = np.vstack([block_sum((m >= 2) & (m <= 5)) for m in mat])
    del mat
    ss_cov, ss_mm = [], []
    for s in ss_sources:
        if s is None:
            ss_cov.append(np.zeros(nb, dtype=np.int32))
            ss_mm.append(np.zeros(nb, dtype=np.int32))
            continue
        c = extract.load_codes(workdir, s, chrom["name"], L)
        ss_cov.append(block_sum((c >= 1) & (c <= 5)))
        ss_mm.append(block_sum((c >= 2) & (c <= 5)))
    tmp = out + ".tmp.npz"
    np.savez(tmp, pos=pos, key=key, cov=cov, mm=mm, ss_cov=np.vstack(ss_cov), ss_mm=np.vstack(ss_mm))
    os.replace(tmp, out)
    return out


def _extract_worker(job):
    chain, source, chrom_len, workdir, tools, threads, target_set = job
    extract.extract_source(chain, source, chrom_len, workdir, tools, threads, target_set)
    return source


# ----------------------------------------------------------------------------- main
def run(subgenome_args, ref_args, outdir, threads=8, verbose=False, keep_work=True):
    os.makedirs(outdir, exist_ok=True)
    work = os.path.join(outdir, "work")
    os.makedirs(work, exist_ok=True)
    log = Log(os.path.join(outdir, "hecall.log"), verbose)
    tools = align.find_tools()

    # ---- prepare
    subs = [inputs.parse_labelled(a, "subgenome") for a in subgenome_args]
    refs = [inputs.parse_labelled(a, "reference") for a in ref_args]
    if len(refs) < PR.MIN_REFS_COVERED:
        raise SystemExit(f"hecall: need at least {PR.MIN_REFS_COVERED} reference genomes (-r)")
    names = [r[0] for r in refs]
    if len(set(names)) != len(names):
        raise SystemExit("hecall: reference names must be unique (use -r NAME=path)")
    chroms, labels = inputs.chromosome_table(subs)
    S, K = len(labels), len(refs)
    log(f"start: {len(chroms)} polyploid chromosomes in {S} subgenomes ({', '.join(labels)}); "
        f"{K} references ({', '.join(names)}); {threads} threads")
    json.dump({"subgenomes": subs, "references": refs, "labels": labels, "threads": threads},
              open(os.path.join(outdir, "run.json"), "w"), indent=1)
    alog = os.path.join(work, "align.log")
    poly = os.path.join(work, "poly")
    if not os.path.exists(poly + ".1gdb"):
        poly_fa = inputs.concat_fastas([p for _, p in subs],
                                       os.path.join(work, "poly.fa.gz" if subs[0][1].endswith(".gz") else "poly.fa"))
        align.make_gdb(poly_fa, poly, tools, alog)
        os.remove(poly_fa)  # the GDB holds the sequence
    ref_gdbs = [align.make_gdb(p, os.path.join(work, f"ref{i}"), tools, alog) for i, (_, p) in enumerate(refs)]
    gaps = inputs.gdb_gaps(poly, tools)
    chrom_len = {c["name"]: c["length"] for c in chroms}

    # ---- align
    tmp = os.path.join(work, "tmp")
    os.makedirs(tmp, exist_ok=True)
    chains = []
    need = [i for i in range(K) if not os.path.exists(os.path.join(work, f"aln_ref{i}.chain.1aln"))]
    need_self = not os.path.exists(os.path.join(work, "aln_self.chain.1aln"))
    if need or need_self:
        if not os.path.exists(poly + ".gix"):
            align.run([tools["GIXmake"], f"-T{threads}", f"-P{tmp}", poly + ".1gdb"], alog, tools)
    for i in range(K):
        log(f"align: polyploid vs {names[i]}", detail=i not in need)
        chains.append(align.align_chain(poly, ref_gdbs[i], os.path.join(work, f"aln_ref{i}"), threads, tmp, tools, alog))
    log("align: polyploid vs itself (homoeologues)", detail=not need_self)
    self_chain = align.align_chain(poly, None, os.path.join(work, "aln_self"), threads, tmp, tools, alog)
    if os.path.exists(poly + ".gix"):
        align.run([tools["GIXrm"], poly + ".gix"], alog, tools)

    # ---- extract
    ref_sources = [f"ref{i}" for i in range(K)]
    jobs = [(chains[i], ref_sources[i], chrom_len, work, tools, 2, None) for i in range(K)]
    for lab in labels:
        tset = {c["name"] for c in chroms if c["label"] == lab}
        jobs.append((self_chain, f"ss_{lab}", chrom_len, work, tools, 2, tset))
    log("extract: per-base comparison with each reference and homoeologue")
    with Pool(min(threads, len(jobs))) as pool:
        for src in pool.imap_unordered(_extract_worker, jobs):
            log(f"  extracted {src}", detail=True)

    # ---- sites
    os.makedirs(os.path.join(work, "sites"), exist_ok=True)
    sjobs = []
    for c in chroms:
        ss = [None if lab == c["label"] else f"ss_{lab}" for lab in labels]
        sjobs.append((work, c, ref_sources, ss))
    log("sites: informative sites per chromosome")
    sj = sorted(sjobs, key=lambda j: -j[1]["length"])
    with Pool(max(1, min(threads, len(sj)))) as pool:
        for f in pool.imap_unordered(_sites_worker, sj):
            pass
    data = {c["name"]: dict(np.load(os.path.join(work, "sites", f"{c['name']}.npz"))) for c in chroms}

    # ---- windows
    dtyp = windows.typical_closest_divergence(data)
    if not dtyp > 0:
        raise SystemExit("hecall: no divergence between the polyploid and its closest references; nothing to tell apart")
    cmin = int(math.ceil(PR.WINDOW_MISMATCHES / dtyp))
    span0 = [np.diff(np.vstack(windows.build(data[c["name"]]["cov"], data[c["name"]]["mm"], cmin, 1 << 60)), axis=0)[0]
             for c in chroms]
    span0 = np.concatenate(span0)
    maxspan = int(PR.WINDOW_MAX_FACTOR * np.median(span0)) if span0.size else 1
    wins = []
    for c in chroms:
        d = data[c["name"]]
        s0, e0 = windows.build(d["cov"], d["mm"], cmin, maxspan)
        cov = windows.window_sums(d["cov"], s0, e0)
        mm = windows.window_sums(d["mm"], s0, e0)
        start = s0 * PR.BIN_BP
        end = np.minimum(e0 * PR.BIN_BP, c["length"])
        ctr = (start + end) / 2.0
        wins.append({"name": c["name"], "length": c["length"], "label": c["label_idx"], "start": start,
                     "end": end, "cov": cov, "mm": mm, "d": np.diff(ctr, prepend=ctr[:1]),
                     "Xraw": model.features(cov, mm) if s0.size else np.zeros((0, K - 1))})
    nwin = sum(w["start"].size for w in wins)
    if nwin == 0:
        raise SystemExit("hecall: no evaluable windows (references do not align to the polyploid)")
    Xall = np.vstack([w["Xraw"] for w in wins])
    proj = model.lda(Xall, np.concatenate([np.full(w["start"].size, w["label"]) for w in wins]), S)
    for w in wins:
        w["X"] = w["Xraw"] @ proj
    log(f"windows: {nwin:,} windows of ~{PR.WINDOW_MISMATCHES} expected mismatches to the closest reference "
        f"(typical divergence {100 * dtyp:.2g}%, >= {cmin:,} compared bases per reference; median span "
        f"{output.fmt_bp(np.median(span0) * PR.BIN_BP)}); {proj.shape[1]}-dimensional origin space")

    # ---- model
    fit = model.run(wins, S, verbose=lambda m: log(m, detail=True))
    log(f"model: segmentation breakpoint rate {fit['rho']:.3g}/bp")
    allkeys = np.unique(np.concatenate([d["key"] for d in data.values()]))
    mchroms = []
    for c, w, a in zip(chroms, wins, fit["assign"]):
        d = data[c["name"]]
        pos = d["pos"].astype(np.int64)
        wi = np.searchsorted(w["start"], pos, side="right") - 1          # window holding each site
        inside = (wi >= 0) & (pos < w["end"][np.clip(wi, 0, None)]) if w["start"].size else np.zeros(pos.size, bool)
        st = np.full(pos.size, c["label_idx"], dtype=np.int64)
        st[inside] = a[wi[inside]]
        mchroms.append({"name": c["name"], "length": c["length"], "label": c["label_idx"], "site_pos": pos,
                        "site_pid": np.searchsorted(allkeys, d["key"]).astype(np.int64), "site_state": st})
    slogp = model.site_profiles(mchroms, S, allkeys.size)
    log(f"sites: {sum(m['site_pos'].size for m in mchroms):,} diagnostic sites place breakpoints")

    # ---- calls
    segs, he, res, rho_c = calls.make_calls(mchroms, chroms, labels, fit, wins, slogp, data, gaps, work, ref_sources, log)
    fit["rho_cal"] = rho_c
    log("model: " + "; ".join(f"{labels[r['donor']]}->{labels[r['recipient']]} min HE {output.fmt_bp(r['min_he_bp'])}"
                             + ("" if r["resolvable"] else " (unresolvable)") for r in res))
    nocall = calls.nocall_regions(chroms, wins, data, refs, np.median(span0) * PR.BIN_BP)
    chrom_rows = calls.chromosome_summary(chroms, labels, segs, he, nocall)

    # ---- outputs
    pattern_names = [sites.key_to_pattern(int(k), K) for k in allkeys]
    output.write_all(outdir, labels, refs, chroms, segs, he, nocall, chrom_rows, fit, wins, res,
                     slogp, pattern_names, data)
    plot.painting(outdir, labels, chroms, segs, he, nocall)
    npass = sum(1 for h in he if h["filter"] == "PASS")
    # per-base codes are K+S bytes per polyploid base: remove them (and their done markers, so a
    # rerun rebuilds them from the chained alignments); keep the small alignment tables
    codes = os.path.join(work, "codes")
    for src in ref_sources + [f"ss_{lab}" for lab in labels]:
        shutil.rmtree(os.path.join(codes, src), ignore_errors=True)
        if os.path.exists(os.path.join(codes, src + ".done")):
            os.remove(os.path.join(codes, src + ".done"))
    from . import __version__
    with open(os.path.join(outdir, "hecall.done"), "w") as fh:
        fh.write(f"hecall {__version__}\t{npass} PASS\t{len(he)} candidates\n")
    log(f"done: {npass} PASS homoeologous exchanges ({len(he)} candidate segments); outputs in {outdir}")
    return 0
