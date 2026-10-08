"""Stage orchestration with checkpoints in OUT/work/.

prepare -> align (FastGA + ALNchain) -> extract (codes per source) -> sites (per chromosome)
-> model (profiles + HMM) -> calls (segments, evidence, flags) -> outputs + figure
"""
import json
import os
import time
from multiprocessing import Pool

import numpy as np

from . import align, extract, inputs, model, sites, calls, output, plot
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
    W = PR.WINDOW_BP
    nw = (L + W - 1) // W

    def win_sum(mask):
        pad = np.zeros(nw * W, dtype=bool)
        pad[:L] = mask
        return pad.reshape(nw, W).sum(axis=1).astype(np.int32)

    cov = np.vstack([win_sum((m >= 1) & (m <= 5)) for m in mat])
    mm = np.vstack([win_sum((m >= 2) & (m <= 5)) for m in mat])
    ncov = ((mat >= 1) & (mat <= 5)).sum(axis=0)
    cov2 = win_sum(ncov >= PR.MIN_REFS_COVERED)
    del mat
    ss_cov, ss_mm = [], []
    for s in ss_sources:
        if s is None:
            ss_cov.append(np.zeros(nw, dtype=np.int32))
            ss_mm.append(np.zeros(nw, dtype=np.int32))
            continue
        c = extract.load_codes(workdir, s, chrom["name"], L)
        ss_cov.append(win_sum((c >= 1) & (c <= 5)))
        ss_mm.append(win_sum((c >= 2) & (c <= 5)))
    nsite = np.bincount(pos // W, minlength=nw).astype(np.int32)
    tmp = out + ".tmp.npz"
    np.savez(tmp, pos=pos, key=key, cov=cov, mm=mm, cov2=cov2, nsite=nsite,
             ss_cov=np.vstack(ss_cov), ss_mm=np.vstack(ss_mm))
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

    # ---- model
    allkeys = np.unique(np.concatenate([d["key"] for d in data.values()]))
    mchroms = []
    for c in chroms:
        d = data[c["name"]]
        mchroms.append({"name": c["name"], "length": c["length"], "label": c["label_idx"],
                        "site_pos": d["pos"].astype(np.int64),
                        "site_pid": np.searchsorted(allkeys, d["key"]).astype(np.int64)})
    nsite = sum(m["site_pos"].size for m in mchroms)
    log(f"model: {nsite:,} informative sites, {allkeys.size} site patterns")
    fit = model.run(mchroms, S, allkeys.size, verbose=lambda m: log(m, detail=True))
    res = model.resolution(fit["logp"], mchroms, S, fit["rho"])
    log(f"model: breakpoint rate {fit['rho']:.3g}/bp; "
        + "; ".join(f"{labels[r['donor']]}->{labels[r['recipient']]} min HE {output.fmt_bp(r['min_he_bp'])}"
                    + ("" if r["resolvable"] else " (unresolvable)") for r in res))

    # ---- calls
    segs, he = calls.make_calls(mchroms, chroms, labels, fit, res, data, gaps, work, ref_sources, log)
    nocall = calls.nocall_regions(mchroms, chroms, data, res, refs)
    chrom_rows = calls.chromosome_summary(chroms, labels, segs, he, nocall)

    # ---- outputs
    pattern_names = calls.pattern_names(allkeys, fit["remap"], fit["n_pat"], K)
    output.write_all(outdir, labels, refs, chroms, segs, he, nocall, chrom_rows, fit, res,
                     pattern_names, data)
    plot.painting(outdir, labels, chroms, segs, he, nocall)
    npass = sum(1 for h in he if h["filter"] == "PASS")
    log(f"done: {npass} PASS homoeologous exchanges ({len(he)} candidate segments); outputs in {outdir}")
    return 0
