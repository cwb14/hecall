"""Whole-genome alignment with FastGA, filtered to syntenic chains with ALNchain.

FastGA (Myers) finds all local alignments >= 100 bp at >= 70% identity from adaptive seeds and is
~5x faster and ~40x lighter in memory than minimap2 -x asm20 on a 1.2 Gb vs 0.6 Gb plant pair
(docs/design.md, section 3). ALNchain keeps, per pair of sequences, the best-scoring collinear
chains (one-to-one), which removes most isolated repeat alignments. Defaults of both tools are
kept: on Coffea they reproduce minimap2's primary-chain coverage and divergence.

Only the chained .1aln is kept on disk; PAF with cs tags is streamed from it on demand.
"""
import os
import shutil
import subprocess

from . import params as PR


class HecallError(Exception):
    """A fatal, explained error (bad input, missing tool, failed external command)."""

TOOLS = ("FastGA", "ALNchain", "ALNtoPAF", "PAFtoALN", "FAtoGDB", "GIXmake", "GIXrm", "GDBshow")


def find_tools():
    """Return dict tool -> path. FastGA calls its helpers through PATH, so the directory holding
    FastGA is prepended to PATH for child processes. HECALL_FASTGA may point to that directory."""
    d = os.environ.get("HECALL_FASTGA")
    paths = {}
    for t in TOOLS:
        p = os.path.join(d, t) if d else shutil.which(t)
        if not p or not os.access(p, os.X_OK):
            raise HecallError(f"hecall: required tool '{t}' not found (put the FASTGA directory on PATH "
                             f"or set HECALL_FASTGA)")
        paths[t] = p
    return paths


def _env(tools):
    env = dict(os.environ)
    env["PATH"] = os.path.dirname(tools["FastGA"]) + os.pathsep + env.get("PATH", "")
    return env


def run(cmd, log, tools, stdout=None):
    with open(log, "a") as lf:
        lf.write("$ " + " ".join(cmd) + "\n")
        lf.flush()
        r = subprocess.run(cmd, stdout=stdout or lf, stderr=lf, env=_env(tools))
    if r.returncode != 0:
        raise HecallError(f"hecall: command failed ({r.returncode}): {' '.join(cmd)}; see {log}")


def make_gdb(fasta, gdb_root, tools, log):
    """FASTA(.gz) -> GDB (2-bit compressed, random access). Idempotent."""
    if not os.path.exists(gdb_root + ".1gdb"):
        run([tools["FAtoGDB"], fasta, gdb_root + ".1gdb"], log, tools)
    return gdb_root


def _prefilter_keep(raw, threads, tools):
    """Indices of raw alignments lying in collinear runs that could reach ALNchain's minimum chain
    score (a superset of what ALNchain keeps).

    Alignments are grouped by (query, target, strand, diagonal band of PREFILTER_BAND_BP; two grids
    shifted by half a band); consecutive alignments of a group within PREFILTER_GAP_BP on the query
    form a run; runs with >= PREFILTER_MIN_BP aligned query bases are kept. Isolated repeat hits
    (the bulk of raw alignments in large TE-rich genomes) fall away, which keeps ALNchain's chaining
    near-linear; on Coffea the chained result is unchanged (99.4-99.6% identical per-base codes)."""
    import numpy as np
    import pandas as pd
    p = subprocess.Popen([tools["ALNtoPAF"], f"-T{threads}", raw], stdout=subprocess.PIPE,
                         stderr=subprocess.DEVNULL, env=_env(tools))
    df = pd.read_csv(p.stdout, sep="\t", header=None, usecols=[0, 2, 3, 4, 5, 7, 8],
                     names=["q", "qs", "qe", "st", "t", "ts", "te"], engine="c",
                     dtype={"q": "category", "t": "category", "st": "category"})
    p.wait()
    if p.returncode != 0:
        raise HecallError(f"hecall: ALNtoPAF failed on {raw}")
    qs, qe = df["qs"].to_numpy(np.int64), df["qe"].to_numpy(np.int64)
    ts, te = df["ts"].to_numpy(np.int64), df["te"].to_numpy(np.int64)
    fw = (df["st"] == "+").to_numpy()
    qi, ti = df["q"].cat.codes.to_numpy(np.int64), df["t"].cat.codes.to_numpy(np.int64)
    off = np.where(fw, ts - qs, ts + qe)
    keep = np.zeros(qs.size, dtype=bool)
    D = PR.PREFILTER_BAND_BP
    for shift in (0, D // 2):
        band = np.floor_divide(off + shift, D)
        o = np.lexsort((qs, band, fw, ti, qi))
        key = np.stack([qi[o], ti[o], fw[o].astype(np.int64), band[o]])
        newkey = np.r_[True, (np.diff(key, axis=1) != 0).any(axis=0)]
        grp = np.cumsum(newkey) - 1
        cm = pd.Series(qe[o]).groupby(grp).cummax().to_numpy()
        newrun = newkey | ((qs[o] - np.r_[0, cm[:-1]]) > PR.PREFILTER_GAP_BP)
        run = np.cumsum(newrun) - 1
        score = np.bincount(run, weights=(qe[o] - qs[o]).astype(float))
        keep[o[score[run] >= PR.PREFILTER_MIN_BP]] = True
    return keep, qs


def align_chain(query_gdb, target_gdb, out_root, threads, tmpdir, tools, log):
    """FastGA query vs target (or query vs itself when target is None) -> chained .1aln.

    raw FastGA alignments -> collinear-run prefilter -> ALNchain (defaults). Returns the chained
    .1aln path. Checkpoints: the chained file; the raw file once FastGA has finished (.done)."""
    chain = out_root + ".chain.1aln"
    if os.path.exists(chain):
        return chain
    raw = out_root + ".1aln"
    os.makedirs(tmpdir, exist_ok=True)
    if not os.path.exists(raw + ".done"):
        cmd = [tools["FastGA"], f"-T{threads}", f"-P{tmpdir}", f"-1:{raw}", query_gdb + ".1gdb"]
        if target_gdb is not None:
            cmd.append(target_gdb + ".1gdb")
        run(cmd, log, tools)
        open(raw + ".done", "w").close()
    keep, qs = _prefilter_keep(raw, threads, tools)
    with open(log, "a") as lf:
        lf.write(f"# prefilter {raw}: kept {int(keep.sum())} of {keep.size} alignments\n")
    pre = out_root + ".pre"
    alnfilter = find_alnfilter()
    if alnfilter:  # copy kept records directly (binary, seconds)
        keep.astype("uint8").tofile(pre + ".mask")
        run([alnfilter, raw, pre + ".mask", pre + ".1aln"], log, tools)
    else:          # fallback: CIGAR PAF of every raw alignment -> filter -> PAFtoALN (slow on big inputs)
        _filter_via_paf(raw, keep, qs, pre, query_gdb, target_gdb, threads, tools, log)
    run([tools["ALNchain"], f"-o{chain}.tmp.1aln", pre + ".1aln"], log, tools)
    os.replace(chain + ".tmp.1aln", chain)
    for f in (raw, raw + ".done", pre + ".paf", pre + ".1aln", pre + ".mask"):
        if os.path.exists(f):
            os.remove(f)
    return chain


def find_alnfilter():
    """The alnfilter helper (csrc/build.sh -> bin/alnfilter in this repository), or on PATH."""
    here = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "bin", "alnfilter")
    if os.access(here, os.X_OK):
        return here
    return shutil.which("alnfilter")


def _filter_via_paf(raw, keep, qs, pre, query_gdb, target_gdb, threads, tools, log):
    p = subprocess.Popen([tools["ALNtoPAF"], "-x", f"-T{threads}", raw], stdout=subprocess.PIPE,
                         stderr=subprocess.DEVNULL, env=_env(tools), bufsize=1 << 20)
    n = 0
    with open(pre + ".paf", "wb") as fo:
        for line in p.stdout:
            if n >= keep.size or int(line.split(b"\t", 3)[2]) != qs[n]:
                p.kill()
                raise HecallError(f"hecall: ALNtoPAF passes over {raw} disagree at line {n + 1}")
            if keep[n]:
                fo.write(line)
            n += 1
    p.wait()
    if p.returncode != 0 or n != keep.size:
        raise HecallError(f"hecall: ALNtoPAF -x on {raw} gave {n} lines, expected {keep.size}")
    srcs = [query_gdb + ".1gdb"] + ([target_gdb + ".1gdb"] if target_gdb is not None else [])
    run([tools["PAFtoALN"], f"-T{threads}", pre + ".paf"] + srcs, log, tools)


def paf_stream(chain_1aln, threads, tools, swap=False, tmpdir=None):
    """Popen streaming PAF lines with short cs tags from a .1aln, grouped by query.

    ALNtoPAF writes records in order of the first genome. With swap (-w) the second genome becomes
    the query, so its records are interleaved; they are then sorted by query name (external sort,
    bounded memory) so that consumers can stream one query sequence at a time."""
    cmd = [tools["ALNtoPAF"], "-s", f"-T{threads}"] + (["-w"] if swap else []) + [chain_1aln]
    if not swap:
        return subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
                                env=_env(tools), bufsize=1 << 20)
    import shlex
    tmp = tmpdir or os.path.dirname(os.path.abspath(chain_1aln))
    pipe = (" ".join(shlex.quote(c) for c in cmd) + " 2>/dev/null | LC_ALL=C sort -t $'\\t' -k1,1 -s "
            f"-S 2G -T {shlex.quote(tmp)}")
    return subprocess.Popen(["bash", "-o", "pipefail", "-c", pipe], stdout=subprocess.PIPE,
                            stderr=subprocess.DEVNULL, text=True, env=_env(tools), bufsize=1 << 20)
