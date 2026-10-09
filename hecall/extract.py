"""Stage 2: chained alignments -> per-chromosome base codes and alignment tables.

For every reference (and for the polyploid's own other subgenomes, from the self alignment) one
worker streams PAF from the chained .1aln and writes, per polyploid chromosome:
    work/codes/<source>/<chrom>.npy   uint8 codes (see sites.py)
and one alignment table per source (no cs): work/codes/<source>.aln.tsv.gz
    qname qs qe strand tname ts te div
"""
import gzip
import os

import numpy as np

from . import align, sites
from .align import HecallError


def _flush(outdir, qname, qlen, recs):
    codes = sites.codes_for_chrom(recs, qlen)
    tmp = os.path.join(outdir, f".{qname}.npy")
    np.save(tmp, codes)
    os.replace(tmp, os.path.join(outdir, f"{qname}.npy"))


def extract_source(chain, source, chrom_len, workdir, tools, threads, target_set=None, swap=False, tag=None):
    """Parse one chained alignment into per-chromosome codes.

    chrom_len: polyploid chromosomes to extract (others in the stream are ignored).
    target_set: if given, only alignments whose target is in this set are used (legacy self
    alignment split by subgenome). swap: read the alignment with genomes swapped (ALNtoPAF -w), i.e.
    the second genome becomes the query. tag: name of this extraction's done marker and alignment
    table (default: source); several extractions may fill one source (pairwise subgenome alignments).
    Writes <workdir>/codes/<source>/<chrom>.npy and skips if <tag>.done exists.
    """
    tag = tag or source
    outdir = os.path.join(workdir, "codes", source)
    done = os.path.join(workdir, "codes", f"{tag}.done")
    if os.path.exists(done):
        return
    os.makedirs(outdir, exist_ok=True)
    alnf = gzip.open(os.path.join(workdir, "codes", f"{tag}.aln.tsv.gz"), "wt", compresslevel=3)
    cur, recs = None, []
    seen = set()
    p = align.paf_stream(chain, threads, tools, swap=swap)
    for line in p.stdout:
        f = line.rstrip("\n").split("\t")
        q, t = f[0], f[5]
        if q not in chrom_len or (target_set is not None and t not in target_set) or q == t:
            continue
        cs = next((x[5:] for x in f[12:] if x.startswith("cs:Z:")), None)
        if cs is None:
            continue
        qs, qe, ts, te = int(f[2]), int(f[3]), int(f[7]), int(f[8])
        div = 1.0 - int(f[9]) / max(int(f[10]), 1)
        alnf.write(f"{q}\t{qs}\t{qe}\t{f[4]}\t{t}\t{ts}\t{te}\t{div:.4f}\n")
        if q != cur:
            if cur is not None:
                if cur in seen:
                    raise HecallError(f"hecall: PAF from {chain} is not grouped by query ({cur})")
                _flush(outdir, cur, chrom_len[cur], recs)
                seen.add(cur)
            cur, recs = q, []
        recs.append((qs, qe, f[4], cs, div))
    if cur is not None:
        _flush(outdir, cur, chrom_len[cur], recs)
    p.wait()
    if p.returncode != 0:
        raise HecallError(f"hecall: ALNtoPAF failed on {chain}")
    alnf.close()
    open(done, "w").close()


def load_codes(workdir, source, chrom, length):
    f = os.path.join(workdir, "codes", source, f"{chrom}.npy")
    return np.load(f) if os.path.exists(f) else np.zeros(length, dtype=np.uint8)


def load_aln(workdir, source):
    """Alignment tables of one source (all its extractions) as dict chrom -> arrays."""
    import glob
    import pandas as pd
    files = sorted(set(glob.glob(os.path.join(workdir, "codes", f"{source}.aln.tsv.gz")) +
                       glob.glob(os.path.join(workdir, "codes", f"{source}.from_*.aln.tsv.gz"))))
    frames = [pd.read_csv(f, sep="\t", header=None, names=["q", "qs", "qe", "st", "t", "ts", "te", "div"],
                          dtype={"q": str, "t": str, "st": str}) for f in files]
    frames = [f for f in frames if len(f)]
    if not frames:
        return {}
    df = pd.concat(frames, ignore_index=True)
    out = {}
    for q, g in df.groupby("q", sort=False):
        out[q] = {"qs": g["qs"].to_numpy(), "qe": g["qe"].to_numpy(), "fwd": (g["st"] == "+").to_numpy(),
                  "t": g["t"].to_numpy(dtype=object), "ts": g["ts"].to_numpy(), "te": g["te"].to_numpy(),
                  "div": g["div"].to_numpy()}
    return out
