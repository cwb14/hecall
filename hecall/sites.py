"""Per-base comparison of the polyploid (PAF query) with each reference (PAF target).

For one polyploid chromosome and one reference, `codes_from_paf` returns a uint8 array over the
chromosome:
    0 uncovered, 1 match, 2..5 mismatch (reference allele A,C,G,T on the polyploid's forward
    strand), 6 ambiguous (covered by more than one alignment of this reference).
Only SNVs are recorded; query bases inside insertions/deletions stay uncovered (0).
"""
import re
import gzip
import numpy as np

UNCOV, MATCH, AMBIG = 0, 1, 6
_BASE = {"a": 2, "c": 3, "g": 4, "t": 5}
_COMP = {"a": "t", "c": "g", "g": "c", "t": "a"}
_CS_OP = re.compile(r":(\d+)|\*([a-z])([a-z])|\+([a-z]+)|-([a-z]+)")
SYMBOLS = ".Mabcd"  # pattern alphabet: uncovered, match, distinct mismatching alleles


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def paf_records(path):
    """Yield (qname, qlen, qs, qe, strand, tname, ts, te, cs, div) for PAF lines with a cs tag.

    div = 1 - matching bases / alignment block length (PAF columns 10, 11)."""
    with open_text(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            cs = None
            for tag in f[12:]:
                if tag.startswith("cs:Z:"):
                    cs = tag[5:]
                    break
            if cs is None:
                continue
            div = 1.0 - int(f[9]) / max(int(f[10]), 1)
            yield (f[0], int(f[1]), int(f[2]), int(f[3]), f[4], f[5], int(f[7]), int(f[8]), cs, div)


def fill_codes(codes, qs, qe, strand, cs):
    """Write one alignment's per-base codes into `codes` (query-forward coordinates)."""
    fwd = strand == "+"
    q = qs if fwd else qe  # '+': next base is q; '-': next base is q-1, walking left
    mm_pos, mm_code = [], []
    for m in _CS_OP.finditer(cs):
        n, x, y, ins, dele = m.groups()
        if n is not None:
            n = int(n)
            if fwd:
                codes[q:q + n] = MATCH
                q += n
            else:
                codes[q - n:q] = MATCH
                q -= n
        elif x is not None:
            p = q if fwd else q - 1
            if x in _BASE and y in _BASE:
                mm_pos.append(p)
                mm_code.append(_BASE[x if fwd else _COMP[x]])
            q = q + 1 if fwd else q - 1
        elif ins is not None:  # bases present in query only: left uncovered
            q = q + len(ins) if fwd else q - len(ins)
        # deletion (bases in target only): query does not advance
    if mm_pos:
        codes[np.asarray(mm_pos, dtype=np.int64)] = np.asarray(mm_code, dtype=np.uint8)


def codes_for_chrom(records, qlen):
    """Combine all alignments of one reference on one polyploid chromosome.

    records: (qs, qe, strand, cs, div). Where alignments overlap on the polyploid, each base takes
    its code from the least divergent alignment: the orthologue (or homoeologue) is the closest
    copy, older paralogues and repeat copies are more diverged.
    """
    codes = np.zeros(qlen, dtype=np.uint8)
    claimed = np.zeros(qlen, dtype=bool)  # spans of better alignments, incl. their indel bases
    for (qs, qe, strand, cs, div) in sorted(records, key=lambda r: r[4]):
        free = ~claimed[qs:qe]
        if not free.any():
            continue
        tmp = np.zeros(qe - qs, dtype=np.uint8)
        fill_codes(tmp, 0, qe - qs, strand, cs)
        codes[qs:qe][free] = tmp[free]
        claimed[qs:qe] = True
    return codes


def informative_sites(code_mat, min_covered=2):
    """code_mat: (K, L) uint8 codes of K references on one chromosome.

    Returns (positions, pattern_keys). A site is informative when >= min_covered references cover
    it and at least one covered reference mismatches. The pattern key encodes one symbol per
    reference (index into SYMBOLS) in base 6, with mismatching alleles relabelled by first
    appearance so that 'aM' and 'bM' are the same pattern.
    """
    K, L = code_mat.shape
    covered = (code_mat >= 1) & (code_mat <= 5)
    n_cov = covered.sum(axis=0)
    any_mm = ((code_mat >= 2) & (code_mat <= 5)).any(axis=0)
    pos = np.flatnonzero((n_cov >= min_covered) & any_mm)
    sub = code_mat[:, pos].astype(np.int16)
    sub[(sub == AMBIG)] = 0
    n = pos.size
    label = np.full((n, 4), -1, dtype=np.int16)  # allele (A,C,G,T) -> canonical letter index
    nxt = np.zeros(n, dtype=np.int16)
    key = np.zeros(n, dtype=np.int64)
    idx = np.arange(n)
    for j in range(K):
        c = sub[j]
        sym = np.zeros(n, dtype=np.int64)          # '.'
        sym[c == MATCH] = 1                       # 'M'
        mm = c >= 2
        if mm.any():
            a = c[mm] - 2
            ii = idx[mm]
            lab = label[ii, a]
            new = lab < 0
            # assign new letters in order of reference index (first appearance)
            lab_new = nxt[ii[new]]
            label[ii[new], a[new]] = lab_new
            nxt[ii[new]] += 1
            lab = label[ii, a]
            sym[mm] = 2 + lab
        key += sym * (6 ** j)
    return pos, key


def key_to_pattern(key, K):
    s = []
    for _ in range(K):
        s.append(SYMBOLS[key % 6])
        key //= 6
    return "".join(s)
