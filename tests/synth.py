"""Tiny synthetic allotetraploid with planted homoeologous exchanges.

ancestor -> donor A, donor B (divergence D)
donor X  -> reference X (within-species divergence d_ref), polyploid subgenome X (d_poly)
Planted on the polyploid (truth):
  chr1: non-reciprocal HE  A[100k,140k) <- B   (B chr1 unchanged: both copies B-derived)
  chr2: reciprocal HE      A[200k,260k) <-> B[200k,260k)
  chr2: small non-reciprocal HE A[50k,56k) <- B
  chr3: clean, plus a decoy: a 3-kb piece of donor-B chr3 inserted into A chr3 at 300k
        (inter-subgenome transposition; not homoeologous; must not be called)
"""
import random

BASES = "ACGT"
CHROM_LEN = 400_000
TRUTH = [  # (chrom label suffix, subgenome carrying it, start, end, donor)
    ("chr1", "A", 100_000, 140_000, "B"),
    ("chr2", "A", 200_000, 260_000, "B"),
    ("chr2", "B", 200_000, 260_000, "A"),
    ("chr2", "A", 50_000, 56_000, "B"),
]


def mutate(seq, rate, rng, indel_frac=0.1):
    """SNVs at `rate`, plus 1-3 bp indels at rate*indel_frac."""
    out = []
    i = 0
    n = len(seq)
    while i < n:
        r = rng.random()
        if r < rate:
            b = seq[i]
            out.append(rng.choice([x for x in BASES if x != b]))
            i += 1
        elif r < rate * (1 + indel_frac):
            k = rng.randint(1, 3)
            if rng.random() < 0.5:
                i += k                      # deletion
            else:
                out.append(seq[i])
                out.extend(rng.choice(BASES) for _ in range(k))  # insertion
                i += 1
        else:
            out.append(seq[i])
            i += 1
    return "".join(out)


def mutate_snv(seq, rate, rng):
    """SNVs only (keeps coordinates fixed, so planted HEs keep their coordinates)."""
    s = list(seq)
    for i in range(len(s)):
        if rng.random() < rate:
            s[i] = rng.choice([x for x in BASES if x != s[i]])
    return "".join(s)


def write_fa(path, recs):
    with open(path, "w") as fh:
        for name, seq in recs:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")


def build(outdir, D=0.06, d_ref=0.005, d_poly=0.001, d_ref_B=None, seed=1):
    """Write polyA.fa, polyB.fa, refA.fa, refB.fa into outdir; return truth list."""
    rng = random.Random(seed)
    d_ref_B = d_ref if d_ref_B is None else d_ref_B
    chroms = ["chr1", "chr2", "chr3"]
    anc = {c: "".join(rng.choice(BASES) for _ in range(CHROM_LEN)) for c in chroms}
    # donors diverge by SNVs only so homoeologous coordinates stay identical (truth is simple)
    donA = {c: mutate_snv(anc[c], D / 2, rng) for c in chroms}
    donB = {c: mutate_snv(anc[c], D / 2, rng) for c in chroms}
    refA = [(f"rA_{c}", mutate(donA[c], d_ref, rng)) for c in chroms]
    refB = [(f"rB_{c}", mutate(donB[c], d_ref_B, rng)) for c in chroms]
    pA = {c: list(donA[c]) for c in chroms}
    pB = {c: list(donB[c]) for c in chroms}
    # planted HEs (on donor sequence, before polyploid-specific mutation)
    pA["chr1"][100_000:140_000] = donB["chr1"][100_000:140_000]
    pA["chr2"][200_000:260_000] = donB["chr2"][200_000:260_000]
    pB["chr2"][200_000:260_000] = donA["chr2"][200_000:260_000]
    pA["chr2"][50_000:56_000] = donB["chr2"][50_000:56_000]
    polyA = {c: mutate_snv("".join(pA[c]), d_poly, rng) for c in chroms}
    polyB = {c: mutate_snv("".join(pB[c]), d_poly, rng) for c in chroms}
    # decoy: transposed piece of donor-B chr3 [50k,53k) inserted into A chr3 at 300k
    te = donB["chr3"][50_000:53_000]
    polyA["chr3"] = polyA["chr3"][:300_000] + te + polyA["chr3"][300_000:]
    write_fa(f"{outdir}/polyA.fa", [(f"pA_{c}", polyA[c]) for c in chroms])
    write_fa(f"{outdir}/polyB.fa", [(f"pB_{c}", polyB[c]) for c in chroms])
    write_fa(f"{outdir}/refA.fa", refA)
    write_fa(f"{outdir}/refB.fa", refB)
    return [(f"p{sub}_{c}", s, e, donor) for c, sub, s, e, donor in TRUTH]


if __name__ == "__main__":
    import sys
    print(build(sys.argv[1]))
