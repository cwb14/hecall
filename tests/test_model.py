"""HMM on simulated site patterns (no aligner needed)."""
import numpy as np

from hecall import model

# pattern ids: 0 = matches ref1 only, 1 = matches ref2 only, 2 = private, 3 = triallelic
PROFILE = {0: [0.90, 0.07, 0.02, 0.01],   # subgenome 0 (donor 1)
           1: [0.07, 0.90, 0.02, 0.01]}   # subgenome 1 (donor 2)


def sim_chrom(name, label, length, density, rng, foreign=()):
    n = rng.poisson(length * density)
    pos = np.sort(rng.choice(length, size=n, replace=False))
    state = np.full(n, label)
    for a, b, s in foreign:
        state[(pos >= a) & (pos < b)] = s
    pid = np.array([rng.choice(4, p=PROFILE[s]) for s in state])
    return {"name": name, "length": length, "label": label, "site_pos": pos.astype(np.int64),
            "site_pid": pid.astype(np.int64)}


def test_recovers_planted_segment_and_null_chromosome():
    rng = np.random.default_rng(3)
    chroms = [sim_chrom("A1", 0, 2_000_000, 0.01, rng, foreign=[(800_000, 830_000, 1)]),
              sim_chrom("A2", 0, 2_000_000, 0.01, rng),
              sim_chrom("B1", 1, 2_000_000, 0.01, rng),
              sim_chrom("B2", 1, 2_000_000, 0.01, rng)]
    fit = model.run(chroms, 2, 4)
    found = []
    for c, b, p in zip(chroms, fit["binned"], fit["post"]):
        for s in model.segments(c, b, p, fit["logp"]):
            if s["state"] != c["label"]:
                found.append((c["name"], s["start"], s["end"], model.qual(s["post_mean"])))
    assert len(found) == 1
    name, a, b, q = found[0]
    assert name == "A1" and abs(a - 800_000) < 2_000 and abs(b - 830_000) < 2_000 and q >= 20


def test_indistinguishable_profiles_unresolvable():
    rng = np.random.default_rng(4)
    global PROFILE
    saved = dict(PROFILE)
    PROFILE = {0: [0.5, 0.48, 0.01, 0.01], 1: [0.49, 0.49, 0.01, 0.01]}
    try:
        chroms = [sim_chrom("A1", 0, 1_000_000, 0.002, rng), sim_chrom("B1", 1, 1_000_000, 0.002, rng)]
        fit = model.run(chroms, 2, 4)
        res = model.resolution(fit["logp"], chroms, 2, fit["rho"])
        assert not any(r["resolvable"] for r in res)
    finally:
        PROFILE = saved


def test_refine_breakpoint_finds_switch():
    pos = np.arange(0, 1000, 10)
    left = np.where(pos < 500, -0.1, -2.0)     # log P under left state
    right = np.where(pos < 500, -2.0, -0.1)
    ml_lo, ml_hi, ci_lo, ci_hi = model.refine_breakpoint(pos, left, right)
    assert ml_lo <= 500 <= ml_hi and ci_lo <= ml_lo and ci_hi >= ml_hi
