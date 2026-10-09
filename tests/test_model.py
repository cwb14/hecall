"""Window-level origin HMM on simulated windows (no aligner needed)."""
import numpy as np

from hecall import model


def sim(label, n, rng, foreign=(), mu=(1.2, -1.2), sd=0.6, span=2000):
    """1-D origin feature per window: state 0 ~ N(mu0, sd), state 1 ~ N(mu1, sd)."""
    st = np.full(n, label)
    for a, b, s in foreign:
        st[a:b] = s
    x = rng.normal(np.take(mu, st), sd)
    start = np.arange(n) * span
    return {"label": label, "X": x[:, None], "d": np.full(n, float(span)), "start": start,
            "end": start + span, "length": n * span}


def test_recovers_planted_run_and_null_chromosomes():
    rng = np.random.default_rng(1)
    wins = [sim(0, 1000, rng, foreign=[(400, 420, 1)]), sim(0, 1000, rng), sim(1, 1000, rng), sim(1, 1000, rng)]
    fit = model.run(wins, 2)
    found = []
    for i, (w, p) in enumerate(zip(wins, fit["post"])):
        st = p.argmax(axis=1)
        for a, b, s in model._runs(st):
            if s != w["label"]:
                found.append((i, a, b, model.qual(p[a:b, s].mean())))
    assert len(found) == 1
    i, a, b, q = found[0]
    assert i == 0 and abs(a - 400) <= 2 and abs(b - 420) <= 2 and q >= 20
    res = model.resolution(fit, wins, 2, fit["rho"], fit["rho_ac"])
    assert all(r["resolvable"] for r in res)


def test_indistinguishable_origins_unresolvable():
    rng = np.random.default_rng(2)
    wins = [sim(0, 500, rng, mu=(0.0, 0.05)), sim(1, 500, rng, mu=(0.0, 0.05))]
    fit = model.run(wins, 2)
    assert not any(r["resolvable"] for r in model.resolution(fit, wins, 2, fit["rho"], fit["rho_ac"]))


def test_refine_breakpoint_finds_switch():
    pos = np.arange(0, 1000, 10)
    left = np.where(pos < 500, -0.1, -2.0)
    right = np.where(pos < 500, -2.0, -0.1)
    _, ml_lo, ml_hi, ci_lo, ci_hi = model.refine_breakpoint(pos, left, right)
    assert ml_lo <= 500 <= ml_hi and ci_lo <= ml_lo and ci_hi >= ml_hi


def test_segment_cuts_ordered_and_positive():
    """Short foreign runs next to each other: tiles must stay ordered with positive length."""
    rng = np.random.default_rng(3)
    w = sim(0, 300, rng, foreign=[(100, 102, 1), (104, 106, 1)])
    w["start"], w["end"] = np.arange(300) * 2000, np.arange(300) * 2000 + 2000
    pos = np.arange(0, 600_000, 50, dtype=np.int64)
    chrom = {"length": 600_000, "label": 0, "site_pos": pos, "site_pid": rng.integers(0, 2, pos.size)}
    post = np.zeros((300, 2))
    post[:, 0] = 1
    post[[100, 101, 104, 105]] = [0.2, 0.8]
    binned = {"E": np.zeros((300, 2))}
    slogp = np.log(np.array([[0.8, 0.2], [0.2, 0.8]]))
    segs = model.segments(chrom, w, post, binned, slogp)
    tiles = [(s["tile_start"], s["tile_end"]) for s in segs]
    assert all(b > a for a, b in tiles)
    assert all(t1[0] == t0[1] for t0, t1 in zip(tiles[:-1], tiles[1:]))


def test_tau_n_grows_with_segment_length_under_correlation():
    rho = np.zeros(200)
    rho[0] = 1
    rho[1:50] = 0.5 * np.exp(-np.arange(1, 50) / 10)
    assert model.tau_n(rho, 1) == 1.0
    assert 1.0 < model.tau_n(rho, 5) < model.tau_n(rho, 100)
    assert model.tau_n(np.r_[1.0, np.zeros(10)], 50) == 1.0


def test_state_density_learned_from_native_windows_only():
    """Windows assigned to state 1 on a state-0 chromosome (candidate exchanges) must not shape
    state 1's density; otherwise ambiguous windows get absorbed round after round (v0.4.0)."""
    rng = np.random.default_rng(4)
    w0 = sim(0, 2000, rng, foreign=[(500, 700, 1)])
    w0["X"][500:700, 0] = rng.normal(-0.4, 0.05, 200)   # ambiguous block, nearer state 1's tail
    w1 = sim(1, 2000, rng)
    fit = model.run([w0, w1], 2)
    native1 = w1["X"][fit["assign"][1] == 1]
    ref = model.KDE(native1, fit["lo"], fit["hi"])
    grid = np.linspace(-3, 3, 61)[:, None]
    assert np.allclose(fit["kdes"][1].logpdf(grid), ref.logpdf(grid))
