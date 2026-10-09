import numpy as np

from hecall import calls


def aln(rows):
    """rows: (qs, qe, fwd, target, ts, te)"""
    return {"qs": np.array([r[0] for r in rows]), "qe": np.array([r[1] for r in rows]),
            "fwd": np.array([r[2] for r in rows]), "t": np.array([r[3] for r in rows], dtype=object),
            "ts": np.array([r[4] for r in rows]), "te": np.array([r[5] for r in rows]),
            "div": np.zeros(len(rows))}


def test_projection_follows_collinear_chain_not_scattered_repeats():
    # a 9-Mb segment aligned collinearly to Y:51-60 Mb (offset +6 Mb), plus repeat hits all over Y
    rows = [(q, q + 400_000, True, "Y", q + 6_000_000, q + 6_400_000) for q in range(45_000_000, 54_000_000, 500_000)]
    rows += [(45_100_000 + i * 700_000, 45_150_000 + i * 700_000, True, "Y", 2_000_000 + i * 3_000_000,
              2_050_000 + i * 3_000_000) for i in range(12)]
    p = calls.project(aln(rows), 45_000_000, 54_000_000)
    assert p["target"] == "Y"
    assert 50_500_000 <= p["lo"] and p["hi"] <= 60_500_000


def test_synteny_classes_homoeologous_translocation_transposed_repeat():
    chain = [(q, q + 90_000, True, "Y", q, q + 90_000) for q in range(0, 3_000_000, 100_000)]
    cls, _ = calls.syntenic(aln(chain), 1_000_000, 1_200_000, 1_000_000, 1_200_000)
    assert cls == "homoeologous"
    # a 200-kb block whose only match is a collinear block on another chromosome: translocation
    b = aln([r for r in chain if not (1_000_000 <= r[0] < 1_200_000)] +
            [(1_000_000 + i * 20_000, 1_018_000 + i * 20_000, True, "Z", 40_000_000 + i * 20_000,
              40_018_000 + i * 20_000) for i in range(10)])
    cls, tgt = calls.syntenic(b, 1_000_000, 1_200_000, 1_000_000, 1_200_000)
    assert cls == "translocation" and tgt.startswith("Z:")
    # a 6-kb segment matching a repeat copy elsewhere (transposed TE): nonsyntenic
    c = aln([r for r in chain if not (1_000_000 <= r[0] < 1_100_000)] +
            [(1_000_000, 1_006_000, True, "Z", 40_000_000, 40_006_000)])
    cls, _ = calls.syntenic(c, 1_000_000, 1_006_000, 1_000_000, 1_006_000)
    assert cls == "nonsyntenic"


def test_label_coherence_and_filter():
    """A label painted mostly with other origins is not a subgenome: its calls are filtered."""
    from hecall import calls
    chroms = [{"name": "c1", "label_idx": 0}, {"name": "c2", "label_idx": 0}, {"name": "c3", "label_idx": 1}]
    seg = lambda a, b, s: {"tile_start": a, "tile_end": b, "state": s}
    segs = {"c1": [seg(0, 100, 1), seg(100, 130, 0)],          # label 0 chromosomes mostly state 1
            "c2": [seg(0, 100, 2), {"tile_start": 100, "tile_end": 900, "state": 0, "nosites": True}],
            "c3": [seg(0, 95, 1), seg(95, 100, 0)]}
    coh = calls.label_coherence(chroms, segs, 3)
    assert abs(coh[0] - 30 / 230) < 1e-9 and abs(coh[1] - 0.95) < 1e-9 and coh[2] != coh[2]   # nan
    h = {"qual": 99, "label_idx": 1, "donor_idx": 0, "syntenic": "homoeologous"}
    res = {(1, 0): True}
    assert calls._filters(h, res, incoherent={0}) == "incoherent_label"
    assert calls._filters(h, res, incoherent=set()) == "PASS"
