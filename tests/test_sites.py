import numpy as np

from hecall import sites


def test_cs_forward_codes():
    codes = np.zeros(12, dtype=np.uint8)
    # query bases 0-2 match, base 3 is 'a' where ref has 'g', 1-bp insertion in query, 2 matches,
    # a 2-bp deletion from query (no advance), 1 match
    sites.fill_codes(codes, 0, 8, "+", ":3*ga+t:2-cc:1")
    assert list(codes[:8]) == [1, 1, 1, sites._BASE["g"], 0, 1, 1, 1]


def test_cs_reverse_codes_use_complement():
    codes = np.zeros(4, dtype=np.uint8)
    # minus strand: cs walks the target forward, the query backward from qe
    sites.fill_codes(codes, 0, 4, "-", ":1*ga:2")
    # second cs base -> query position 2; ref allele 'g' on target strand = 'c' on query strand
    assert list(codes) == [1, 1, sites._BASE["c"], 1]


def test_best_alignment_wins_overlap():
    recs = [(0, 10, "+", ":10", 0.05), (5, 15, "+", ":2*ac:7", 0.01)]
    c = sites.codes_for_chrom(recs, 15)
    assert c[7] == sites._BASE["a"]          # from the less divergent alignment
    assert (c[:5] == 1).all() and (c[10:15] == 1).all()


def test_informative_patterns_canonical():
    # 3 refs, 5 positions
    M, A, C, G = 1, 2, 3, 4
    mat = np.array([[M, A, A, 0, M],
                    [A, M, C, M, M],
                    [A, M, A, A, M]], dtype=np.uint8)
    pos, key = sites.informative_sites(mat, 2)
    # col 2: polyploid matches no reference (private) -> not diagnostic; col 4: all match
    assert list(pos) == [0, 1, 3]
    pats = [sites.key_to_pattern(int(k), 3) for k in key]
    assert pats == ["Maa", "aMM", ".Ma"]
