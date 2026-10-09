"""End-to-end on a tiny synthetic allotetraploid (needs FastGA on PATH or HECALL_FASTGA)."""
import csv
import os
import shutil

import pytest

import synth
from hecall import pipeline

pytestmark = pytest.mark.skipif(shutil.which("FastGA") is None and not os.environ.get("HECALL_FASTGA"),
                                reason="FastGA not available")


def _calls(out):
    with open(os.path.join(out, "he.tsv")) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _match(calls, truth, tol=500):
    hit = []
    for chrom, a, b, donor in truth:
        m = [c for c in calls if c["chrom"] == chrom and c["donor"] == donor and c["filter"] == "PASS"
             and abs(int(c["start"]) - a) <= tol and abs(int(c["end"]) - b) <= tol]
        hit.append(len(m) == 1)
    return hit


def test_planted_exchanges_types_and_decoy(tmp_path):
    truth = synth.build(str(tmp_path))
    out = str(tmp_path / "out")
    d = str(tmp_path)
    pipeline.run([f"A={d}/polyA.fa", f"B={d}/polyB.fa"], [f"refA={d}/refA.fa", f"refB={d}/refB.fa"], out, threads=2)
    calls = _calls(out)
    assert all(_match(calls, truth))
    passed = [c for c in calls if c["filter"] == "PASS"]
    assert len(passed) == len(truth)                       # decoy and clean chromosome: no calls
    types = {(c["chrom"], int(c["start"]) // 10000): c["type"] for c in passed}
    assert types[("pA_chr1", 9)] == "non-reciprocal"
    flags = {(c["chrom"], int(c["start"]) // 10000): c["flags"] for c in passed}
    assert "homoeolog_not_closer" not in flags[("pA_chr1", 9)]   # a true copy of the donor homoeologue
    assert types[("pA_chr2", 19)] == "reciprocal"
    assert types[("pB_chr2", 20)] == "reciprocal"
    for f in ("painting.pdf", "painting.png", "segments.bed", "nocall.bed", "chromosomes.tsv", "resolution.tsv"):
        assert os.path.exists(os.path.join(out, f))
    res = list(csv.DictReader((l for l in open(os.path.join(out, "resolution.tsv")) if not l.startswith("#")),
                              delimiter="\t"))
    assert all(0.03 < float(r["ss_div_baseline"]) < 0.15 for r in res)   # planted donor divergence 6%


def test_stand_in_reference(tmp_path):
    """Donor B represented only by a distant stand-in (3% from the donor): calls still found."""
    truth = synth.build(str(tmp_path), d_ref_B=0.03, seed=2)
    out = str(tmp_path / "out")
    d = str(tmp_path)
    pipeline.run([f"A={d}/polyA.fa", f"B={d}/polyB.fa"], [f"refA={d}/refA.fa", f"refB={d}/refB.fa"], out, threads=2)
    calls = _calls(out)
    assert all(_match(calls, truth, tol=1500))
    assert sum(c["filter"] == "PASS" for c in calls) == len(truth)


def test_missing_homoeologue_is_reported(tmp_path):
    """A subgenome chromosome without its homoeologue (LOH, aneuploidy) is stated in chromosomes.tsv."""
    synth.build(str(tmp_path))
    d = str(tmp_path)
    recs = open(f"{d}/polyB.fa").read().split(">")[1:]
    with open(f"{d}/polyB_lost.fa", "w") as fh:
        fh.write("".join(">" + r for r in recs if not r.startswith("pB_chr3")))
    out = str(tmp_path / "out")
    pipeline.run([f"A={d}/polyA.fa", f"B={d}/polyB_lost.fa"], [f"refA={d}/refA.fa", f"refB={d}/refB.fa"], out,
                 threads=2)
    rows = {r["chrom"]: r for r in csv.DictReader(open(os.path.join(out, "chromosomes.tsv")), delimiter="\t")}
    assert float(rows["pA_chr3"]["pct_homoeolog_B"]) < 5
    assert float(rows["pA_chr1"]["pct_homoeolog_B"]) > 80
    assert rows["pA_chr1"]["pct_homoeolog_A"] == "."

