import os
import stat

from hecall import align


def test_swapped_stream_is_grouped_by_query(tmp_path):
    """ALNtoPAF -w keeps first-genome order, so swapped queries arrive interleaved; paf_stream must
    regroup them (unsorted queries made the extractor fail and hung the worker pool, v0.3.3-0.3.5)."""
    paf = tmp_path / "x.paf"
    rows = ["B2", "B1", "B2", "B3", "B1", "B2"]
    paf.write_text("".join(f"{q}\t100\t0\t10\t+\tA1\t100\t0\t10\t10\t10\t60\tcs:Z::10\n" for q in rows))
    fake = tmp_path / "ALNtoPAF"
    fake.write_text(f"#!/bin/bash\ncat {paf}\n")
    fake.chmod(fake.stat().st_mode | stat.S_IEXEC)
    tools = {t: str(fake) for t in align.TOOLS}
    p = align.paf_stream("chain.1aln", 1, tools, swap=True, tmpdir=str(tmp_path))
    q = [l.split("\t", 1)[0] for l in p.stdout]
    p.wait()
    assert p.returncode == 0
    assert q == sorted(q, key=lambda x: x) and len(q) == len(rows)


def test_errors_are_exceptions_not_systemexit():
    assert issubclass(align.HecallError, Exception)
