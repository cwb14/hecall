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


def test_alnchain_many_tied_end_positions(tmp_path):
    """Thousands of alignments ending at one target position (a repeat array cut by a contig edge)
    made upstream ALNchain's quickselect quadratic (Digitaria: > 80 min). The tie-safe build must
    chain them in seconds."""
    import random
    import subprocess
    import pytest
    try:
        tools = align.find_tools()
    except align.HecallError:
        pytest.skip("FASTGA tools not on PATH")
    if align._repo_bin("ALNchain") is None:
        pytest.skip("tie-safe ALNchain not built (csrc/build.sh)")
    rng = random.Random(1)
    unit = "".join(rng.choice("ACGT") for _ in range(200))
    n = 15000
    (tmp_path / "q.fa").write_text(">q\n" + unit * n + "\n")
    tgt = "".join(rng.choice("ACGT") for _ in range(9800)) + unit
    (tmp_path / "t.fa").write_text(">t\n" + tgt + "\n")
    paf = "".join(f"q\t{200 * n}\t{200 * i}\t{200 * i + 200}\t+\tt\t10000\t9800\t10000\t200\t200\t60"
                  f"\tcg:Z:200=\n" for i in range(n))
    (tmp_path / "x.paf").write_text(paf)
    env = align._env(tools)
    for f in ("q", "t"):
        subprocess.run([tools["FAtoGDB"], str(tmp_path / f"{f}.fa")], check=True, env=env)
    subprocess.run([tools["PAFtoALN"], "-T1", str(tmp_path / "x.paf"), str(tmp_path / "q.1gdb"),
                    str(tmp_path / "t.1gdb")], check=True, env=env, cwd=tmp_path)
    r = subprocess.run([tools["ALNchain"], f"-o{tmp_path / 'c.1aln'}", str(tmp_path / "x.1aln")],
                       env=env, timeout=60, capture_output=True)
    assert r.returncode == 0
