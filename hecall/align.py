"""Whole-genome alignment with FastGA, filtered to syntenic chains with ALNchain.

FastGA (Myers) finds all local alignments >= 100 bp at >= 70% identity from adaptive seeds and is
~5x faster and ~40x lighter in memory than minimap2 -x asm20 on a 1.2 Gb vs 0.6 Gb plant pair
(docs/design.md, section 3). ALNchain keeps, per pair of sequences, the best-scoring collinear
chains (one-to-one), which removes most isolated repeat alignments. Defaults of both tools are
kept: on Coffea they reproduce minimap2's primary-chain coverage and divergence.

Only the chained .1aln is kept on disk; PAF with cs tags is streamed from it on demand.
"""
import os
import shutil
import subprocess

TOOLS = ("FastGA", "ALNchain", "ALNtoPAF", "FAtoGDB", "GIXmake", "GIXrm", "GDBshow")


def find_tools():
    """Return dict tool -> path. FastGA calls its helpers through PATH, so the directory holding
    FastGA is prepended to PATH for child processes. HECALL_FASTGA may point to that directory."""
    d = os.environ.get("HECALL_FASTGA")
    paths = {}
    for t in TOOLS:
        p = os.path.join(d, t) if d else shutil.which(t)
        if not p or not os.access(p, os.X_OK):
            raise SystemExit(f"hecall: required tool '{t}' not found (put the FASTGA directory on PATH "
                             f"or set HECALL_FASTGA)")
        paths[t] = p
    return paths


def _env(tools):
    env = dict(os.environ)
    env["PATH"] = os.path.dirname(tools["FastGA"]) + os.pathsep + env.get("PATH", "")
    return env


def run(cmd, log, tools, stdout=None):
    with open(log, "a") as lf:
        lf.write("$ " + " ".join(cmd) + "\n")
        lf.flush()
        r = subprocess.run(cmd, stdout=stdout or lf, stderr=lf, env=_env(tools))
    if r.returncode != 0:
        raise SystemExit(f"hecall: command failed ({r.returncode}): {' '.join(cmd)}; see {log}")


def make_gdb(fasta, gdb_root, tools, log):
    """FASTA(.gz) -> GDB (2-bit compressed, random access). Idempotent."""
    if not os.path.exists(gdb_root + ".1gdb"):
        run([tools["FAtoGDB"], fasta, gdb_root + ".1gdb"], log, tools)
    return gdb_root


def align_chain(query_gdb, target_gdb, out_root, threads, tmpdir, tools, log):
    """FastGA query vs target (or query vs itself when target is None) -> chained .1aln.

    Returns path of the chained .1aln. Skips work whose output exists (checkpoint)."""
    chain = out_root + ".chain.1aln"
    if os.path.exists(chain):
        return chain
    raw = out_root + ".1aln"
    os.makedirs(tmpdir, exist_ok=True)
    cmd = [tools["FastGA"], f"-T{threads}", f"-P{tmpdir}", f"-1:{raw}", query_gdb + ".1gdb"]
    if target_gdb is not None:
        cmd.append(target_gdb + ".1gdb")
    run(cmd, log, tools)
    run([tools["ALNchain"], f"-o{chain}.tmp.1aln", raw], log, tools)
    os.replace(chain + ".tmp.1aln", chain)
    os.remove(raw)
    return chain


def paf_stream(chain_1aln, threads, tools, swap=False):
    """Popen streaming PAF lines with short cs tags from a .1aln (query = first genome)."""
    cmd = [tools["ALNtoPAF"], "-s", f"-T{threads}"] + (["-w"] if swap else []) + [chain_1aln]
    return subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
                            env=_env(tools), bufsize=1 << 20)
