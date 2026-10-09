"""Input parsing: polyploid subgenome FASTAs with labels, reference FASTAs, chromosome table."""
import gzip
import os
import re
import subprocess

from . import align
from .align import HecallError


def parse_labelled(arg, kind):
    """'LABEL=path' -> (LABEL, path); for references the label is optional (default: file stem)."""
    if "=" in arg and not os.path.exists(arg):
        lab, path = arg.split("=", 1)
    else:
        lab, path = None, arg
    if not os.path.exists(path):
        raise HecallError(f"hecall: {kind} file not found: {path}")
    if lab is None:
        if kind == "subgenome":
            raise HecallError(f"hecall: subgenome needs a label, e.g. -s A={path}")
        lab = re.sub(r"\.(fa|fasta|fna)(\.gz)?$", "", os.path.basename(path))
    return lab, os.path.abspath(path)


def fasta_lengths(path):
    """[(name, length)] from .fai if present, else by streaming the FASTA."""
    fai = path + ".fai"
    if os.path.exists(fai):
        with open(fai) as fh:
            return [(l.split("\t")[0], int(l.split("\t")[1])) for l in fh if l.strip()]
    out, name, n = [], None, 0
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    out.append((name, n))
                name, n = line[1:].split()[0], 0
            else:
                n += len(line.strip())
    if name is not None:
        out.append((name, n))
    return out


def chromosome_table(subgenomes):
    """subgenomes: [(label, path)] -> list of dict(name, length, label, label_idx, file)."""
    labels = []
    for lab, _ in subgenomes:
        if lab not in labels:
            labels.append(lab)
    chroms, seen = [], set()
    for lab, path in subgenomes:
        for name, n in fasta_lengths(path):
            if name in seen:
                raise HecallError(f"hecall: chromosome name '{name}' occurs twice among subgenome files")
            seen.add(name)
            chroms.append({"name": name, "length": n, "label": lab, "label_idx": labels.index(lab),
                           "file": path})
    if len(labels) < 2:
        raise HecallError("hecall: need at least two subgenomes (-s) to call homoeologous exchanges")
    return chroms, labels


def concat_fastas(paths, out):
    """Concatenate FASTAs into one file FastGA can read (gzip members concatenate validly)."""
    if os.path.exists(out):
        return out
    gz = [p.endswith(".gz") for p in paths]
    tmp = out + ".tmp"
    with open(tmp, "wb") as fo:
        for p in paths:
            if all(gz) or not any(gz):
                with open(p, "rb") as fi:
                    while True:
                        b = fi.read(1 << 24)
                        if not b:
                            break
                        fo.write(b)
            else:
                raise HecallError("hecall: give subgenome FASTAs either all gzipped or all plain")
    os.replace(tmp, out)
    return out


def gdb_gaps(gdb_root, tools):
    """Assembly gaps (runs of N between contigs) per scaffold, from GDBshow -h."""
    r = subprocess.run([tools["GDBshow"], "-h", gdb_root + ".1gdb"], capture_output=True, text=True,
                       env=align._env(tools))
    if r.returncode != 0:
        raise HecallError(f"hecall: GDBshow failed on {gdb_root}")
    contigs = {}
    pat = re.compile(r"[<\[](\d+),(\d+)[>\]] :: Contig")
    for line in r.stdout.splitlines():
        if not line.startswith(">"):
            continue
        name = line[1:].split()[0]
        m = pat.search(line)
        if m:
            contigs.setdefault(name, []).append((int(m.group(1)), int(m.group(2))))
    gaps = {}
    for name, cl in contigs.items():
        cl.sort()
        gaps[name] = [(a[1], b[0]) for a, b in zip(cl[:-1], cl[1:]) if b[0] > a[1]]
    return gaps
