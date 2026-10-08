"""hecall command line."""
import argparse

from . import __version__


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="hecall",
        description="Call homoeologous exchanges in a polyploid assembly from its subgenome FASTAs and "
                    "reference genomes of the donor lineages (or stand-ins).")
    ap.add_argument("-s", "--subgenome", action="append", required=True, metavar="LABEL=FASTA",
                    help="polyploid subgenome FASTA with its label; give once per subgenome (>= 2)")
    ap.add_argument("-r", "--ref", action="append", required=True, metavar="[NAME=]FASTA",
                    help="reference genome FASTA (diploid progenitor, stand-in, or one subgenome of a "
                         "polyploid parent); give once per reference (>= 2)")
    ap.add_argument("-o", "--outdir", required=True, help="output directory (resumes if it exists)")
    ap.add_argument("-t", "--threads", type=int, default=8, help="threads [8]")
    ap.add_argument("-v", "--verbose", action="store_true", help="per-step progress")
    ap.add_argument("--version", action="version", version=f"hecall {__version__}")
    a = ap.parse_args(argv)
    from . import pipeline
    return pipeline.run(a.subgenome, a.ref, a.outdir, threads=a.threads, verbose=a.verbose)
