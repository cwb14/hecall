"""Output tables: he.tsv, he.gff3, segments.bed, nocall.bed, chromosomes.tsv, profiles.tsv,
resolution.tsv, bins.tsv.gz."""
import gzip
import math
import os

import numpy as np

from . import params as PR


def fmt_bp(x):
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "inf"
    x = float(x)
    for unit, f in (("Mb", 1e6), ("kb", 1e3)):
        if x >= f:
            return f"{x / f:.1f} {unit}"
    return f"{x:.0f} bp"


def _f(x, nd=4):
    if x is None:
        return "."
    if isinstance(x, bool):
        return "yes" if x else "no"
    if isinstance(x, float):
        return "." if not math.isfinite(x) else f"{x:.{nd}g}"
    return str(x)


def write_all(outdir, labels, refs, chroms, segs, he, nocall, chrom_rows, fit, res, pattern_names, data):
    rnames = [r[0] for r in refs]
    # he.tsv
    cols = ["chrom", "start", "end", "length", "label", "donor", "qual", "filter", "type", "llr",
            "n_sites", "n_support", "n_against", "core_start", "core_end", "bp_left_ci", "bp_right_ci",
            "syntenic", "closest_ref", "ref_target"] + [f"div_{n}" for n in rnames] + \
           ["homoeolog", "homoeolog_frac_recipient_origin", "ss_div", "ss_div_baseline", "ss_cover", "flags"]
    with open(os.path.join(outdir, "he.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for h in sorted(he, key=lambda h: (h["chrom"], h["start"])):
            row = dict(h)
            row["closest_ref"] = rnames[h["closest_ref"]]
            row["syntenic"] = "." if h["syntenic"] is None else ("yes" if h["syntenic"] else "no")
            for n, v in zip(rnames, h["ref_div"]):
                row[f"div_{n}"] = v
            fh.write("\t".join(_f(row.get(c)) if not isinstance(row.get(c), str) else row[c] for c in cols) + "\n")
    # he.gff3 (all candidates; FILTER in attributes)
    with open(os.path.join(outdir, "he.gff3"), "w") as fh:
        fh.write("##gff-version 3\n")
        for i, h in enumerate(sorted(he, key=lambda h: (h["chrom"], h["start"])), 1):
            attr = (f"ID=HE{i};donor={h['donor']};recipient={h['label']};qual={h['qual']};"
                    f"filter={h['filter']};type={h['type']};n_sites={h['n_sites']};flags={h['flags']}")
            fh.write(f"{h['chrom']}\thecall\thomoeologous_exchange\t{h['start'] + 1}\t{h['end']}\t"
                     f"{h['qual']}\t.\t.\t{attr}\n")
    # segments.bed: full origin painting
    with open(os.path.join(outdir, "segments.bed"), "w") as fh:
        fh.write("#chrom\tstart\tend\torigin\tscore\tstrand\tlabel\tqual\n")
        for c in chroms:
            for s in segs[c["name"]]:
                origin = "." if s.get("nosites") else labels[s["state"]]
                fh.write(f"{c['name']}\t{s['tile_start']}\t{s['tile_end']}\t{origin}\t"
                         f"{min(1000, int(s['qual'] * 10))}\t.\t{c['label']}\t{s['qual']:.1f}\n")
    with open(os.path.join(outdir, "nocall.bed"), "w") as fh:
        fh.write("#chrom\tstart\tend\treason\n")
        for c in chroms:
            for a, b, r in nocall[c["name"]]:
                fh.write(f"{c['name']}\t{a}\t{b}\t{r}\n")
    with open(os.path.join(outdir, "chromosomes.tsv"), "w") as fh:
        keys = list(chrom_rows[0].keys())
        fh.write("\t".join(keys) + "\n")
        for r in chrom_rows:
            fh.write("\t".join(_f(r[k]) for k in keys) + "\n")
    # profiles
    p = np.exp(fit["logp"])
    with open(os.path.join(outdir, "profiles.tsv"), "w") as fh:
        fh.write("# site pattern: one symbol per reference (" + ",".join(rnames) +
                 "); M=matches polyploid, .=uncovered, a/b/c=distinct other alleles\n")
        fh.write("pattern\t" + "\t".join(f"freq_{l}" for l in labels) + "\t" +
                 "\t".join(f"count_{l}" for l in labels) + "\n")
        order = np.argsort(-fit["counts"].sum(axis=0))
        for j in order:
            fh.write(pattern_names[j] + "\t" + "\t".join(f"{p[k, j]:.5f}" for k in range(len(labels))) + "\t" +
                     "\t".join(f"{int(fit['counts'][k, j])}" for k in range(len(labels))) + "\n")
    with open(os.path.join(outdir, "resolution.tsv"), "w") as fh:
        fh.write("recipient\tdonor\tkl_per_site\tsites_per_kb\tllr_per_kb\tllr_needed\tmin_he_bp\t"
                 "median_chrom_bp\tresolvable\n")
        for r in res:
            fh.write(f"{labels[r['recipient']]}\t{labels[r['donor']]}\t{r['kl_per_site']:.4f}\t"
                     f"{r['sites_per_kb']:.3f}\t{r['llr_per_kb']:.4f}\t{r['llr_needed']:.2f}\t"
                     f"{_f(r['min_he_bp'], 6)}\t{r['median_chrom_bp']:.0f}\t{'yes' if r['resolvable'] else 'no'}\n")
        fh.write(f"# breakpoint rate (per bp, EM): {fit['rho']:.4g}\n")
    # bins.tsv.gz: 10-kb windows
    W = PR.WINDOW_BP
    with gzip.open(os.path.join(outdir, "bins.tsv.gz"), "wt", compresslevel=3) as fh:
        fh.write("chrom\tstart\tend\tn_sites\tcov_2refs\t" +
                 "\t".join(f"cov_{n}\tident_{n}" for n in rnames) + "\t" +
                 "\t".join(f"ss_cov_{l}\tss_ident_{l}" for l in labels) + "\t" +
                 "\t".join(f"post_{l}" for l in labels) + "\n")
        for c, binned, post in zip(chroms, fit["binned"], fit["post"]):
            d = data[c["name"]]
            nw = d["nsite"].size
            S = len(labels)
            psum = np.zeros((nw, S))
            pn = np.zeros(nw)
            if binned["bins"].size:
                wi = (binned["bins"] * PR.BIN_BP) // W
                np.add.at(psum, wi, post)
                np.add.at(pn, wi, 1)
            for w in range(nw):
                a, b = w * W, min((w + 1) * W, c["length"])
                cells = [c["name"], str(a), str(b), str(d["nsite"][w]), str(d["cov2"][w])]
                for k in range(len(rnames)):
                    cv, mm = d["cov"][k, w], d["mm"][k, w]
                    cells += [str(cv), f"{1 - mm / cv:.4f}" if cv else "."]
                for k in range(S):
                    cv, mm = d["ss_cov"][k, w], d["ss_mm"][k, w]
                    cells += [str(cv), f"{1 - mm / cv:.4f}" if cv else "."]
                cells += [f"{psum[w, k] / pn[w]:.3f}" if pn[w] else "." for k in range(S)]
                fh.write("\t".join(cells) + "\n")
