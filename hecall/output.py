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


def write_all(outdir, labels, refs, chroms, segs, he, nocall, chrom_rows, fit, wins, res, slogp,
              pattern_names, data):
    rnames = [r[0] for r in refs]
    S = len(labels)
    cols = ["chrom", "start", "end", "length", "label", "donor", "qual", "filter", "type", "core", "core_qual",
            "scale", "llr", "tau_n", "p_identity", "p_donor_fit", "log_odds", "n_windows", "n_sites", "n_support", "n_against", "bp_left_ci", "bp_right_ci",
            "syntenic", "closest_ref", "ref_target"] + [f"div_{n}" for n in rnames] + \
           ["homoeolog", "homoeolog_frac_recipient_origin", "ss_div", "ss_div_baseline", "ss_cover", "flags"]
    with open(os.path.join(outdir, "he.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for h in sorted(he, key=lambda h: (h["chrom"], h["start"])):
            row = dict(h)
            row["closest_ref"] = rnames[h["closest_ref"]]
            row["syntenic"] = h["syntenic"] or "."
            for n, v in zip(rnames, h["ref_div"]):
                row[f"div_{n}"] = v
            fh.write("\t".join(row[c] if isinstance(row.get(c), str) else _f(row.get(c)) for c in cols) + "\n")
    with open(os.path.join(outdir, "he.gff3"), "w") as fh:
        fh.write("##gff-version 3\n")
        for i, h in enumerate(sorted(he, key=lambda h: (h["chrom"], h["start"])), 1):
            attr = (f"ID=HE{i};donor={h['donor']};recipient={h['label']};qual={h['qual']};"
                    f"filter={h['filter']};type={h['type']};n_windows={h['n_windows']};flags={h['flags']}")
            fh.write(f"{h['chrom']}\thecall\thomoeologous_exchange\t{h['start'] + 1}\t{h['end']}\t"
                     f"{h['qual']}\t.\t.\t{attr}\n")
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
    p = np.exp(slogp)
    with open(os.path.join(outdir, "profiles.tsv"), "w") as fh:
        fh.write("# diagnostic-site patterns (used to place breakpoints): one symbol per reference ("
                 + ",".join(rnames) + "); M=matches polyploid, .=uncovered, a/b=distinct other alleles\n")
        fh.write("pattern\t" + "\t".join(f"freq_{l}" for l in labels) + "\n")
        for j in np.argsort(-p.max(axis=0)):
            fh.write(pattern_names[j] + "\t" + "\t".join(f"{p[k, j]:.5f}" for k in range(S)) + "\n")
    with open(os.path.join(outdir, "resolution.tsv"), "w") as fh:
        fh.write("recipient\tdonor\tllr_per_window\twindow_bp\ttau_at_min\tmin_he_bp\tmedian_chrom_bp\tresolvable"
                 "\tss_div_baseline\trecipient_own_origin\n")
        for r in res:
            fh.write(f"{labels[r['recipient']]}\t{labels[r['donor']]}\t{r['llr_per_window']:.3f}\t"
                     f"{r['window_bp']:.0f}\t{_f(r['tau_at_min'])}\t{_f(r['min_he_bp'], 6)}\t"
                     f"{r['median_chrom_bp']:.0f}\t{'yes' if r['resolvable'] else 'no'}\t"
                     f"{_f(r.get('ss_div_baseline'))}\t{_f(r.get('recipient_own_origin'))}\n")
        fh.write(f"# breakpoint rate (per bp): segmentation HMM {fit['rho']:.4g}; calibrated prior {fit['rho_cal']:.4g}\n")
        fh.write("# serial correlation of window evidence (native windows), lag: rho = "
                 + ", ".join(f"{k}:{fit['rho_ac'][k]:.2f}" for k in (1, 2, 5, 10, 20, 50, 100) if k < fit['rho_ac'].size) + "\n")
    with gzip.open(os.path.join(outdir, "windows.tsv.gz"), "wt", compresslevel=3) as fh:
        dim = wins[0]["X"].shape[1] if wins else 0
        fh.write("chrom\tstart\tend\t" + "\t".join(f"compared_{n}\tdiv_{n}" for n in rnames) + "\t"
                 + "\t".join(f"origin_axis{j + 1}" for j in range(dim)) + "\t"
                 + "\t".join(f"post_{l}" for l in labels) + "\tcalled\n")
        for w, post in zip(wins, fit["post"]):
            for i in range(w["start"].size):
                cells = [w["name"], str(int(w["start"][i])), str(int(w["end"][i]))]
                for k in range(len(rnames)):
                    cv, mm = int(w["cov"][i, k]), int(w["mm"][i, k])
                    cells += [str(cv), f"{mm / cv:.5f}" if cv else "."]
                cells += [f"{x:.4f}" for x in w["X"][i]]
                cells += [f"{x:.4f}" for x in post[i]]
                cells.append(labels[int(post[i].argmax())])
                fh.write("\t".join(cells) + "\n")
