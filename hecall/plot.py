"""Chromosome painting with homoeologous exchanges: painting.pdf / painting.png (+ legend .md)."""
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle, Patch  # noqa: E402

from . import params as PR  # noqa: E402

# Okabe-Ito order; first four pass all-pairs CVD dE >= 8 (dataviz validator). Further hues are
# only legal with direct labels (column headers name each subgenome).
ORIGIN_COLORS = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9"]
UNCERTAIN = "#dcdcdc"


def _tint(hexc, f):
    c = np.array(matplotlib.colors.to_rgb(hexc))
    return tuple(c + (1 - c) * f)


def painting(outdir, labels, chroms, segs, he, nocall, max_rows=None):
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
                         "pdf.fonttype": 42, "svg.fonttype": "none", "font.size": 7})
    S = len(labels)
    cols = {lab: [c for c in chroms if c["label"] == lab] for lab in labels}
    nrow = max(len(v) for v in cols.values())
    maxlen = max(c["length"] for c in chroms)
    color = {i: ORIGIN_COLORS[i % len(ORIGIN_COLORS)] for i in range(S)}
    he_by = {}
    for h in he:
        he_by.setdefault(h["chrom"], []).append(h)
    row_h = 0.22
    fig_h = 0.9 + nrow * row_h
    fig_w = 7.2
    fig, axes = plt.subplots(1, S, figsize=(fig_w, fig_h), sharey=True, squeeze=False)
    axes = axes[0]
    bar = 0.5
    for j, lab in enumerate(labels):
        ax = axes[j]
        for r, c in enumerate(cols[lab]):
            y = nrow - 1 - r
            for s in segs[c["name"]]:
                a, b = s["tile_start"], s["tile_end"]
                if s.get("nosites"):
                    fc = "white"
                elif s["qual"] < PR.QUAL_PASS or s["state"] != c["label_idx"]:
                    fc = UNCERTAIN  # foreign stretches are drawn saturated only where a PASS call covers them
                else:
                    fc = _tint(color[s["state"]], 0.62)
                ax.add_patch(Rectangle((a, y - bar / 2), b - a, bar, fc=fc, ec="none", lw=0))
            for h in he_by.get(c["name"], []):
                d = labels.index(h["donor"])
                if h["filter"] == "PASS":
                    ax.add_patch(Rectangle((h["start"], y - bar / 2), h["end"] - h["start"], bar,
                                           fc=color[d], ec="none", lw=0, zorder=3))
                    ax.plot([(h["start"] + h["end"]) / 2], [y + bar / 2 + 0.17], marker="v", ms=3.2,
                            color=color[d], mec="none", zorder=4, clip_on=False)
            for a, b, _ in nocall[c["name"]]:
                ax.add_patch(Rectangle((a, y - bar / 2), b - a, bar, fc="white", ec="none", zorder=2.5))
                ax.add_patch(Rectangle((a, y - bar / 2), b - a, bar, fc="none", ec="#c8c8c8", lw=0,
                                       hatch="xxxx", zorder=2.6))
            ax.add_patch(Rectangle((0, y - bar / 2), c["length"], bar, fc="none", ec="#4d4d4d", lw=0.4, zorder=5))
            ax.text(-maxlen * 0.015, y, c["name"], ha="right", va="center", fontsize=5.5)
        ax.set_xlim(0, maxlen * 1.01)
        ax.set_ylim(-0.8, nrow - 0.2 + 0.25)
        ax.set_yticks([])
        for sp in ("left", "right", "top"):
            ax.spines[sp].set_visible(False)
        ax.spines["bottom"].set_color("#4d4d4d")
        ax.spines["bottom"].set_linewidth(0.5)
        ticks = np.arange(0, maxlen + 1, _tick_step(maxlen))
        ax.set_xticks(ticks)
        ax.set_xticklabels([f"{t / 1e6:g}" for t in ticks], fontsize=6)
        ax.set_xlabel("Position (Mb)", fontsize=7)
        ax.tick_params(axis="x", length=2, width=0.5, color="#4d4d4d")
        ax.text(0, nrow - 0.2 + 0.25, f"Subgenome {lab}", fontsize=7.5, ha="left", va="bottom")
    handles = [Patch(fc=_tint(color[i], 0.62), ec="none", label=f"{lab} origin") for i, lab in enumerate(labels)]
    handles += [Patch(fc=color[i], ec="none", label=f"exchange from {lab} (PASS)") for i, lab in enumerate(labels)]
    handles += [Patch(fc=UNCERTAIN, ec="none", label="origin uncertain / not PASS"),
                Patch(fc="white", ec="#c8c8c8", hatch="xxxx", lw=0.4, label="not evaluable")]
    fig.legend(handles=handles, loc="lower center", ncol=min(4, len(handles)), fontsize=6, frameon=False,
               bbox_to_anchor=(0.5, 0.0), handlelength=1.6, columnspacing=1.2)
    fig.subplots_adjust(left=0.12, right=0.99, top=1 - 0.35 / fig_h, bottom=(0.55 + 0.1 * len(handles) / 4) / fig_h,
                        wspace=0.35)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"painting.{ext}"), dpi=300)
    plt.close(fig)
    _legend(outdir, labels, chroms, he, nocall)


def _tick_step(L):
    for s in (1e6, 2e6, 5e6, 10e6, 20e6, 50e6, 100e6, 200e6):
        if L / s <= 8:
            return s
    return 500e6


def _legend(outdir, labels, chroms, he, nocall):
    p = [h for h in he if h["filter"] == "PASS"]
    tot = sum(c["length"] for c in chroms)
    nc = sum(b - a for v in nocall.values() for a, b, _ in v)
    lens = sorted(h["length"] for h in p)
    med = f"{np.median(lens) / 1e3:.1f} kb (range {lens[0] / 1e3:.1f}–{lens[-1] / 1e3:.1f} kb)" if lens else "n/a"
    types = {t: sum(1 for h in p if h["type"] == t) for t in ("reciprocal", "non-reciprocal", "unpaired")}
    with open(os.path.join(outdir, "painting_legend.md"), "w") as fh:
        fh.write(
            "# painting\n\n"
            "Origin = the subgenome donor whose sequence a stretch of chromosome most resembles, judged from its "
            "divergence to every reference genome in windows of equal information, against the spread seen in "
            "the polyploid's own subgenomes. A homoeologous exchange (HE) is a stretch "
            "whose origin differs from the chromosome's own subgenome. PASS = posterior probability ≥ 0.99 "
            "(QUAL ≥ 20), syntenic with the donor-like reference, and a donor/recipient pair the data can tell "
            "apart.\n\n"
            f"Each polyploid chromosome is drawn to scale, grouped by subgenome ({', '.join(labels)}). Pale fill "
            "shows the inferred origin along the chromosome; saturated blocks are PASS HEs, coloured by donor, "
            "each marked by a triangle so that exchanges too small to see at this scale remain visible. Light "
            "grey marks stretches whose origin is uncertain or whose exchange call did not pass; grey cross-hatching marks regions "
            "too sparsely aligned to evaluate. "
            f"{len(p)} PASS HEs ({types['reciprocal']} reciprocal, {types['non-reciprocal']} non-reciprocal, "
            f"{types['unpaired']} without an assignable homoeologue); median length {med}. "
            f"{100 * nc / tot:.1f}% of the assembly was not evaluable.\n")
