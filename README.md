# hecall: homoeologous exchanges in polyploid genome assemblies

hecall finds stretches of a polyploid chromosome whose sequence comes from the "wrong" subgenome.
These are homoeologous exchanges (HEs): a piece of an A chromosome replaced by its B homoeolog,
or swapped with it. Each call gets:

- a quality score
- breakpoints at the finest resolution the sequence allows
- an independent check against the homoeologous chromosome

hecall also says where it cannot tell.

## What it needs

- **The polyploid assembly, split into subgenomes:** one FASTA per subgenome, each with a label.
- **Reference genomes for the donor lineages:** two or more FASTAs.
  - Diploid progenitors are best.
  - Stand-in relatives work.
  - A polyploid parent can be used once split into its subgenomes.
  - References don't need labels: hecall never assumes reference X *is* donor X.
- **Software:**
  - FastGA suite on PATH (or `HECALL_FASTGA=/path/to/FASTGA`)
  - Python ≥ 3.9 with numpy, scipy, pandas and matplotlib
- **Recommended for multi-Gb genomes:** build the small `alnfilter` helper against the FASTGA
  source tree: `bash csrc/build.sh /path/to/FASTGA`, which writes `bin/alnfilter`.
  - Without it, hecall falls back to a slower PAF round trip.

## How it works (one paragraph)

Every reference genome is aligned to the polyploid with FastGA, keeping only syntenic chains
(ALNchain). Along each polyploid chromosome, hecall builds windows that each hold the same
information: enough compared bases for ~10 mismatches to the closest reference. In each window it
measures how much closer the polyploid is to one reference than to another. How these
measurements are spread for "A-donor" and "B-donor" sequence is learned from the polyploid's own A
and B subgenomes. Native quirks (incomplete lineage sorting, divergent reference haplotypes,
hard-to-align regions) are therefore part of what is expected, not evidence of exchange. A hidden
Markov model proposes stretches whose origin differs from their chromosome's subgenome. Their
breakpoints are then placed between single diagnostic bases. Each candidate's evidence is
discounted for the correlation between neighbouring windows that native sequence shows, before it
gets a quality score. The polyploid's own homoeologue gives an independent check: a swapped piece
(reciprocal), or two near-identical copies (non-reciprocal).

## Run

```bash
export PATH=/path/to/FASTGA:$PATH
python -m hecall -s A=polyploid_A.fa.gz -s B=polyploid_B.fa.gz \
                 -r diploidA.fa.gz -r diploidB.fa.gz -o out/ -t 16
```

Reruns resume from checkpoints in `out/work/`.

## Outputs

| file | content |
|---|---|
| `he.tsv`, `he.gff3` | every candidate HE: position, donor, QUAL, FILTER, type, evidence |
| `segments.bed` | full origin painting of every chromosome |
| `nocall.bed` | regions where an HE could not have been detected, and why |
| `chromosomes.tsv` | per chromosome: % of each origin, % not evaluable, number of HEs, label conflicts |
| `resolution.tsv` | per subgenome pair: evidence per window, window size, **smallest detectable HE**; serial correlation |
| `profiles.tsv` | diagnostic-site pattern frequencies per subgenome (used to place breakpoints) |
| `windows.tsv.gz` | every window: compared bases and divergence to each reference, origin axis, posterior |
| `painting.pdf/.png` + `painting_legend.md` | chromosome painting with every PASS HE marked |

### Key columns of `he.tsv`

- **`qual`**: Phred-scaled posterior that the segment has the called donor origin, after
  discounting correlated windows (`llr / tau_n`) and a breakpoint prior. 20 = 99%.
- **`filter`**:
  - `PASS`
  - `LowQual` (QUAL < 20)
  - `unresolvable_pair` (the data cannot tell these two donors apart)
  - `nonsyntenic` (the matching reference sequence is not at the homoeologous position:
    transposed TE, translocation or misassembly)
- **`type`**:
  - `reciprocal`: the homoeologue carries the swapped piece
  - `non-reciprocal`: both copies are from the donor; `ss_div` ≪ `ss_div_baseline` confirms it
  - `unpaired`: no alignable homoeologue
- **`llr`, `tau_n`, `n_windows`**: summed window evidence, its serial-correlation discount, and the
  number of windows.
- **`n_sites`, `n_support`, `n_against`**: diagnostic sites inside the segment and how many favour
  the donor.
- **`p_identity`**: for non-reciprocal calls, the fraction of same-length genome blocks whose
  homoeologues are at least as similar. Small values back a recent non-reciprocal exchange; the
  value is conservative.
- **`bp_left_ci`, `bp_right_ci`**: 95% credible interval of each breakpoint, or `terminal`.
- **`flags`**:
  - `gap_at_breakpoint`: an assembly gap lies within a breakpoint interval (a possible misjoin)
  - `terminal`: the segment reaches a chromosome end

## Limits

- HEs are detected relative to the references. Exchanges older than the split between the
  references are invisible.
- One sequence per subgenome means a heterozygous HE is seen only if it sits in the assembled
  haplotype.
- A misjoin between homoeologous chromosomes looks like a reciprocal HE. Check calls flagged
  `gap_at_breakpoint` with reads.
- Repeat-rich regions that do not align to two or more references are not evaluable (`nocall.bed`).
