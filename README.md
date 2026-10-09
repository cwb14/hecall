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
- **Recommended:** build hecall's two C helpers against the FASTGA source tree:
  `bash csrc/build.sh /path/to/FASTGA`. It writes:
  - `bin/alnfilter`. Without it, hecall falls back to a slower PAF round trip.
  - `bin/ALNchain`, the upstream ALNchain with a small speed fix (`csrc/alnchain_ties.patch`).
    Upstream ALNchain can run for hours when thousands of repeat alignments end at the same base,
    e.g. at a contig edge.

## How it works (one paragraph)

Every reference genome is aligned to the polyploid with FastGA, keeping only syntenic chains
(ALNchain); the polyploid's subgenomes are aligned to each other the same way. Along each polyploid
chromosome, hecall builds windows that each hold the same information. A window is long enough
that native sequence and sequence from the other subgenome are expected to differ clearly in how
far they are from the references. In each window it measures how much closer the polyploid is to
one reference than to another. How these measurements are spread for each origin is learned only
from that subgenome's own chromosomes. Native quirks (incomplete lineage sorting, divergent
reference haplotypes, hard-to-align regions) are therefore part of what is expected, not evidence
of exchange.

A hidden Markov model proposes stretches whose origin differs from their chromosome's subgenome.
Their breakpoints are then placed between single diagnostic bases. Each candidate's evidence is
discounted for the correlation between neighbouring windows that native sequence shows, before it
gets a quality score. A candidate must also look like sequence of the donor's origin against every
reference, not merely be unlike its own. The polyploid's own homoeologue gives an independent
check: a swapped piece (reciprocal), or two near-identical copies (non-reciprocal).

hecall also says when no call can be made: regions that align to too few references (`nocall.bed`),
subgenome pairs it cannot tell apart, and labels that are not real subgenomes (e.g. phased
haplotypes).

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
| `chromosomes.tsv` | per chromosome: % of each origin, % not evaluable, % aligned to each other subgenome (`pct_homoeolog_<label>`; low = no homoeologue in the assembly: LOH, aneuploidy, lost or unassembled arm), number of HEs, label conflicts |
| `resolution.tsv` | per subgenome pair: evidence per window, window size, **smallest detectable HE**, genome-wide homoeologue SNV divergence (`ss_div_baseline`, reported even when there are no calls); serial correlation |
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
  - `donor_misfit`: the segment is farther from every reference than donor-origin sequence ever
    is at that length, so it is not a copy of the donor. Typical causes are paralogous or
    repeat alignments, or a haplotype missing from the references. `p_donor_fit` is the test's
    p-value; the null is runs of the same number of native windows of the donor subgenome.
  - `incoherent_label`: the recipient or donor label is less than half its own origin, so it is not
    a subgenome. Typical causes are haplotypes numbered independently per chromosome, or
    subgenomes the input mis-assigned. resolution.tsv gives each label's own-origin share
    (`recipient_own_origin`), and the log warns.
  - `nonsyntenic` (the matching reference sequence is neither at the homoeologous position nor
    one long collinear block elsewhere: transposed repeat, scattered hits or misassembly)
- **`type`**:
  - `reciprocal`: the homoeologue carries the swapped piece
  - `non-reciprocal`: both copies are from the donor; `ss_div` ≪ `ss_div_baseline` confirms it
  - `unpaired`: no alignable homoeologue
- **`syntenic`**:
  - `homoeologous`: the donor-like sequence is collinear with the native flanks, i.e. an exchange
    at the homoeologous position
  - `translocation`: a long collinear block from another chromosome or position, i.e. an
    intergenomic exchange between non-homoeologous sites, such as the tobacco Nt21/Nt22
    translocation
  - `nonsyntenic`
- **`core`, `core_qual`**: the most concentrated well-supported stretch inside the call. Long
  exchanges are reported whole; `core` gives where the evidence peaks.
- **`scale`**: the segmentation that proposed the call:
  - `x1`: raw window evidence, which finds small exchanges
  - `xN`: window evidence divided by N, which keeps long exchanges with noisy windows whole
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
  - `homoeolog_not_closer`: a non-reciprocal call whose segment is no more similar to the
    polyploid's own donor homoeologue than homoeologues are genome-wide. A copy of that homoeologue
    should be much closer, so the reference signal may have another source: a donor copy missing
    from the assembly (other homologues in a higher polyploid), a stand-in reference, or paralogy.

## Limits

- HEs are detected relative to the references. Exchanges older than the split between the
  references are invisible.
- One sequence per subgenome means a heterozygous HE is seen only if it sits in the assembled
  haplotype.
- A misjoin between homoeologous chromosomes looks like a reciprocal HE. Check calls flagged
  `gap_at_breakpoint` with reads.
- Repeat-rich regions that do not align to two or more references are not evaluable (`nocall.bed`).
