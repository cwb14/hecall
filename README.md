# hecall: homoeologous exchanges in polyploid genome assemblies

hecall finds stretches of a polyploid chromosome whose sequence comes from the "wrong" subgenome.
These are homoeologous exchanges (HEs): a piece of an A chromosome replaced by its B homoeolog,
or swapped with it. Each call gets:

- a quality score
- breakpoints at the finest resolution the sequence allows
- an independent check against the homoeologous chromosome

hecall also says where it cannot tell.

## Install

```bash
git clone <this repository> hecall && cd hecall
mamba env create -f environment.yml && mamba activate hecall
pip install -e .        # optional: adds the `hecall` command; otherwise run `python -m hecall` here
pytest                  # 25 tests, ~3 min; the end-to-end tests need FastGA on PATH
```

`environment.yml` installs Python (≥ 3.9; numpy, scipy, pandas, matplotlib, pytest) and the FASTGA
suite from bioconda (`fastga=1.5.20260729`, i.e. FASTGA commit 178afcb, the version hecall was
validated with). hecall calls FastGA, ALNchain, ALNtoPAF, PAFtoALN, FAtoGDB, GIXmake, GIXrm and
GDBshow from PATH, or from the directory named by `HECALL_FASTGA`.

**Recommended:** build hecall's two C helpers against the FASTGA source of the same version
(a C compiler and zlib headers are needed):

```bash
git clone https://github.com/thegenemyers/FASTGA /path/to/FASTGA
git -C /path/to/FASTGA checkout 178afcb
bash csrc/build.sh /path/to/FASTGA
```

It writes into `bin/` of this repository:

- `bin/alnfilter`. Without it, hecall falls back to a slower PAF round trip.
- `bin/ALNchain`, the upstream ALNchain with a small speed fix (`csrc/alnchain_ties.patch`).
  Upstream ALNchain can run for hours when thousands of repeat alignments end at the same base,
  e.g. at a contig edge. Results are the same up to ties between equal-scoring chains.

hecall looks for `bin/` next to its package directory, so run it from the clone (`python -m hecall`
or `pip install -e .`). The log says when a helper is missing.

## What it needs

- **The polyploid assembly, split into subgenomes:** one FASTA per subgenome, each with a label.
  - Chromosome names must be unique across all subgenome files.
  - Give the subgenome FASTAs either all gzipped or all plain.
- **Reference genomes for the donor lineages:** two or more FASTAs (`.gz` is fine).
  - Diploid progenitors are best.
  - Stand-in relatives work.
  - A polyploid parent can be used once split into its subgenomes.
  - References don't need labels: hecall never assumes reference X *is* donor X.

## Run

```bash
hecall -s A=polyploid_A.fa.gz -s B=polyploid_B.fa.gz \
       -r diploidA.fa.gz -r diploidB.fa.gz -o out/ -t 16
```

- `-r NAME=path` names a reference (default: the file name); names appear in the output columns.
- Reruns resume from checkpoints in `out/work/`; a rerun with a new hecall version reuses the
  alignments.
- Alignment dominates the run time. At 8 threads, < 15 min for genomes under 0.5 Gb, ~30 min for
  *Coffea arabica* (1.2 Gb, two references), ~3 h for tobacco (4 Gb). Reruns from existing
  alignments took 3 min (*Coffea*) to 20 min (8-Gb *Aegilops*).
- Memory: the largest process peaked at a few GB for 1-Gb genomes and at 18–19 GB for 4–11-Gb
  genomes (8 threads).
- `out/work/` keeps the genome databases and chained alignments: about 0.8 GB for *Coffea* and
  2.7 GB for tobacco. Delete it when no rerun is needed.

## How it works

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

The model, every fixed constant and the reason for each change between versions are in
[`docs/design.md`](docs/design.md).

## Outputs

| file | content |
|---|---|
| `he.tsv`, `he.gff3` | every candidate HE: position, donor, QUAL, FILTER, type, evidence |
| `segments.bed` | full origin painting of every chromosome |
| `nocall.bed` | regions where an HE could not have been detected, and why |
| `chromosomes.tsv` | per chromosome: % of each origin, % not evaluable, % aligned to each other subgenome (`pct_homoeolog_<label>`; low = no homoeologue in the assembly: LOH, aneuploidy, lost or unassembled arm), number of HEs, label conflicts |
| `resolution.tsv` | per subgenome pair: evidence per window, window size, **smallest detectable HE**, genome-wide homoeologue SNV divergence (`ss_div_baseline`, reported even when there are no calls), the recipient label's own-origin share; serial correlation |
| `profiles.tsv` | diagnostic-site pattern frequencies per subgenome (used to place breakpoints) |
| `windows.tsv.gz` | every window: compared bases and divergence to each reference, origin axis, posterior |
| `painting.pdf/.png` + `painting_legend.md` | chromosome painting with every PASS HE marked |
| `hecall.log`, `run.json`, `hecall.done` | run log, inputs, and the version and call counts of the finished run |
| `work/` | checkpoints: genome databases, chained alignments, alignment tables, per-chromosome sums |

### Key columns of `he.tsv`

- **`qual`**: Phred-scaled posterior that the segment has the called donor origin, after
  discounting correlated windows (`llr / tau_n`) and a breakpoint prior. 20 = 99%.
- **`filter`**:
  - `PASS`
  - `LowQual` (QUAL < 20)
  - `unresolvable_pair` (the data cannot tell these two donors apart)
  - `donor_misfit`: the segment's divergence to some reference is outside what donor-origin
    sequence of that length shows (two-sided, every reference). Too far from all references means
    paralogous or repeat alignments, or a haplotype missing from the references; too close to the
    wrong reference means native sequence pulled into a broad foreign state. `p_donor_fit` is the
    test's p-value; the null is runs of the same number of native windows of the donor subgenome.
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
- **`div_<reference>`, `closest_ref`, `ref_target`**: the segment's SNV divergence to each
  reference, the closest one, and where the matching sequence lies in it.
- **`homoeolog`, `homoeolog_frac_recipient_origin`, `ss_div`, `ss_cover`**: the homoeologous
  interval, the share of it painted with the recipient's origin (≥ 0.5 = reciprocal), and the
  divergence to it and the fraction of the segment it covers.
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

## Validation

- **Synthetic genome (pytest).** These are planted:
  - reciprocal and non-reciprocal exchanges of 6–60 kb, found within 0.5 kb of their true ends
    (1.5 kb with a stand-in reference) and typed correctly
  - a 3-kb inter-subgenome transposition (decoy) and a clean chromosome, which give no call
  - a stand-in reference 3% from its donor, which still gives every call
  - a missing homoeologue, which is reported
- **41 polyploid assemblies from 38 complexes** (36 Mb to 11 Gb; diploid hybrid to 12x; plants,
  fish, an oyster hybrid and a fungal hybrid), run with v0.4.10:
  - **30 of 31 published events reproduced:**
    - the *A. suecica* chr1 tip swap (Jiang et al. 2021)
    - the chr1 B/C swaps in *Oryza minuta* and *O. malampuzhaensis* (Fornasiero et al. 2025)
    - the *Aegilops ventricosa* 1N/3D arm exchange (Liu et al. 2025)
    - both *Solanum demissum* events, with the right reciprocal/non-reciprocal type (Hosaka & Hosaka 2026)
    - four tobacco junctions (Zan et al. 2025)
    - the *Mimulus sookensis* read-depth blocks (Whitener et al. 2024)
    - the *Saccharum* chr02 *S. officinarum* block (Wang et al. 2026)
    - 9 of 10 banana segments (Xie et al. 2024). The miss is *Musa* 'Silk' chr03B, only 14% called.
  - **Traps:** a *Glycine* chr8 rearrangement inherited from the diploid parent is not called. The
    phased haplotypes of *Rosa* 'Carola' give no PASS call, and hecall reports that its labels are
    not subgenomes.
  - **Published negatives** carry no large exchange:
    - *Camelina microcarpa*: 3 PASS calls, ≤ 51 kb
    - *Brachypodium hybridum*: 1 call, 32 kb
    - *Prunus cerasus*: 3 calls, ≤ 144 kb
    - *A. kamchatica*: 21 calls, ≤ 0.5 Mb, 0.5% of the genome. Its published negative comes from
      read coverage.
  - **Totals:**
    - *Arachis hypogaea* has 1.9% of its genome in PASS calls (Yin et al. 2020: 2.54%).
    - In synthetic *A. thaliana* × *A. arenosa* Allo738, direction and amounts are within 2× of
      Jiang et al. 2021.
  - **Resolution.** The smallest detectable HE has a median of 28 kb across genomes (quartiles
    10–237 kb). It reaches Mb only where the subgenomes are barely farther apart than the
    references (*Saccharum*, kiwifruit, *Magallana*, *Avena agadiriana*).
  - These genomes were also used to develop hecall, so this is a regression check, not an
    independent test.
- **Where it is weakest:**
  - Fragmented assemblies: *A. monticola* (2018, 42,618 gaps) gets 884 PASS calls against 26 in
    its near-identical sister *A. hypogaea*.
  - Stand-in references equidistant from both subgenomes (*Avena agadiriana*).
  - Clonal triploids with many small calls (*Musa* 'Silk').
  - Check such calls with reads.

## Limits

- HEs are detected relative to the references. Exchanges older than the split between the
  references are invisible.
- One sequence per subgenome means a heterozygous HE is seen only if it sits in the assembled
  haplotype.
- A misjoin between homoeologous chromosomes looks like a reciprocal HE. Check calls flagged
  `gap_at_breakpoint` with reads.
- Repeat-rich regions that do not align to two or more references are not evaluable (`nocall.bed`).
- Only SNVs are used; indels and structural variants carry no origin evidence here.

## Credits

hecall runs the FASTGA suite (FastGA, ALNchain, ALNtoPAF and others; ONElib) by Gene Myers,
Chenxi Zhou and Richard Durbin: <https://github.com/thegenemyers/FASTGA> (BSD-style licence).
`csrc/alnchain_ties.patch` modifies FASTGA's `ALNchain.c`, and `csrc/alnfilter.c` is compiled
against FASTGA sources.

## License

MIT (see `LICENSE`). The FASTGA suite has its own BSD-style licence.
