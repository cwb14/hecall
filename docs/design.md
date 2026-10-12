# hecall design note

Describes hecall 0.4.11: the model, the fixed choices and their reasons, and what changed in each
version and why. User-facing usage and outputs are in the README.

## 1. Goal

Call homoeologous exchanges (HEs) in a polyploid assembly at the finest resolution the data allow,
with a confidence value on every call, and say explicitly where no call can be made.

- An HE here = a stretch of a chromosome of subgenome *s* whose sequence descends from the donor of
  another subgenome *t*.
- HEs are detectable only relative to the references given. Exchanges older than the split between
  the references are invisible (e.g. ancient strawberry HEs shared by all octoploids).

## 2. Model (the one idea)

Each polyploid region is compared with every reference genome. The HMM states are the polyploid's
own subgenomes. What "sequence from donor *t*" looks like is learned from the polyploid's subgenome
*t* itself, not from the reference labelled *t*.

This is needed because references are often not the true donor: stand-ins, divergent accessions,
polyploid parents, missing donors. A raw rule "closer to reference A than to B" fails in those
cases. A learned profile does not need the donor to be present; it needs the references to *tell
the donors apart*.

### 2.1 Alignment and per-base comparison

- Each reference is aligned to the whole polyploid, and each pair of subgenomes to each other
  (the homoeologues), with FastGA. See section 3 for the benchmark.
- **Collinear-run prefilter:** raw alignments are grouped by sequence pair, strand and 10-kb
  diagonal band. Runs of alignments ≤ 10 kb apart with ≥ 5 kb aligned are kept. This is a superset
  of what ALNchain can keep, so only isolated repeat hits are removed. Without it, ALNchain
  stalls on TE-rich multi-Gb genomes.
- **ALNchain** (defaults) keeps one-to-one syntenic chains. A transposed TE between subgenomes
  forms a short isolated chain and is mostly removed here.
- **Per-base codes:** match, mismatch (with the reference allele) or uncovered, from the PAF `cs`
  tag. Where alignments of one reference overlap, each base takes its code from the least
  divergent alignment: the orthologue or homoeologue is the closest copy, older paralogues are
  more diverged. Indel bases inside a better alignment are not filled from worse ones.
- Only SNVs are used. Codes are summed per 1-kb block (compared bases and mismatches per
  reference and per homoeologue).

### 2.2 Windows of equal information

- A window grows over blocks until every reference has ≥ c compared bases:
  c = WINDOW_Z2 · (d1 + dh) / (dh − d1)².
  - d1 = typical divergence to the closest reference.
  - dh = typical divergence between the polyploid's homoeologues. An exchanged segment sits about
    dh from its native reference.
  - With WINDOW_Z2 = 10, native and foreign windows differ by ~3.2 Poisson SDs in mismatch count.
  - "Typical" = the median over blocks that start at 10 kb and double until the median is
    non-zero.
- If dh ≤ d1, the subgenomes are no farther apart than the references. Origin cannot be told
  apart: windows fall back to ~10 mismatches to the closest reference for painting, and every
  pair is unresolvable.
- A window needing more than 10× the median span is not formed. Such regions go to `nocall.bed`
  (reason `no_alignment` or `sparse_alignment(<references>)`).

### 2.3 Features and emissions

- **Feature:** centred log-divergences to the K references (Helmert contrasts, K−1 dimensions).
  Centring removes the local mutation-rate or constraint level shared by all references.
  - When K−1 > S−1 (S = number of subgenomes), LDA projects onto S−1 dimensions.
- **Emission:** a gridded Gaussian KDE per state.
  - State *s*'s density is fitted only on windows on *s*-labelled chromosomes that are currently
    assigned to *s*. Candidate exchanges are tested, never used to train.
  - Rounds: fit densities → HMM → reassign windows → refit (up to 3 re-estimation rounds).
  - Native heterogeneity is therefore inside each state's density: incomplete lineage sorting,
    introgression between progenitors, divergent reference haplotypes, hard-to-align regions.
- **Diagnostic sites** are used only to place breakpoints between windows. A diagnostic site has
  ≥ 2 covered references that disagree, and the polyploid matches at least one of them.

### 2.4 HMM and segmentation

- **Transitions:** continuous-time uniform switching at rate ρ (closed form over the distance
  between window centres). ρ and the start distribution are fitted by EM.
- **Segments** are runs of the max-posterior state. Each boundary is refined at diagnostic-site
  level inside the two flanking windows:
  - P(cut) ∝ exp(cumulative site log-likelihood)
  - the output is the ML cut and a 95% credible interval
  - cuts are processed left to right and never cross
- **Multi-scale proposals:** the segmentation is repeated with window evidence divided by √τ and
  τ, where τ is the integrated autocorrelation time of native-window evidence (when ≥ 2).
  - Raw evidence finds small exchanges.
  - Tempered evidence keeps long exchanges with serially correlated noisy windows whole.
  - Candidates from all scales are scored identically. Selection: PASS candidates longest first
    (the full supported extent), then the rest by log-odds; no overlaps. Each call keeps its
    best-scoring candidate inside it as `core`.
- **Terminal segments** stop at the outermost evaluated window.

### 2.5 Calibrated quality

Window evidence is serially correlated: shared-ancestry blocks and regional divergence span many
windows. In *Coffea* the correlation is 0.60 at lag 1, 0.20 at 10 and 0.06 at 100 windows. The
HMM therefore over-counts evidence and is used only to propose candidates. Each candidate is
scored with:

- **LLR:** Σ over its n windows of [log f_t − log f_s] (KDE densities).
- **τₙ:** 1 + 2 Σₖ₌₁ⁿ⁻¹ (1 − k/n) ρₖ, from the autocorrelation of native windows' evidence.
  - Truncation is Geyer's initial positive sequence alone, with no lag cap.
  - This is the composite-likelihood / effective-sample-size correction for a sum of n
    correlated terms: ≈3 for 5 windows and ≈15 for 100 in *Coffea*.
- **log posterior odds:** LLR/τₙ + n_breakpoints · log(ρ_c · window span / (S−1)).
  - ρ_c is empirical Bayes: the breakpoint rate implied by the PASS calls, iterated to a fixed
    point from the lower bound. Several self-consistent rates can exist; this finds the smallest
    one, i.e. the most conservative rate the calls support.
- **QUAL** = Phred(1 − posterior), capped at 99.
- **Resolution:** the smallest n with n·KL/τₙ + 2·log(ρ_c·span/(S−1)) ≥ ln 99. KL is the mean
  evidence per window of donor sequence. A pair is unresolvable if no n up to a median chromosome
  achieves this.

### 2.6 Filters (anything but PASS)

- `LowQual`: QUAL < 20 (posterior 0.99, a conventional evidence level).
- `unresolvable_pair`: section 2.5, or dh ≤ d1 (section 2.2).
- `incoherent_label`: a label less than half of whose evaluable sequence is painted with its own
  origin is not a subgenome; calls involving it are filtered.
  - The model assumes each label is one subgenome of one origin, with exchanges a minority.
  - With two labels the share rarely falls below one half by construction (the majority origin
    becomes the label's state). It falls when a label mixes several origins or the states
    overlap heavily, e.g. phased haplotypes numbered per chromosome.
- `donor_misfit`: the call must look like sequence of the donor's origin against every reference.
  - Null: runs of the same number of consecutive native windows of the donor subgenome (serial
    correlation kept), 4999 draws.
  - Per reference, the candidate's divergence is compared with the null in both directions.
    Bonferroni over 2K tests; p < 0.01 filters.
  - It catches windows far from every reference (paralogous or repeat alignments, haplotypes
    missing from the references) and native windows pulled into a broad foreign state. Emission
    densities see only contrasts, so both can land in a foreign state's tail.
- `nonsyntenic`: section 2.7.

### 2.7 Independent evidence, flags and no-call statements

- **Homoeologue typing.** Each candidate is projected onto its homoeologous interval in
  subgenome *t*, through the pairwise subgenome alignment. The projection uses the dominant
  collinear chain only, and robust aligned-base-weighted percentiles. It is classified as:
  - `reciprocal`: that interval is mostly of origin *s*
  - `non-reciprocal`: it keeps origin *t*. Both copies come from one donor, and `ss_div` ≪
    `ss_div_baseline` confirms it.
  - `unpaired`: no homoeologous alignment on subgenome *t* covers ≥ 25% of the segment
- **`p_identity`** = fraction of same-span genome blocks at least as similar between homoeologues.
  - It is reported, not added to QUAL: near-identical homoeologous blocks are common genome-wide
    (1% of 10-kb *Coffea* blocks have 0 mismatches), most likely conversion tracts themselves, so
    that null is contaminated.
- **Synteny.** The closest reference's sequence matching the candidate is compared with the
  sequence matching the native flanks of its non-native block (flanks of max(block, 100 kb), at
  most 5 Mb). Positions are aligned-base-weighted medians along the dominant chain.
  - `homoeologous`: on the same reference chromosome as at least one flank, within
    1 Mb + 3 × the query distance (local TE turnover and small rearrangements).
  - `translocation`: not collinear with a flank, but one collinear block elsewhere spanning at
    least half the segment and ≥ 100 kb. Reported, not filtered.
  - `nonsyntenic`: neither (transposed repeat, scattered hits, misassembly). Filtered.
- **Flags** (annotation, not filters):
  - `gap_at_breakpoint`: an assembly gap lies inside a breakpoint CI. This may be a misjoin, which
    looks like a reciprocal exchange.
  - `terminal`
  - `homoeolog_not_closer`: non-reciprocal, the donor homoeologue covers ≥ 25% of the segment, and
    its divergence from the segment is ≥ the genome-wide divergence between those subgenomes.
    A copy should be much closer; the donor copy may be absent from the assembly.
- **No-call statements:**
  - `nocall.bed` (section 2.2)
  - `unresolvable_pair` and the smallest detectable HE per pair (`resolution.tsv`)
  - each label's own-origin share (`recipient_own_origin`)
  - `pct_homoeolog_<label>` in `chromosomes.tsv`: % of the chromosome in collinear chains with
    each other subgenome. Low cover means there is no homoeologue in the assembly (LOH in
    trio-binned or collapsed assemblies, aneuploidy, lost or unassembled arms); descriptive, no
    threshold.

## 3. Alignment (benchmarked 2026-10-08)

Benchmark: *Coffea arabica* polyploid C+E (1.2 Gb) vs *C. canephora* (0.6 Gb), 8 threads.

| aligner | wall | peak RAM | chr01C covered by own parent (best hit/base) | SNV divergence |
|---|---|---|---|---|
| minimap2 2.30 -x asm20 -c --cs --secondary=no | 56.5 min | 112 GB | 61.6% | 1.92% |
| FastGA (defaults) + ALNchain (defaults) | 6.3 min (+0.5 min chain/PAF) | 2.4 GB | 60.6% | 1.94% |
| FastGA + ALNchain -g1M -p0.001 -z10k -s5000 | same | same | 92.0% | 2.90% |

- **Choice: FastGA + ALNchain with defaults.** Coverage and divergence are the same as minimap2's
  primary chains, at about 1/8 of the time and 1/45 of the memory.
- Relaxed chaining inflates coverage with repeat alignments. The extra coverage is spread over all
  *canephora* chromosomes at a median 90% identity, so it is not orthologous.
- About 40% of chr01C does not align to its own parent in any setting. These are TE-rich regions
  turned over since the species split; they become nocall regions.
- Not used:
  - Aligning only syntenic chunks is unnecessary at FastGA's speed.
  - sassy does approximate pattern search, not genome alignment.
  - mm2-plus was not benchmarked.

## 4. Fixed constants and remaining heuristics

Every fixed constant is in `hecall/params.py` with a one-line reason. Values learned from each
input are the emission densities, breakpoint rates, window size, serial correlation and
resolution. None of the constants was tuned to make known answers come out right.

Constants set by convention or by a tool's own defaults:

- QUAL 20
- donor-fit p < 0.01
- 95% breakpoint intervals
- the prefilter, which mirrors ALNchain's defaults
- WINDOW_Z2 = 10

Heuristics that remain (reasoned, not derived):

- synteny tolerance 1 Mb + 3× distance; flanks 100 kb–5 Mb; translocation block ≥ 100 kb
- dominant-chain tolerance max(1 Mb, half the interval)
- homoeologue cover ≥ 25% to type a call
- majority rules at one half (reciprocal typing, label coherence)
- the window span cap (10× median)
- nocall reasons (< 5% aligned = `no_alignment`; a reference below half the best = sparse)
- tempering only when τ ≥ 2

## 5. Validation

- pytest on a synthetic allotetraploid, among others (`tests/`):
  - planted reciprocal and non-reciprocal HEs of several sizes
  - a TE transposed between subgenomes (decoy)
  - a clean chromosome (true negative)
  - a stand-in reference (divergent donor)
  - a dropped homoeologue
  - the many-ties ALNchain case
  - unit tests of every model stage
- Real data: 41 polyploids from 38 complexes against published events and negatives; the README
  summarises it. These genomes were also used during development, so the check is a regression
  check, not an independent test.

## 6. Version history

Each entry: the change, the evidence that prompted it, and its effect. Changes to the model were
made only on a principle, not to match a known answer, and were checked on at least two complexes.

- **0.1.0** (2026-10-08): per-site i.i.d. HMM over diagnostic and polyploid-private sites.
- **0.2.0** (2026-10-08): window model learned from the polyploid's own subgenomes.
  - Why: v0.1 failed on *Coffea*: 5,671 PASS calls, and the breakpoint rate sat at its upper
    bound.
    - Native C sequence matches the E reference at 20% of diagnostic sites, in correlated blocks
      (1-kb lag-1 autocorrelation 0.73; variance 8–130× binomial).
    - Polyploid-private sites pile up in fast-evolving regions.
  - Context classes (conditioning on closest-reference divergence) were tried and dropped: with a
    stand-in reference the class becomes confounded with origin (synthetic stand-in test failed).
- **0.3.0**: length-specific τₙ per candidate and an empirical-Bayes breakpoint prior, replacing
  one global τ (116 windows in *Coffea*) that over-tempered short segments. *Coffea*: 149 → 15
  PASS. `p_identity` added as an annotation.
- **0.3.1**: projections use the dominant collinear chain. Repeat hits no longer smear homoeologue
  intervals; the *Oryza* chr1 swap is now typed reciprocal.
- **0.3.2**: collinear-run prefilter before ALNchain; raw alignments checkpointed. Equivalent on
  *Coffea* (99.4–99.6% identical per-base codes). It removes the ALNchain bottleneck on TE-rich
  multi-Gb genomes.
- **0.3.3**: homoeologues from pairwise subgenome alignments instead of one self alignment of the
  whole polyploid, which took ~45% of the run time.
- **0.3.4**: `alnfilter` (C) copies the prefiltered alignment records directly. Before, every raw
  alignment was converted to CIGAR PAF and back: tens of minutes on large raw files, now seconds.
- **0.3.5**: synteny classes `homoeologous` / `translocation` / `nonsyntenic`.
  - Why: tobacco's published reciprocal intergenomic translocation (Nt21/Nt22 = chr06T / chr09S,
    32 + 81 Mb) was found at QUAL 99 but filtered as nonsyntenic. The *Ae. ventricosa* 1N/3D arm
    exchange is the same class, so the change rests on 2 complexes.
- **0.3.6**: fixed a hang in pairwise mode.
  - Swapped (`ALNtoPAF -w`) PAF arrives interleaved; it is now sorted by query before extraction.
  - Fatal errors are exceptions, so a failing worker stops the run instead of hanging it.
- **0.3.7**: process pools for extraction and sites: a worker that dies (e.g. out of memory) stops
  the run with an explained error.
- **0.3.8**: multi-scale proposals (section 2.4) and the conservative (lowest fixed point)
  empirical-Bayes rate.
  - Why: *Saccharum* chr02S's O-derived block (66–97 Mb; published *S. officinarum* ~82–95 Mb) was
    shattered into failing fragments at raw scale (only 2 Mb PASS), because 30% of its windows
    look S-like in ~100-kb runs. Now it is called whole.
  - The *A. suecica* chr1 swap was unchanged. Kiwifruit went from 318 to 199 PASS; the remaining
    calls look donor-like (donor-reference divergence 0.7% vs native 0.9%).
- **0.3.9**: tie-safe ALNchain (speed only).
  - Upstream ALNchain builds a KD-tree per scaffold pair by quickselect with a two-way partition.
    When thousands of alignments share an end coordinate, every partition step removes one
    element, and the build becomes quadratic or worse.
  - *Digitaria sanguinalis* 6x C vs E: 79,597 prefiltered alignments end at one base of a repeat
    array cut by a contig edge. ALNchain ran > 80 min on a 28-MB input.
  - Fix: `csrc/alnchain_ties.patch`, a three-way partition in which equal keys form a middle
    block. *Digitaria* now takes 1.2 s, and *Coffea* polyploid vs *C. canephora* 4 s (was 23 s).
  - Output is the same up to tie-breaks among equal-scoring chains: 0.06% of aligned bp differ on
    *Coffea*, all between tandem copies.
- **0.4.0**: each state's emission density is learned from native windows only (section 2.3).
  - Before: the density of state *s* was fitted on every window assigned to *s*, including
    windows on other subgenomes' chromosomes. Those are the candidate exchanges being tested, so
    they were used both to train and to test.
  - Failure: windows about equally far from all references (repeats aligned to paralogues) sit
    midway on the origin axis. On *Digitaria* 4x, chr04D 21.7–25.9 Mb is 13.3% from the D
    reference and 13.9% from E. These windows went to E, E's density grew there, and the block
    ended at posterior E = 1.0.
  - Effect, as % of window bp confidently foreign: *Digitaria* 4x 1.56 → 1.05 (that block's
    own-origin posterior 0.00 → 0.80); *Poa* 0.66 → 0.56; *B. napus* 1.90 → 1.83; GR-A3n 15.6 →
    14.9; seven other genomes unchanged (±0.1). All validated exchanges were kept.
  - Not chosen: fitting densities on labels alone (no re-estimation). It lowers foreign painting
    further, but every real exchange then also enters its recipient's density, which costs power
    in exchange-rich genomes (kiwifruit, 16%).
- **0.4.1**: flag `homoeolog_not_closer`.
  - In validated complexes, non-reciprocal PASS calls are much closer to their donor homoeologue
    than the baseline. Median ratio: *Musa* 0.01, *Arachis* 0.01, tobacco 0.02, *Coffea* 0.05,
    *B. napus* 0.20.
  - *Sphagnum troendelagicum* (median 1.04) and kiwifruit (22/136 ≥ 1) are exceptions, so the
    flag marks them.
- **0.4.2**: `pct_homoeolog_<label>` in `chromosomes.tsv`. Missing homoeologues are stated, e.g.
  *Cryptococcus* AD chr03A and chr08D (LOH in a trio-binned assembly).
- **0.4.3**: `ss_div_baseline` per pair in `resolution.tsv`, also when there are no calls.
- **0.4.4**: extraction robustness.
  - Multi-threaded `ALNtoPAF -s` segfaults intermittently (an upstream race, FASTGA 178afcb).
    *Rosa* hap3 vs hap4 crashes with -T2 to -T8; -T1 always works and gives the same records.
    A failed threaded stream is now redone with one thread, with a warning.
  - ALNtoPAF gets at most half of each extraction worker's share of `-t`. Before, it had a fixed
    2 threads per worker, i.e. up to ~3× oversubscription.
- **0.4.5**: label coherence, `incoherent_label` (section 2.6). Own-origin share per label:
  - *Rosa* 'Carola' (phased haplotypes, numbered per chromosome): 0.22–0.46 for all four labels
  - *Magallana* (subgenomes not assigned to species): A 0.27, G 0.85
  - GR-A3n (auratus/gibelio homologue recombination): Au 0.30, Cu 0.94, G 0.91
  - every other complex: ≥ 0.71 for every label. Kiwifruit (polysomic) is 0.53–0.64: reported,
    not filtered.
- **0.4.6**: no lag cap on the autocorrelation of native-window evidence.
  - The 2000-window cap was arbitrary. Where native evidence wanders over megabases, it
    understated the correlation of long blocks and overstated their significance. Examples:
    stand-in references, introgression, recombination-poor cereal chromosomes.
  - *Avena agadiriana* (proxy references ~4.5% from both subgenomes): Geyer's cut moved from 2000
    to 13,328 lags, and the integrated τ from 281 to 652. Large-call log-odds fell 2–4×; 4 of 11
    PASS calls fell below 0.
  - Accepted cost: the *Saccharum* chr02S block was expected to fall to LowQual (a defensible rule
    over a tuned one). With v0.4.8's coarser windows for *Saccharum* it is PASS again (65.5–97.1
    Mb, QUAL 42).
- **0.4.7**: the typical closest-reference divergence uses blocks that double from 10 kb until the
  median is non-zero.
  - Why: in the *Avena magna* × *A. longiglumis* synthetic hexaploid the references are its
    parents, so most 10-kb blocks hold no mismatch. The median was 0, and the run stopped with
    "nothing to tell apart", which is the opposite of the truth.
  - No other genome changed.
- **0.4.8**: window size from the native-vs-homoeologous contrast (section 2.2).
  - Before: ~10 expected mismatches to the closest reference. Windows were huge where the
    references are the parents (*Avena magna*: d1 = 1.6e-6, 29-Mb windows) or very close
    (*A. hypogaea*: 80-kb windows). Those are the cases where origin is clearest.
  - Label-free: a second-nearest-reference version failed with several references per subgenome
    (sister references), so the homoeologue divergence is used instead.
  - Effect (compared bases per window, old → new):
    - finer where the references are close: *Avena magna* 6.4M → 225, *A. hypogaea* 33,324 → 191,
      *Poa* 2,564 → 188, *Brachypodium* 2,271 → 162
    - coarser where subgenomes are barely farther apart than the references: *Magallana* 523 →
      10,288, *Avena agadiriana* 280 → 9,336, kiwifruit 660 → 3,942, *Saccharum* 463 → 2,801
  - Validated swaps unchanged (*A. suecica*, *O. minuta*, *Musa* Plantain, *Fragaria*).
  - dh ≤ d1 (*Rosa* haplotypes: 0.22% between, 0.51% to the references) → unresolvable.
- **0.4.9**: donor fit, one-sided (farther than donor sequence from the closest reference).
  - Why: *Avena magna* × *A. longiglumis* (references = parents). 7 calls are 0.00005–0.25% from
    the donor parent, textbook exchanges; about 36 are 3–14% from every reference. Seen earlier
    in *Digitaria* chr04D and *Raphanobrassica* C3.
- **0.4.10**: donor fit checks every reference in both directions (section 2.6).
  - Why: *Avena magna* D chromosomes kept "A-origin" calls closest to the D reference (1.4% to D
    vs 4.1% to A). The A reference is not the exact A parent, so A-native sequence is itself 3–5%
    from every reference, and the one-sided test let these through.
- **0.4.11**: cleanup for release; calls unchanged.
  - Removed resume support for pre-0.3.3 work directories that hold a whole-polyploid self
    alignment.
  - Alignments whose polyploid and reference sequence share a name are no longer skipped
    (a leftover of the self alignment).
  - The log says when the C helpers are not built.
  - The painting legend lists every PASS condition.
  - Added `pyproject.toml` and `environment.yml`.
  - Regression: byte-identical outputs (legend text aside) on *Cryptococcus* AD and synthetic
    Allo738. *A. suecica*, rerun from pairwise subgenome alignments instead of its v0.4.10 self
    alignment, gives the same calls and QUALs with slightly different homoeologue statistics.
