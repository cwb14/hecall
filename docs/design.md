# hecall design note

Date: 2026-10-08. Status: approved design (chat, 2026-10-08); model revised after real-data testing (v0.2 windows, v0.3 calibration; section 6).

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

### 2.1 Windows of equal information (v0.2)

- Per-base codes are summed in 1-kb blocks: compared bases and SNV mismatches per reference.
- A window grows over blocks until every reference has ≥ cmin compared bases. cmin gives ~10
  expected mismatches to the closest reference at the genome's typical divergence, so the window's
  log-divergence ratio has SE ~0.3–0.45.
- A window needing more than 10× the median span is not formed. Such regions go to nocall.
- **Feature:**
  - centred log-divergences to the K references (Helmert contrasts, K−1 dimensions)
  - projected by LDA onto ≤ S−1 dimensions when K−1 > S−1
  - Centring removes the local mutation-rate or constraint level shared by all references.
- **Emission:** a gridded KDE of the feature, one per state, fitted to windows assigned to that
  state. The fit starts from chromosome labels; windows called non-native are reassigned and the
  KDEs re-fitted (a few rounds).
  - Native heterogeneity is therefore inside each state's density: incomplete lineage sorting,
    introgression between progenitors, divergent reference haplotypes, hard-to-align regions.
- **Diagnostic sites** are used only to place breakpoints between windows. A diagnostic site has
  ≥ 2 covered references that disagree, and the polyploid matches at least one of them.

### 2.2 HMM and segmentation

- **Transitions:** continuous-time uniform switching at rate ρ (closed form over the distance
  between window centres).
- **Fitting:** ρ and the start distribution are fitted by EM. Forward-backward gives window
  posteriors.
- **Segments** are runs of the max-posterior state. Each boundary is refined at diagnostic-site
  level inside the two flanking windows:
  - P(cut) ∝ exp(cumulative site log-likelihood)
  - the output is the ML cut and a 95% credible interval
  - cuts are processed left to right and never cross
- **Terminal segments** stop at the outermost evaluated window.

### 2.3 Calibrated quality (v0.3)

Window evidence is serially correlated: shared-ancestry blocks and regional divergence span many
windows. In *Coffea* the correlation is 0.60 at lag 1, 0.20 at 10 and 0.06 at 100 windows. The
HMM therefore over-counts evidence and is used only to propose candidates. Each candidate is
scored with:

- **LLR:** Σ over its n windows of [log f_t − log f_s] (KDE densities).
- **τₙ:** 1 + 2 Σₖ₌₁ⁿ⁻¹ (1 − k/n) ρₖ, from native windows' evidence autocorrelation (Geyer
  initial-positive truncation). This is the composite-likelihood / effective-sample-size
  correction for a sum of n correlated terms: ≈3 for 5 windows, ≈15 for 100 in *Coffea*.
- **log posterior odds:** LLR/τₙ + n_breakpoints · log(ρ_c · window span / (S−1)).
  - ρ_c is empirical Bayes: the breakpoint rate implied by the PASS calls, iterated to a fixed
    point.
- **QUAL** = Phred(1 − posterior). PASS requires all of:
  - QUAL ≥ 20
  - a resolvable donor/recipient pair
  - not `nonsyntenic`
- **Resolution:** the smallest n with n·KL/τₙ + 2·log(ρ_c·span/(S−1)) ≥ ln 99. KL is the mean
  evidence per window of donor sequence. A pair is unresolvable if no n up to a median chromosome
  achieves this.

### 2.4 Independent evidence

- **Homoeologue (self alignment).** Each candidate is projected onto its homoeologous interval in
  subgenome *t* (robust weighted percentiles of aligned bases) and classified as:
  - `reciprocal`: that interval has origin *s*
  - `non-reciprocal`: it keeps origin *t*. Both copies come from one donor, and `ss_div` ≪
    `ss_div_baseline` confirms it.
  - `unpaired`: no homoeologous alignment
- `p_identity` = fraction of same-span genome blocks at least as similar between homoeologues.
  - It is reported, not added to QUAL: near-identical homoeologous blocks are common genome-wide
    (1% of 10-kb *Coffea* blocks have 0 mismatches), most likely conversion tracts themselves, so
    that null is contaminated.
- **Synteny.** The closest reference's sequence matching the candidate must be collinear with the
  sequence matching at least one native flank of its non-native block. Positions are
  aligned-base-weighted medians.
  - A transposed TE fails, since neither flank agrees.
  - So does a non-homoeologous (translocation-like) segment.
  - An exchange next to a rearrangement between assemblies passes, because one flank agrees.

### 2.5 Artifact flags and no-calls

- **Flags:**
  - `gap_at_breakpoint`: an assembly gap lies inside a breakpoint CI. This may be a misjoin, which
    looks like a reciprocal exchange.
  - `terminal`
  - `homoeolog_not_closer` (v0.4.1): non-reciprocal, the donor homoeologue covers >= 25% of the
    segment, and its divergence from the segment is >= the genome-wide divergence between those
    subgenomes. It is an evidence flag, not a filter: the donor copy may be absent from the assembly.
    - In validated complexes, non-reciprocal PASS calls are much closer to their donor homoeologue
      than the baseline. Median ratio: Musa 0.01, Arachis 0.01, Nicotiana 0.02, Coffea 0.05,
      B. napus 0.20.
    - Reciprocal swaps sit near 1, as expected.
    - Sphagnum troendelagicum non-reciprocal calls: median 1.04 (10/18 >= 1). Actinidia: 22/136 >= 1.
- `label_conflict` appears in chromosomes.tsv.
- `pct_homoeolog_<label>` (v0.4.2) in chromosomes.tsv: % of the chromosome in collinear chains with
  each other subgenome. Low cover means there is no homoeologue in the assembly, so no exchange can
  be evaluated there (descriptive, no threshold). Causes include:
  - LOH in trio-binned or collapsed assemblies (Cryptococcus AD: chr03A and chr08D absent)
  - aneuploidy
  - lost or unassembled arms
- `nocall.bed`: stretches covered by no window, at least one median window span long. The reason
  is `no_alignment` or `sparse_alignment(<refs>)`.
- Subgenome pairs the data cannot tell apart are `unresolvable`: their calls are never PASS.

## 3. Alignment (benchmarked 2026-10-08)

Benchmark: *Coffea arabica* polyploid C+E (1.2 Gb) vs *C. canephora* (0.6 Gb), 8 threads.

| aligner | wall | peak RAM | chr01C covered by own parent (best hit/base) | SNV divergence |
|---|---|---|---|---|
| minimap2 2.30 -x asm20 -c --cs --secondary=no | 56.5 min | 112 GB | 61.6% | 1.92% |
| FastGA (defaults) + ALNchain (defaults) | 6.3 min (+0.5 min chain/PAF) | 2.4 GB | 60.6% | 1.94% |
| FastGA + ALNchain -g1M -p0.001 -z10k -s5000 | same | same | 92.0% | 2.90% |

- **Choice: FastGA + ALNchain with defaults.** Coverage and divergence are the same as minimap2's
  primary chains, at about 1/6 of the time and 1/40 of the memory.
- Relaxed chaining inflates coverage with repeat alignments. The extra coverage is spread over all
  *canephora* chromosomes at a median 90% identity, so it is not orthologous.
- About 40% of chr01C does not align to its own parent in any setting. These are TE-rich regions
  turned over since the species split; they become nocall regions.
- The user's alternatives:
  - Syntenic-chunk alignment is unnecessary at FastGA's speed.
  - sassy does pattern search, not genome alignment.
  - mm2-plus is not installed.
- Per base, the least divergent chained alignment supplies the reference allele. The orthologue or
  homoeologue is the closest copy; older paralogues are more diverged. Indel bases inside a better
  alignment are not filled from worse ones.
- Inter-subgenome transposed TEs are removed twice:
  - by ALNchain, as short isolated chains (synthetic decoy test)
  - by the per-HE synteny check (FILTER `nonsyntenic`)
- The self alignment (FastGA on the polyploid alone) reports both directions. It provides the
  homoeologue of every base.

## 4. Inputs / outputs

Input:
- polyploid subgenome FASTAs, each with a label: `-s A=polyA.fa.gz -s B=polyB.fa.gz`
- reference FASTAs: `-r refA.fa.gz -r refB.fa.gz …`, any number and any lineage, including polyploid
  parents split into subgenomes

Output (one directory per polyploid):

- `he.tsv`: every non-native segment with QUAL, LLR, site counts, breakpoint intervals, type and FILTER
- `he.gff3`: the same as GFF3
- `segments.bed`: the full origin painting, native segments included
- `nocall.bed`: unevaluable regions and the reason for each
- `chromosomes.tsv`: per chromosome:
  - label and majority origin
  - % evaluable
  - number of HEs and HE bp
- `profiles.tsv`: the emission profiles
- `resolution.tsv`: for each subgenome pair, the KL divergence, site density, L* and whether the pair
  is resolvable
- `bins.tsv.gz`: per 10-kb window, identity to each reference, site counts and posterior (for
  plotting and checks)
- `painting.pdf` / `painting.png`:
  - each chromosome drawn to scale, coloured by origin
  - every PASS HE marked so that small HEs stay visible
  - nocall regions in grey

## 5. Validation

- pytest on a synthetic genome with:
  - planted reciprocal and non-reciprocal HEs of several sizes
  - a TE transposed between subgenomes (decoy)
  - a clean chromosome (true negative)
  - a stand-in reference (divergent donor)
- one real-data check against published HE coordinates (literature: `_hecall_batch/literature/known_HEs.md`)
- one quick run on a dissimilar genome

## 6. Decisions log

- 2026-10-08 v0.1 → v0.2: the per-site i.i.d. HMM over diagnostic and private sites failed on
  *Coffea*: 5,671 PASS calls, and the breakpoint rate sat at its upper bound.
  - Native C sequence matches the E reference at 20% of diagnostic sites, in correlated blocks
    (1-kb lag-1 autocorrelation 0.73; variance 8–130× binomial).
  - Polyploid-private sites pile up in fast-evolving regions.
  - Context classes (conditioning on closest-reference divergence) were tried and dropped: with a
    stand-in reference the class becomes confounded with origin (synthetic stand-in test failed).
  - The window model learned from the polyploid's own subgenomes replaced it.
- 2026-10-08 v0.3.8: multi-scale proposals.
  - Segmentation runs with window evidence divided by 1, sqrt(tau) and tau, where tau is the
    integrated autocorrelation time of native-window evidence. All scales are scored identically
    (llr / tau_n + breakpoint prior).
  - Selection: PASS candidates longest first (the full supported extent), then the rest by
    log-odds; no overlaps. Each call keeps its best-scoring core.
  - The empirical-Bayes breakpoint rate now starts from the lower bound, so it converges to the
    smallest self-consistent rate.
  - Why: Saccharum chr02S's O-derived block (66–97 Mb; published So ~82–95) was shattered into
    failing fragments at raw scale (only 2 Mb PASS), because 30% of its windows look S-like in
    ~100-kb runs. Now it is called whole (QUAL 24, core 71–73 Mb at QUAL 31).
  - Regression checks:
    - A. suecica chr1 swap is unchanged; chr01Ar is now one block, and the small weak calls are gone.
    - Actinidia: 318 -> 199 PASS. The calls look donor-like (donor-ref divergence 0.7% vs native
      0.9%).
- 2026-10-09 v0.3.9: tie-safe ALNchain (speed only).
  - Upstream ALNchain builds a KD-tree per scaffold pair by quickselect with a two-way partition.
    When thousands of alignments share an end coordinate, every partition step removes one element,
    and the build becomes quadratic or worse.
  - Digitaria sanguinalis 6x C vs E: 79,597 prefiltered alignments end at chr02E:36,854,183, a
    repeat array cut by a contig edge. ALNchain ran > 80 min on a 28-MB input.
  - Fix: csrc/alnchain_ties.patch (three-way partition: equal keys form a middle block), built by
    csrc/build.sh into bin/ALNchain and preferred by find_tools. Digitaria: 1.2 s. Coffea
    polyploid vs C. canephora: 23 s -> 4 s.
  - Output is the same up to tie-breaks among equal-scoring chains: 0.06% of aligned bp differ on
    Coffea, all between tandem copies.
  - Regression test: 15,000 alignments ending at one target base (upstream > 120 s; patched < 1 s).
- 2026-10-09 v0.4.0: each state's emission density is learned from native windows only.
  - Before: in each round, the density of state s was fitted on every window assigned to s,
    including windows on other subgenomes' chromosomes. Those windows are the candidate exchanges
    being tested, so they were used both to train and to test.
  - Failure: windows about equally far from all references (repeats aligned to paralogues) sit midway
    on the origin axis. On Digitaria 4x, chr04D 21.7–25.9 Mb is 13.3% from the D reference and 13.9%
    from E. In round 0 these went to whichever state had the heavier tail (E). E's density then grew
    there, and the block ended at posterior E = 1.0 (segments.bed QUAL 31; 10% of chr04D listed as E).
  - Now: s's density uses windows on s-labelled chromosomes that are currently assigned to s. The
    labelled subgenomes are the training set; candidate exchanges are only tested.
  - Effect, from refitting the model on each batch run's windows (% of window bp confidently foreign:
    all / windows whose closest reference is within 20% of the second closest):

    | polyploid | before | after |
    |---|---|---|
    | Digitaria 4x | 1.56 / 18.3 | 1.05 / 10.2 (chr04D block: own-origin posterior 0.00 -> 0.80) |
    | Poa annua | 0.66 / 10.1 | 0.56 / 6.2 |
    | B. napus | 1.90 / 7.1 | 1.83 / 6.5 |
    | GR-A3n | 15.6 / 30.1 | 14.9 / 29.4 |
    | A. suecica, O. minuta, Musa Plantain, Musa Silk, Coffea, Magallana, Actinidia | unchanged (±0.1) | |

  - All validated exchanges are kept: A. suecica chr1 swap, O. minuta chr1 swap, Musa Plantain 6/6,
    Musa Silk 4/4, Saccharum chr02S, Poa chr07B, Digitaria chr08D.
  - Not changed: a no-re-estimation variant (labels only) lowers foreign painting further, but every
    real exchange then also enters its recipient's density, which costs power in exchange-rich
    genomes (Actinidia, 16%).
- 2026-10-09 v0.4.4: extraction robustness and thread budget.
  - Multi-threaded `ALNtoPAF -s` segfaults intermittently (an upstream race, FASTGA 178afcb). Rosa
    hap3 vs hap4 crashes with -T2 to -T8; -T1 always works and gives the same records. A failed
    threaded stream is now redone with one thread, with a warning.
  - Each extraction worker is a parser process plus ALNtoPAF. ALNtoPAF now gets at most half a
    worker's share of `-t` (previously a fixed 2 threads per worker, so up to ~3x oversubscription).
- 2026-10-09 v0.4.5: label coherence (issue #2).
  - The model assumes each label is one subgenome of one origin, with exchanges a minority.
    Coherence is the share of a label's evaluable sequence painted with its own origin. Below one
    half, the label is not a subgenome, and calls involving it get FILTER `incoherent_label`.
  - Batch evidence (own-origin share per label):
    - Rosa 'Carola' (phased haplotypes, numbered per chromosome): 0.22–0.46 for all four labels
    - Magallana (subgenomes not assigned to species): A 0.27, G 0.85
    - GR-A3n (auratus/gibelio homologue recombination): Au 0.30, Cu 0.94, G 0.91
    - every other complex: >= 0.71 for every label. Actinidia (polysomic) is borderline at
      0.53–0.64 and is reported, not filtered.
  - With two labels the share rarely falls below 0.5 by construction (the majority origin becomes
    the label's state). It falls in practice when a label mixes several origins or the states
    overlap heavily.
- 2026-10-09 v0.4.6: no lag cap on the autocorrelation of native-window evidence.
  - The 2000-window cap (TAU_MAX_LAG) was arbitrary. Where native evidence wanders over megabases,
    the cap understated the correlation of long blocks and overstated their significance.
    Examples: stand-in references, introgression, recombination-poor cereal chromosomes. Geyer's
    initial positive sequence now decides the truncation alone (FFT, so any lag is cheap).
  - Avena agadiriana (proxy references ~4.5% from both subgenomes): Geyer cut 2000 -> 13,328 lags,
    integrated tau 281 -> 652. Large-call log-odds fall 2–4×; 4 of 11 PASS fall below 0.
  - Cost: Saccharum chr02S 66–97 Mb (published So block, QUAL 24 before) is expected to fall to
    log-odds ~ -0.3 (LowQual). Kept anyway: a defensible rule over a tuned one. With Saccharum's
    long-range native wander, the data alone do not separate a 31-Mb block from it.
  - Calls spanning <= 2000 windows are unchanged.
- 2026-10-09 v0.4.7: typical closest-reference divergence uses blocks that double from 10 kb
  until the median is non-zero.
  - Avena magna x longiglumis synthetic hexaploid: its references are its parents, so most 10-kb
    blocks hold no mismatch to the closest reference. The median was 0 and the run stopped with
    "nothing to tell apart", which is the opposite of the truth.
  - No other genome changes (all have a non-zero 10-kb median), so v0.4.6 outputs equal v0.4.7.
- 2026-10-09 v0.4.8: window size from the native-vs-homoeologous contrast.
  - Before: ~10 expected mismatches to the closest reference. Windows were huge when the references
    are the parents (Avena magna x longiglumis: d1 = 1.6e-6, 29-Mb windows) or very close to them
    (A. hypogaea: 80-kb windows). Those are the cases where origin is clearest.
  - Now: c = WINDOW_Z2 (d1 + dh) / (dh - d1)^2 compared bases per reference. d1 is the typical
    closest-reference divergence; dh is the typical divergence between the polyploid's own
    homoeologues (an exchanged segment sits ~dh from its native reference). Native and foreign
    windows then differ by ~3.2 Poisson SDs at typical divergences.
  - Label-free: a second-nearest-reference version failed with several references per subgenome
    (sister references), so the homoeologue divergence is used instead.
  - Effect (compared bases per window, old -> new):
    - finer where the references are close: Avena magna 6.4M -> 225, A. hypogaea 33,324 -> 191,
      Poa 2,564 -> 188, Brachypodium 2,271 -> 162
    - coarser where subgenomes are barely farther apart than the references: Magallana 523 -> 10,288,
      Avena agadiriana 280 -> 9,336, Actinidia 660 -> 3,942, Saccharum 463 -> 2,801
    - little change for most others (1-kb block floor)
  - Validated swaps are unchanged offline (A. suecica, O. minuta, Musa Plantain, Fragaria).
  - dh <= d1 (Rosa haplotypes: 0.22% between, 0.51% to the references): origin cannot be told apart.
    Windows fall back to the closest-reference rule for painting, and every pair is unresolvable.
- 2026-10-09 v0.4.10: donor fit checks every reference in both directions (Bonferroni over 2K),
  with 4999 draws.
  - In v0.4.9 only "farther than donor sequence from the closest reference" was tested. Avena magna
    D chromosomes then kept "A-origin" calls that are closest to the D reference (1.4% to D vs 4.1%
    to A). The A reference is not the exact A parent, so A-native sequence is itself 3–5% from
    every reference, and the one-sided test let these through.
- 2026-10-09 v0.4.9: donor fit.
  - A call must look like sequence of the donor's origin. Null: runs of the same number of
    consecutive native windows of the donor subgenome (serial correlation kept). Statistic: the
    divergence to the closest reference. p < 0.01 gives FILTER `donor_misfit`.
  - Why: emission densities only see contrasts, so windows far from every reference (paralogues,
    repeats, haplotypes absent from the references) can land in a foreign state's tail.
    - Avena magna x longiglumis (references = parents): 7 calls are 0.00005–0.25% from the donor
      parent, i.e. textbook exchanges. About 36 are 3–14% from every reference.
    - Seen earlier: Digitaria chr04D, Raphanobrassica C3.
  - A crude per-window p99 version flagged 41/1331 PASS calls, concentrated in Avena magna (22/43)
    and A. hypogaea (17/54).
- 2026-10-08 v0.3.5: synteny check has three classes:
  - `homoeologous`
  - `translocation`: one collinear chain on another chromosome or position, spanning at least half
    the segment and >= 100 kb
  - `nonsyntenic`, filtered
  - Why: *N. tabacum*'s published reciprocal intergenomic translocation (Nt21/Nt22 = our chr06T /
    chr09S, 32 + 81 Mb) was found at QUAL 99 but filtered. *Ae. ventricosa*'s 1N/3D arm exchange is
    the same class, so the change rests on 2 complexes.
- 2026-10-08 v0.2 → v0.3: a single global τ (integrated autocorrelation, 116 windows in *Coffea*)
  over-tempered short segments. It was replaced by the length-specific τₙ per candidate and an
  empirical-Bayes breakpoint prior. *Coffea*: 149 → 15 PASS.

- 2026-10-08:
  - States = polyploid subgenomes, not references (handles missing and stand-in donors).
  - SNVs only.
  - Profiles are anchored to labels (no free EM).
  - ρ is fitted by EM.
  - The user suggested aligning only syntenic sequence for speed, and pointed to FastGA, sassy and
    mm2-plus. sassy does approximate pattern search, not genome alignment, so it was not considered.
    mm2-plus is not installed. FastGA is built at ~/bin/src/FASTGA and is being benchmarked.
