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
- `label_conflict` appears in chromosomes.tsv.
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
