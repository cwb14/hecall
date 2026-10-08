# hecall design note

Date: 2026-10-08. Status: approved design (chat, 2026-10-08); implemented in hecall 0.1.0.

## 1. Goal

Call homoeologous exchanges (HEs) in a polyploid assembly at the finest resolution the data allow,
with a confidence value on every call, and say explicitly where no call can be made.

- An HE here = a stretch of a chromosome of subgenome *s* whose sequence descends from the donor of
  another subgenome *t*.
- HEs are detectable only relative to the references given. Exchanges older than the split between
  the references are invisible (e.g. ancient strawberry HEs shared by all octoploids).

## 2. Model (the one idea)

Each polyploid base is compared with every reference genome. The states of an HMM are the
polyploid's own subgenomes, and what "sequence from donor *t*" looks like is learned from the
polyploid's subgenome *t* itself, not from the reference labelled *t*.

Why: references are often not the true donor (stand-ins, divergent accessions, polyploid parents,
missing donors). A raw rule "closer to reference A than to B" fails in those cases. Sphagnum B,
for example, is closer to the T reference than to its own B reference. A learned profile does not
need the donor to be present. It only needs the references to tell the donors apart.

### 2.1 Sites

- Every reference genome is aligned to the whole polyploid, which is the query. Each polyploid base
  then gets at most one aligned base from each reference.
- Bases covered by more than one alignment from the same reference are dropped as ambiguous.
- An **informative site** is a polyploid base that ≥ 2 references cover and that at least one covered
  reference does not match.
- Only SNVs are used. Indels are not, because a 5-kb TE indel would otherwise count as 5,000
  correlated sites.
- **Pattern**: one symbol per reference:
  - `M`: matches the polyploid
  - `.`: uncovered
  - `a`, `b`, …: the distinct non-matching alleles, in order of appearance
- With two references, the pattern is one of:
  - `Ma`: matches reference 1 only
  - `aM`: matches reference 2 only
  - `aa`: polyploid-private
  - `ab`: triallelic

### 2.2 Emission profiles

- **P(pattern | state t)** = the frequency of that pattern among informative sites on the chromosomes
  labelled *t*, with a pseudocount.
- Refinement: HE segments contaminate the native profile of the subgenome that carries them. After
  the first pass, bins called non-native are reassigned to their called state and the profiles are
  re-estimated. This is repeated until the calls stop changing (at most a few rounds).
- Profiles are never free-fitted (full Baum-Welch). With free fitting, a state could drift to model
  gene-rich versus gene-poor regions instead of origin.
- Patterns too rare to estimate are pooled into one `rare` category (threshold in params.py).

### 2.3 HMM

- Sites are pooled into fixed-size bins. A bin's log emission is the sum over its sites. Empty bins
  emit nothing.
- Transition model: a continuous-time switch at rate ρ per bp, to any other state with equal chance.
  The transition over a gap of *d* bp uses the closed form:
  - P(stay) = 1/S + (S−1)/S · e^(−Sρd/(S−1))
  - This gives the exact cost of crossing long site-free gaps.
- ρ is estimated by EM (expected switches / total length) within fixed bounds.
- Forward-backward gives a posterior per bin. Segments are runs of the maximum-posterior state.

### 2.4 Calls and quality

For each non-native segment (state ≠ chromosome label):

- **QUAL** = Phred(1 − mean posterior of the called state over the segment's bins), capped at 99.
  PASS requires QUAL ≥ 20 (99%).
- **LLR** = Σ over sites of [log P(pattern|t) − log P(pattern|s)]. This is the evidence independent
  of ρ.
- n_sites, n_support (sites favouring t), n_against.
- **Breakpoints** are refined to site resolution. Within the boundary bins, P(cut after site i) ∝
  exp(cumulative LLR). The output is the maximum-likelihood cut plus a 95% credible interval, and
  the interval between the flanking informative sites is the resolution floor.

### 2.5 Independent evidence (homoeolog check)

The subgenomes are aligned to each other. Each HE segment on chromosome X is projected to its
homoeologous interval on the subgenome-*t* chromosome, and classified as:

- **reciprocal**: the homoeolog carries origin *s* there (the segments swapped, 2:2).
- **non-reciprocal**: the homoeolog keeps origin *t*. Both copies come from *t* (4:0 dosage), so the
  identity between the two copies across the segment should rise to within-donor level. That
  identity is reported against the genome-wide baseline for the pair.
- **unpaired**: there is no homoeologous alignment, so the type cannot be told.

Pre-polyploid introgression into a donor shows neither signature.

### 2.6 Artifact flags

- `gap_bp`: an assembly gap (N run) lies inside a breakpoint credible interval. This may be a misjoin,
  which cannot be told apart from a reciprocal HE without reads.
- `terminal`: the segment reaches a chromosome end (information only; crossover-derived HEs are often
  distal).
- `label_conflict`: the chromosome's majority origin differs from its label.
- Synteny: only alignments within collinear chains are used. TEs that transposed between subgenomes
  after polyploidy would otherwise mimic tiny HEs.

### 2.7 What is not called (and said so)

- `nocall.bed`: runs with no informative sites that are longer than the minimum detectable HE. The
  reason is one of no reference coverage, only one reference, or no differences.
- Resolution per subgenome pair:
  - expected LLR per kb = KL(profile_t ‖ profile_s) × site density
  - minimum detectable length L* = LLR needed for QUAL 20 at the fitted ρ ÷ LLR per kb
  - If L* exceeds the median chromosome length of *s*, the pair is `unresolvable`: its calls are
    reported with FILTER = unresolvable_pair, never PASS.
- No subgenome structure (e.g. Rosa 'Carola' haplotypes) shows up as unresolvable pairs; the tool
  does not force a call.

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

- 2026-10-08:
  - States = polyploid subgenomes, not references (handles missing and stand-in donors).
  - SNVs only.
  - Profiles are anchored to labels (no free EM).
  - ρ is fitted by EM.
  - The user suggested aligning only syntenic sequence for speed, and pointed to FastGA, sassy and
    mm2-plus. sassy does approximate pattern search, not genome alignment, so it was not considered.
    mm2-plus is not installed. FastGA is built at ~/bin/src/FASTGA and is being benchmarked.
