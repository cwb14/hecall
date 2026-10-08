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
  - Python ≥ 3.9 with numpy and matplotlib

## How it works (one paragraph)

Every reference genome is aligned to the polyploid with FastGA, keeping only syntenic chains
(ALNchain). At each polyploid base that two or more references cover, hecall records which
references carry the same base, e.g. "matches the first reference only". The frequencies of these
patterns are learned on each subgenome of the polyploid itself. "What B-donor sequence looks like"
is therefore measured on the polyploid's own B subgenome. This is why stand-in or distant
references still work: what matters is that the references *tell the donors apart*, not that they
*are* the donors. A hidden Markov model walks along each chromosome and assigns every stretch to
the donor whose pattern it matches. The breakpoint rate is estimated from the data. Where the
origin differs from the chromosome's subgenome, that stretch is an HE.

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
| `resolution.tsv` | per subgenome pair: information per site, sites per kb, **smallest detectable HE** |
| `profiles.tsv` | the learned site-pattern profiles of each subgenome |
| `bins.tsv.gz` | 10-kb windows: identity to every reference and homoeologue, sites, posterior |
| `painting.pdf/.png` + `painting_legend.md` | chromosome painting with every PASS HE marked |

### Key columns of `he.tsv`

- **`qual`**: Phred-scaled posterior that the segment has the called donor origin. 20 = 99%.
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
- **`llr`, `n_sites`, `n_support`, `n_against`**: the raw evidence, independent of the HMM prior.
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
