"""Fixed internal constants, one place, each with its reason.

Nothing here is tuned to a dataset. Values that can be learned from the input (emission
densities, breakpoint rate, window size, resolution) are learned in model.py / windows.py.
"""

# --- sites -----------------------------------------------------------------------------------
# A diagnostic site needs >= 2 covered references that disagree, and the polyploid must match at
# least one of them: only then does the base say which reference it resembles. (Polyploid-private
# bases carry no origin information but pile up in fast-evolving regions.)
MIN_REFS_COVERED = 2

# Pseudocount per site pattern and state (Jeffreys-like: no log(0), rare patterns not inflated).
PSEUDOCOUNT = 0.5

# --- alignment prefilter (before ALNchain) ----------------------------------------------------------
# Raw FastGA alignments are kept if they lie in a collinear run (same query/target/strand, diagonal
# band of BAND bp, consecutive alignments <= GAP bp apart on the query) with >= MIN_BP aligned query
# bases. BAND and GAP equal ALNchain's default max gap (10 kb); MIN_BP = 5 kb because ALNchain's
# default minimum chain score (10,000) counts both genomes' extents, i.e. ~5 kb of alignment. The
# prefilter is thus a superset of ALNchain's output and only removes isolated repeat hits.
PREFILTER_BAND_BP = 10_000
PREFILTER_GAP_BP = 10_000
PREFILTER_MIN_BP = 5_000

# --- windows ---------------------------------------------------------------------------------------
# Per-base comparisons are summed in blocks of this size; windows are built from whole blocks.
BIN_BP = 1000

# Information per window: enough compared bases c that native-origin and other-origin sequence are
# expected to differ by sqrt(WINDOW_Z2) ~ 3.2 Poisson SDs in mismatch count against the native
# reference: c = WINDOW_Z2 (d1 + dh) / (dh - d1)^2, with d1 the typical divergence to the closest
# reference and dh the typical divergence between the polyploid's homoeologues (an exchanged segment
# sits ~dh from its native reference). Set by statistics, not by any genome (v0.4.8; before:
# ~10 mismatches to the closest reference, which fails when the references are the parents).
WINDOW_Z2 = 10

# A window that would need more than this multiple of the median window span to collect its
# information is not formed: the region is too sparsely aligned to be evaluated at the genome's
# own resolution (reported in nocall.bed).
WINDOW_MAX_FACTOR = 10

# Compared bases needed in a 10-kb block before its divergence enters the typical-divergence
# estimate (~10 mismatches at 1% divergence).
WINDOW_MIN_COMPARED = 1000

# --- emission densities (KDE) ----------------------------------------------------------------------
# Grid cells per dimension (dimension = min(#references - 1, #subgenomes - 1) <= 3).
KDE_GRID = {1: 400, 2: 100, 3: 40}
# Grid range: 0.1-99.9% quantiles of all windows, padded by 1 log unit (e-fold divergence ratio).
KDE_RANGE_Q = 0.001
KDE_PAD = 1.0


# --- HMM ---------------------------------------------------------------------------------------
# Bounds for the fitted breakpoint rate (per bp): 1e-9 = one switch per Gb, 1e-5 = one per
# 100 kb. Outside these the model is either flat or fitting noise.
RHO_MIN, RHO_MAX = 1e-9, 1e-5
RHO_INIT = 1e-7

# EM stopping: relative log-likelihood change and iteration cap; density re-estimation rounds
# after reassigning windows to their called origin (removes HE contamination of native densities).
EM_TOL = 1e-4
EM_MAX_ITER = 30
PROFILE_ROUNDS = 3

# Multi-scale proposals: extra segmentations with window evidence divided by sqrt(tau) and tau
# (integrated autocorrelation time of native-window evidence), when that factor is >= MIN_TEMPER
# (below it, tempering changes little). Raw segmentation shatters long exchanges whose windows are
# serially correlated (Saccharum chr02S, 2026-10-08); all scales are scored identically.
MIN_TEMPER = 2.0

# --- calls -------------------------------------------------------------------------------------
# PASS threshold: QUAL 20 = posterior 0.99 (conventional evidence level).
QUAL_PASS = 20.0

# Donor fit (v0.4.9): a call must look like sequence of the donor's origin against every reference.
# Null = runs of the same number of consecutive native windows of the donor subgenome (serial
# correlation kept); per reference, both tails; Bonferroni over references x 2. Conventional
# evidence level; 4999 draws keep power with up to ~25 references (min p 2e-4 x 2K < 0.01).
DONOR_FIT_P = 0.01
DONOR_FIT_DRAWS = 4999
QUAL_CAP = 99.0

# Breakpoint credible interval mass.
BREAKPOINT_CI = 0.95

# --- evidence checks ----------------------------------------------------------------------------
# Synteny check of an HE against the closest reference: the alignment covering the segment must
# land on the same reference chromosome as the nearest flanking alignments, and at a distance
# compatible with collinearity: |dt| <= SLACK + FOLD * |dq|. FOLD 3 allows a 3-fold local
# expansion/contraction (TE turnover); SLACK 1 Mb absorbs small local rearrangements. Positions are
# aligned-base-weighted medians over the segment and over flank windows of max(segment, FLANK_MIN)
# bp (FLANK_MIN = 100 kb: several genes' worth of collinear sequence even in large genomes), capped
# at FLANK_MAX (farther sequence cannot vouch for local synteny).
SYNTENY_SLACK_BP = 1_000_000
SYNTENY_FOLD = 3
SYNTENY_FLANK_MIN_BP = 100_000
SYNTENY_FLANK_MAX_BP = 5_000_000

# Collinear chain for projections: pieces whose diagonal offset lies within max(1 Mb, half the
# projected interval) of the weighted mode (indel drift across a segment stays well inside this;
# repeat hits elsewhere on the chromosome fall outside).
CHAIN_TOL_MIN_BP = 1_000_000

# A homoeologue is assigned to an HE segment only if the self alignment covers at least this
# fraction of it (otherwise 'unpaired').
HOMOEOLOG_MIN_COVER = 0.25
# Identity null: blocks (and segments) need homoeologue-compared bases over at least this fraction
# of their span for their divergence to be compared.
IDENTITY_MIN_COVER = 0.25
