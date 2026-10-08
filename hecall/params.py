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

# --- windows ---------------------------------------------------------------------------------------
# Per-base comparisons are summed in blocks of this size; windows are built from whole blocks.
BIN_BP = 1000

# Information per window: enough compared bases that ~10 mismatches to the closest reference are
# expected at this genome's typical divergence (Poisson CV ~0.3, so a window's log-divergence
# ratio has SE ~0.3-0.45, small against the separation between origins). Set by statistics, not by
# any genome.
WINDOW_MISMATCHES = 10

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

# Longest lag (windows) for the autocorrelation-time estimate; Geyer's sequence normally stops much
# earlier (first non-positive pair of lags).
TAU_MAX_LAG = 2000

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

# --- calls -------------------------------------------------------------------------------------
# PASS threshold: QUAL 20 = posterior 0.99 (conventional evidence level).
QUAL_PASS = 20.0
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

# A homoeologue is assigned to an HE segment only if the self alignment covers at least this
# fraction of it (otherwise 'unpaired').
HOMOEOLOG_MIN_COVER = 0.25
# Identity null: blocks (and segments) need homoeologue-compared bases over at least this fraction
# of their span for their divergence to be compared.
IDENTITY_MIN_COVER = 0.25
