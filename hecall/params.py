"""Fixed internal constants, one place, each with its reason.

Nothing here is tuned to a dataset. Values that can be learned from the input (emission
profiles, breakpoint rate, resolution) are learned in model.py.
"""

# --- sites -----------------------------------------------------------------------------------
# Minimum references covering a base for it to be informative: with one reference a mismatch
# cannot be attributed to either lineage.
MIN_REFS_COVERED = 2

# Pseudocount added to every pattern count per state (Jeffreys-like: avoids log(0) for patterns
# never seen in one state without inflating rare ones).
PSEUDOCOUNT = 0.5

# A pattern seen fewer times than this over all states is pooled into one 'rare' category:
# below ~30 counts a frequency estimate has >18% relative error (1/sqrt(30)).
MIN_PATTERN_COUNT = 30

# --- HMM ---------------------------------------------------------------------------------------
# Bin size (bp) for pooling sites. Below the shortest detectable HE in practice (>= ~10 sites);
# breakpoints are refined to single-site resolution afterwards, so bins do not limit resolution.
BIN_BP = 1000

# Bounds for the fitted breakpoint rate (per bp): 1e-9 = one switch per Gb, 1e-5 = one per
# 100 kb. Outside these the model is either flat or fitting noise.
RHO_MIN, RHO_MAX = 1e-9, 1e-5
RHO_INIT = 1e-7

# EM stopping: relative log-likelihood change, and iteration caps.
EM_TOL = 1e-4
EM_MAX_ITER = 30
PROFILE_ROUNDS = 3  # profile re-estimation rounds after removing called HEs (see design 2.2)

# --- calls -------------------------------------------------------------------------------------
# PASS threshold: QUAL 20 = posterior 0.99 (conventional evidence level).
QUAL_PASS = 20.0
QUAL_CAP = 99.0

# Breakpoint credible interval mass.
BREAKPOINT_CI = 0.95

# --- reporting -----------------------------------------------------------------------------------
# Window for identity / site-count tracks (bins.tsv.gz) and figures.
WINDOW_BP = 10_000

# --- evidence checks ----------------------------------------------------------------------------
# Synteny check of an HE against the closest reference: the alignment covering the segment must
# land on the same reference chromosome as the nearest flanking alignments, and at a distance
# compatible with collinearity: |dt| <= SLACK + FOLD * |dq|. FOLD 3 allows a 3-fold local
# expansion/contraction (TE turnover); SLACK 1 Mb absorbs small local rearrangements. Flanks
# farther than FLANK_MAX are not used (too far to vouch for local synteny).
SYNTENY_SLACK_BP = 1_000_000
SYNTENY_FOLD = 3
SYNTENY_FLANK_MAX_BP = 5_000_000

# A homoeologue is assigned to an HE segment only if the self alignment covers at least this
# fraction of it (otherwise 'unpaired').
HOMOEOLOG_MIN_COVER = 0.25
