#!/bin/bash
# Build hecall's C helpers against a FASTGA source tree: bash build.sh /path/to/FASTGA
#   ../bin/alnfilter  copies the .1aln records kept by the prefilter
#   ../bin/ALNchain   upstream ALNchain with alnchain_ties.patch: a three-way quickselect partition, so
#                     thousands of alignments ending at one position (e.g. at a contig edge in a repeat
#                     array) no longer make KD-tree construction quadratic (Digitaria: > 80 min -> 1 s).
#                     If the patch does not apply to this FASTGA version, it is skipped with a warning
#                     and hecall uses the upstream ALNchain.
set -euo pipefail
SRC=${1:?usage: build.sh /path/to/FASTGA}
HERE=$(cd "$(dirname "$0")" && pwd)
BIN="$HERE/../bin"
LIBS=("$SRC/GDB.c" "$SRC/alncode.c" "$SRC/gene_core.c" "$SRC/ONElib.c" -lpthread -lm -lz)
mkdir -p "$BIN"
gcc -O3 -Wall -fno-strict-aliasing -I"$SRC" -o "$BIN/alnfilter" "$HERE/alnfilter.c" "${LIBS[@]}"
echo "built $BIN/alnfilter"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
cp "$SRC/ALNchain.c" "$TMP/"
if patch -s -d "$TMP" -p0 < "$HERE/alnchain_ties.patch"; then
  gcc -O3 -fno-strict-aliasing -Wno-unused-result -I"$SRC" -o "$BIN/ALNchain" "$TMP/ALNchain.c" "${LIBS[@]}"
  echo "built $BIN/ALNchain (tie-safe)"
else
  rm -f "$BIN/ALNchain"
  echo "warning: alnchain_ties.patch does not apply to $SRC/ALNchain.c; hecall will use the upstream ALNchain" >&2
fi
