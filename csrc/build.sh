#!/bin/bash
# Build alnfilter against a FASTGA source tree: bash build.sh /path/to/FASTGA  -> ../bin/alnfilter
set -euo pipefail
SRC=${1:?usage: build.sh /path/to/FASTGA}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$HERE/../bin"
gcc -O3 -Wall -fno-strict-aliasing -I"$SRC" -o "$HERE/../bin/alnfilter" "$HERE/alnfilter.c" \
    "$SRC/GDB.c" "$SRC/alncode.c" "$SRC/gene_core.c" "$SRC/ONElib.c" -lpthread -lm -lz
echo "built $HERE/../bin/alnfilter"
