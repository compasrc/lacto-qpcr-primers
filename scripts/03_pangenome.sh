#!/usr/bin/env bash
set -euo pipefail

OUTDIR="results/panaroo"
THREADS=$(nproc)

mkdir -p "$OUTDIR"

mapfile -t gffs < <(find data/annotations -name "*.gff" | sort)
echo "Found ${#gffs[@]} GFF files"

panaroo \
    -i "${gffs[@]}" \
    -o "$OUTDIR" \
    --clean-mode moderate \
    --threads "$THREADS" \
    -a core \
    --core_threshold 0.95

echo "Done. Outputs in $OUTDIR"
