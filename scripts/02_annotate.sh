#!/usr/bin/env bash
set -euo pipefail

INDIR="data/genomes"
OUTDIR="data/annotations"
THREADS=$(nproc)

mkdir -p "$OUTDIR"

for fna in "$INDIR"/*.fna; do
    base=$(basename "$fna" .fna)
    outdir="$OUTDIR/$base"

    if [[ -s "$outdir/$base.gff" ]]; then
        echo "[skip] $base"
        continue
    fi

    genus=$(echo "$base" | cut -d_ -f1)
    species=$(echo "$base" | cut -d_ -f2)
    strain=$(echo "$base" | awk -F'__' '{print $2}')

    echo "[run ] $base"
    prokka \
        --outdir "$outdir" \
        --prefix "$base" \
        --genus "$genus" \
        --species "$species" \
        --strain "$strain" \
        --kingdom Bacteria \
        --cpus "$THREADS" \
        --force \
        --quiet \
        "$fna"
done

echo
echo "Annotated: $(ls -1d "$OUTDIR"/*/ 2>/dev/null | wc -l) genomes"
