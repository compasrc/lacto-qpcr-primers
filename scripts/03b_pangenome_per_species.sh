#!/usr/bin/env bash
set -euo pipefail

THREADS=$(nproc)

for sp in Lactobacillus_crispatus Lactobacillus_gasseri Lactobacillus_jensenii \
          Lactobacillus_mulieris Lactobacillus_paragasseri; do

    outdir="results/panaroo_per_species/$sp"
    if [[ -s "$outdir/gene_presence_absence.csv" ]]; then
        echo "[skip] $sp"
        continue
    fi

    mapfile -t gffs < <(find data/annotations -name "${sp}__*.gff" | sort)
    echo "[run ] $sp (${#gffs[@]} genomes)"

    mkdir -p "$outdir"
    panaroo -i "${gffs[@]}" -o "$outdir" \
        --clean-mode moderate \
        --threads "$THREADS" \
        --core_threshold 1.0
done

echo
echo "Done."
