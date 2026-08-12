#!/usr/bin/env bash
set -euo pipefail

SRC="work/download/genomes/ncbi_dataset/data"
DEST="data/genomes"

mkdir -p "$DEST"

tail -n +2 config/strains.csv | while IFS=',' read -r species strain acc; do
    [[ -z "$acc" ]] && continue

    src=$(find "$SRC/$acc" -name "*_genomic.fna" -type f 2>/dev/null | head -1)
    if [[ -z "$src" ]]; then
        echo "  ! missing: $acc ($species $strain)" >&2
        continue
    fi

    dest="$DEST/${species}__${strain}__${acc}.fna"
    cp "$src" "$dest"
    echo "  staged $(basename "$dest")"
done

echo
echo "$(ls -1 "$DEST"/*.fna | wc -l) genomes staged"
