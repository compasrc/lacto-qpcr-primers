#!/usr/bin/env bash
set -euo pipefail

INDIR="work/candidate_seqs"
OUTDIR="work/candidate_aln"

mkdir -p "$OUTDIR"

for f in "$INDIR"/*.fna; do
    base=$(basename "$f" .fna)
    out="$OUTDIR/$base.aln"

    if [[ -s "$out" ]]; then
        echo "[skip] $base"
        continue
    fi

    mafft --auto --quiet --thread "$(nproc)" "$f" > "$out"
    printf "%-45s %s cols\n" "$base" \
        "$(awk '/^>/{if(n)exit} !/^>/{c+=length($0)} END{print c}' "$out")"
done

echo
echo "Aligned $(ls -1 "$OUTDIR"/*.aln | wc -l) candidates"
