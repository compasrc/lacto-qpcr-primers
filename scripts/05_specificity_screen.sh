#!/usr/bin/env bash
set -euo pipefail

SPECIES=(Lactobacillus_crispatus Lactobacillus_gasseri Lactobacillus_jensenii
         Lactobacillus_mulieris Lactobacillus_paragasseri)
THREADS=$(nproc)
OUT="results/specificity"

mkdir -p "$OUT" work/blastdb

for sp in "${SPECIES[@]}"; do
    echo "=== $sp"

    outgroup="work/blastdb/${sp}__outgroup.fna"
    find data/genomes -name "*.fna" ! -name "${sp}__*" -print0 \
        | xargs -0 cat > "$outgroup"

    n_out=$(find data/genomes -name "*.fna" ! -name "${sp}__*" | wc -l)
    echo "  outgroup: $n_out genomes"

    makeblastdb -in "$outgroup" -dbtype nucl \
        -out "work/blastdb/${sp}__outgroup" -logfile /dev/null

    ref="results/panaroo_per_species/$sp/pan_genome_reference.fa"
    blastn \
        -query "$ref" \
        -db "work/blastdb/${sp}__outgroup" \
        -outfmt "6 qseqid qlen sseqid pident length qstart qend evalue bitscore" \
        -task blastn \
        -evalue 1e-3 \
        -num_threads "$THREADS" \
        -max_target_seqs 5000 \
        > "$OUT/${sp}__hits.tsv"

    echo "  queries: $(grep -c '^>' "$ref")"
    echo "  hits:    $(wc -l < "$OUT/${sp}__hits.tsv")"
done
