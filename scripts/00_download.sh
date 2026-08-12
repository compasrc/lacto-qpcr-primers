#!/usr/bin/env bash
set -euo pipefail

mkdir -p work/download

datasets download genome accession \
    --inputfile config/accessions.txt \
    --include genome,seq-report \
    --filename work/download/genomes.zip

unzip -o work/download/genomes.zip -d work/download/genomes

cd work/download/genomes && md5sum -c md5sum.txt && cd -

echo "Downloaded $(find work/download/genomes -name '*_genomic.fna' | wc -l) genomes"
