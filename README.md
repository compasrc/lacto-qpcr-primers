# Strain-Specific qPCR Primer Design for Vaginal *Lactobacillus* Species

Pipeline for designing species-specific qPCR primers (~100 bp amplicons) to
quantify relative abundance of *Lactobacillus* species in co-culture, using
comparative pan-genome analysis of publicly available genome assemblies.

## Target species

- *L. crispatus*
- *L. jensenii*
- *L. mulieris*
- *L. paragasseri*
- *L. gasseri*

## Planned pipeline

1. Genome retrieval from NCBI (`config/genomes.tsv`)
2. Annotation — Prokka
3. Pan-genome analysis — Panaroo
4. Species-specific gene identification (presence/absence filtering)
5. Within-species alignment — MAFFT
6. Primer design — Primer3
7. In silico specificity screening — BLAST+
8. Wet-lab validation — melt curve, standard curve, efficiency

## Design targets

| Parameter | Target |
|---|---|
| Amplicon size | 80–140 bp |
| Primer Tm | 59–62°C |
| ΔTm within pair | ≤1°C |
| GC content | 40–60% |
| Primer length | 18–24 nt |

## Status

Project initialized. No candidates selected yet.

## Reproducing

Environment: `environment.yml`. Documentation: `docs/`.
