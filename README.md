# Species-Specific qPCR Primer Design for Vaginal *Lactobacillus*

A reproducible pipeline for designing species-specific qPCR primers to
quantify the relative abundance of five *Lactobacillus* species in co-culture,
using comparative pan-genome analysis of 14 clinical isolate genomes.

## Problem

Five vaginal *Lactobacillus* species are grown in co-culture and their relative
abundance compared by qPCR. This requires one primer pair per species that
amplifies every strain of its target species and none of the other genomes in
the experiment. Two species pairs make this non-trivial: *L. gasseri* /
*L. paragasseri* were a single species until 2018, and *L. jensenii* /
*L. mulieris* until 2020.

Because the co-culture is a closed system, specificity has a definitive answer
— every primer is screened against all 14 genomes directly rather than against
a reference database.

## Final designs

| Species | Target gene | Forward | Reverse | Amplicon |
|---|---|---|---|---|
| *L. crispatus* | `group_623` | CGCCGCAGTTGTTTCATTCA | TTCATGGCTTGCATTTGGGC | 102 bp |
| *L. gasseri* | `group_97` | ACGCGACCTGGAAGTTGATC | TCCTTGTTGCCTTGGATCTCC | 130 bp |
| *L. jensenii* | `group_214` | CGTTCTAGCTCTGGTGGTGG | TGCAGCAGTATTTGGACCATCT | 102 bp |
| *L. mulieris* | `guaD` | GTGGCGCTTGATAGTCCTCT | GCAGTAGTTGTCCCATGAGCT | 140 bp |
| *L. paragasseri* | `asp2` | TTGCTGATACCGAGGAAAGCT | ACGCATTCCTGTTTGACCCA | 122 bp |

All Tm 59.4–60.4 °C, ΔTm within pair < 0.7 °C. Full candidate set with
alternates in `results/primers/`.

## Input genomes

14 UMB clinical isolates (GenBank assemblies, accessions in
`config/strains.csv`): *L. crispatus* ×4, *L. gasseri* ×2, *L. jensenii* ×3,
*L. mulieris* ×2, *L. paragasseri* ×3. NCBI organism names were verified
against strain labels before analysis; all matched.

## Pipeline

| Script | Step |
|---|---|
| `00_download.sh` | Retrieve assemblies via NCBI Datasets, verify MD5 |
| `01_stage_genomes.sh` | Flatten to `species__strain__accession.fna` |
| `02_annotate.sh` | Prokka annotation (separate conda env) |
| `03b_pangenome_per_species.sh` | Panaroo, run **per species** (see notes) |
| `05_specificity_screen.sh` | Gene-level BLAST vs. out-group genomes |
| `06_rank_candidates.py` | Rank core genes by worst-case out-group similarity |
| `07_extract_candidates.py` | Pull shortlisted genes from every strain |
| `08_align_candidates.sh` | MAFFT alignment within species |
| `09_find_conserved_windows.py` | Locate conserved primer sites |
| `10_design_primers.py` | Primer3 on masked consensus (qPCR constraints) |
| `11_screen_primers.py` | Primer-level BLAST vs. all 14 genomes |

`03_pangenome.sh` (all 14 genomes at once) is retained for the diagnostic
result documented in `docs/notes.md`; it is not part of the working path.

## Design constraints

Amplicon 80–140 bp · primer 18–24 nt · Tm 59–62 °C · ΔTm ≤ 1 °C ·
GC 30–65 % (widened from the usual 40–60 % because these genomes are ~35 % GC) ·
poly-X ≤ 3 · dimer and hairpin penalties tightened for SYBR Green.

Primers are designed on a masked consensus in which any alignment column that
varies across strains, or contains a gap, is replaced with `N`. Primer3 is set
to reject any primer containing `N`, so every design sits in a region conserved
across all strains of its species.

## Results

| Species | Core genes | Zero out-group hits | Pairs designed | Pairs passing screen |
|---|---|---|---|---|
| *L. crispatus* | 1608 | 319 | 25 | 20 |
| *L. gasseri* | 1576 | 16 | 15 | 11 |
| *L. jensenii* | 1337 | 19 | 23 | 18 |
| *L. mulieris* | 1281 | 35 | 25 | 18 |
| *L. paragasseri* | 1535 | 22 | 10 | 2 |

A pair passes only if both primers align near-full-length to every genome of
the target species and to no other genome in the set. 69 of 98 pairs passed.

The candidate counts track phylogeny: *L. crispatus* is the most distinct of
the five, while *L. gasseri* and *L. paragasseri* — separated in 2018 — yield
the fewest specific genes.

## Reproducing

```bash
conda env create -f environment.yml            # analysis tools
conda env create -f environment-prokka.yml     # annotation (isolated)

conda activate lacto-qpcr && bash scripts/00_download.sh
bash scripts/01_stage_genomes.sh
conda activate prokka   && bash scripts/02_annotate.sh
conda activate lacto-qpcr
bash scripts/03b_pangenome_per_species.sh
bash scripts/05_specificity_screen.sh
python3 scripts/06_rank_candidates.py
python3 scripts/07_extract_candidates.py
bash scripts/08_align_candidates.sh
python3 scripts/09_find_conserved_windows.py
python3 scripts/10_design_primers.py
python3 scripts/11_screen_primers.py
```

Genome FASTAs, annotations, and BLAST output are not version-controlled; they
are regenerated from `config/strains.csv` and the scripts above.

## Limitations

- *L. gasseri* and *L. mulieris* are represented by two genomes each, so
  "present in all strains" is a weak constraint for those species. Adequate
  for this closed co-culture; not sufficient to claim species-wide specificity.
- *L. paragasseri* has only one independent passing design (`asp2`). The
  alternative gene, `group_266`, was rejected — see `docs/notes.md`.
- Two assemblies are metagenome-assembled genomes (UMB8929B, UMB7783B) and may
  carry binning artifacts.
- Several assemblies are heavily fragmented (up to 400 contigs), which can
  cause genes spanning contig breaks to be called absent.
- Designs are validated *in silico* only. Wet-lab confirmation requires melt
  curve analysis (single peak), a standard curve across 5–6 log dilutions with
  90–110 % efficiency and R² > 0.99, and cross-reactivity testing against
  genomic DNA from each non-target species.

## Notes

`docs/notes.md` records two methodological findings worth reading before
modifying the pipeline: why Panaroo is run per species rather than across all
14 genomes, and why `blastn` requires an explicit `-task` for cross-species
comparison.
