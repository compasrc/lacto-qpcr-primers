#!/usr/bin/env python3
"""Extract each shortlisted candidate gene's DNA sequence from every strain
of its species, writing one unaligned FASTA per gene."""

import pandas as pd
from pathlib import Path

SHORTLIST = Path("config/shortlist.tsv")
OUTDIR = Path("work/candidate_seqs")
META = ("Gene", "Non-unique Gene name", "Annotation")

OUTDIR.mkdir(parents=True, exist_ok=True)
short = pd.read_csv(SHORTLIST, sep="\t")

for sp, grp in short.groupby("species"):
    base = Path("results/panaroo_per_species") / sp
    pa = pd.read_csv(base / "gene_presence_absence.csv", low_memory=False)
    gd = pd.read_csv(base / "gene_data.csv", low_memory=False)

    seq_of = dict(zip(gd.annotation_id, gd.dna_sequence))
    strain_of = dict(zip(gd.annotation_id, gd.gff_file))
    n_strains = len([c for c in pa.columns if c not in META])

    for gene in grp.gene:
        hit = pa[pa.Gene == gene]
        if hit.empty:
            print(f"  ! {sp} {gene}: not found in matrix")
            continue

        row = hit.iloc[0]
        tags = [str(v).strip() for v in row[list(pa.columns[3:])]
                if str(v).strip() not in ("", "nan")]

        records = []
        for tag in tags:
            for t in tag.split(";"):
                seq = seq_of.get(t)
                if isinstance(seq, str) and seq:
                    strain = Path(str(strain_of.get(t, t))).stem
                    records.append((strain, t, seq))

        out = OUTDIR / f"{sp}__{gene}.fna"
        with open(out, "w") as fh:
            for strain, tag, seq in records:
                fh.write(f">{strain}|{tag}\n{seq}\n")

        lens = [len(s) for _, _, s in records]
        flag = "" if len(records) == n_strains else "  <-- INCOMPLETE"
        print(f"  {sp:<26} {gene:<12} {len(records)}/{n_strains} strains  "
              f"{min(lens)}-{max(lens)} bp{flag}" if lens else
              f"  {sp:<26} {gene:<12} NO SEQUENCES")

print(f"\nWrote FASTAs to {OUTDIR}/")
