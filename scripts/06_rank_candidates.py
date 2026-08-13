#!/usr/bin/env python3
"""Rank per-species core genes by cross-species specificity.

For each gene cluster in a species' pan-genome, summarise the worst-case
BLAST hit against the outgroup genomes: highest percent identity, longest
alignment, and the fraction of the query covered.
"""

import pandas as pd
from pathlib import Path

SPECIES = [
    "Lactobacillus_crispatus", "Lactobacillus_gasseri", "Lactobacillus_jensenii",
    "Lactobacillus_mulieris", "Lactobacillus_paragasseri",
]
COLS = ["qseqid", "qlen", "sseqid", "pident", "length",
        "qstart", "qend", "evalue", "bitscore"]
OUTDIR = Path("results/candidates")
OUTDIR.mkdir(parents=True, exist_ok=True)

summary = []

for sp in SPECIES:
    pa = pd.read_csv(f"results/panaroo_per_species/{sp}/gene_presence_absence.csv",
                     low_memory=False)
    strain_cols = [c for c in pa.columns
                   if c not in ("Gene", "Non-unique Gene name", "Annotation")]

    # Core = present in every strain of this species
    present = pa[strain_cols].notna() & pa[strain_cols].map(
        lambda v: str(v).strip() != "")
    core = pa.loc[present.all(axis=1), ["Gene", "Non-unique Gene name", "Annotation"]]

    hits = pd.read_csv(f"results/specificity/{sp}__hits.tsv", sep="\t", names=COLS)

    # Worst-case outgroup similarity per query
    agg = hits.groupby("qseqid").agg(
        max_pident=("pident", "max"),
        max_aln_len=("length", "max"),
        n_hits=("pident", "size"),
        qlen=("qlen", "first"),
    ).reset_index()
    agg["max_cov"] = (agg["max_aln_len"] / agg["qlen"]).round(3)

    merged = core.merge(agg, left_on="Gene", right_on="qseqid", how="left")
    merged["n_hits"] = merged["n_hits"].fillna(0).astype(int)
    merged["max_pident"] = merged["max_pident"].fillna(0)
    merged["max_aln_len"] = merged["max_aln_len"].fillna(0).astype(int)
    merged["max_cov"] = merged["max_cov"].fillna(0)

    merged = merged.sort_values(["max_pident", "max_cov"]).drop(columns=["qseqid"])
    merged.to_csv(OUTDIR / f"{sp}__ranked.csv", index=False)

    clean = (merged["n_hits"] == 0).sum()
    weak = ((merged["max_pident"] < 85) & (merged["max_aln_len"] < 100)).sum()
    summary.append((sp, len(core), clean, weak))

print(f"{'Species':<28} {'core':>5} {'no hit':>7} {'weak only':>10}")
for row in summary:
    print(f"{row[0]:<28} {row[1]:>5} {row[2]:>7} {row[3]:>10}")
print(f"\nWrote ranked CSVs to {OUTDIR}/")
