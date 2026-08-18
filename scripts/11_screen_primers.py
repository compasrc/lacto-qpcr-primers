#!/usr/bin/env python3
"""Screen primer pairs against all 14 genomes.

A BLAST hit counts as a plausible priming event only if it covers nearly the
full primer with few mismatches. Each pair is then scored on whether both
primers hit all genomes of the target species and nothing else.
"""

import pandas as pd
from pathlib import Path

MIN_COV = 0.90      # fraction of primer length aligned
MAX_MM = 2          # mismatches tolerated within that alignment

COLS = ["qseqid", "sseqid", "pident", "length", "mismatch",
        "gapopen", "qstart", "qend", "evalue"]

hits = pd.read_csv("results/specificity/primer_hits.tsv", sep="\t", names=COLS)
c2g = pd.read_csv("work/contig2genome.tsv", sep="\t",
                  names=["sseqid", "genome"])
prim = pd.read_csv("results/primers/primer_candidates.csv")

# Primer lengths, keyed by the FASTA id used in the BLAST query
plen = {}
for _, r in prim.iterrows():
    tag = r.candidate.replace("Lactobacillus_", "L")
    plen[f"{tag}__p{r.pair}__L"] = len(r.left_seq)
    plen[f"{tag}__p{r.pair}__R"] = len(r.right_seq)

hits["plen"] = hits.qseqid.map(plen)
hits = hits.dropna(subset=["plen"])
hits["cov"] = hits["length"] / hits["plen"]

real = hits[(hits["cov"] >= MIN_COV) & (hits.mismatch <= MAX_MM)]
real = real.merge(c2g, on="sseqid", how="left")
real["sp_hit"] = real.genome.str.split("__").str[0]

rows = []
for _, r in prim.iterrows():
    tag = r.candidate.replace("Lactobacillus_", "L")
    target = r.candidate.split("__")[0]
    tgt_short = target

    rec = {"candidate": r.candidate, "pair": r.pair,
           "product_size": r.product_size}
    for side in ("L", "R"):
        q = f"{tag}__p{r.pair}__{side}"
        h = real[real.qseqid == q]
        on = h[h.sp_hit == tgt_short].genome.nunique()
        off = h[h.sp_hit != tgt_short]
        rec[f"{side}_on_target_genomes"] = on
        rec[f"{side}_off_target_genomes"] = off.genome.nunique()
        rec[f"{side}_off_target_species"] = ",".join(sorted(off.sp_hit.dropna().unique()))
    rows.append(rec)

res = pd.DataFrame(rows)

n_expected = {sp: n for sp, n in
              pd.read_csv("config/strains.csv").groupby("species").size().items()}
res["expected_genomes"] = res.candidate.str.split("__").str[0].map(n_expected)

res["pass"] = (
    (res.L_on_target_genomes == res.expected_genomes) &
    (res.R_on_target_genomes == res.expected_genomes) &
    (res.L_off_target_genomes == 0) &
    (res.R_off_target_genomes == 0)
)

res.to_csv("results/primers/primer_screen.csv", index=False)

print(f"{len(res)} pairs screened, {res['pass'].sum()} pass\n")
summary = res.groupby(res.candidate.str.split("__").str[0])["pass"].agg(["sum", "size"])
print(summary.to_string())
