#!/usr/bin/env python3
"""Design qPCR primers in conserved regions of each candidate alignment.

Non-conserved alignment columns are masked to N so Primer3 avoids them.
"""

from Bio import AlignIO
from pathlib import Path
import subprocess, glob, re

OUTDIR = Path("results/primers")
WORK = Path("work/primer3")
OUTDIR.mkdir(parents=True, exist_ok=True)
WORK.mkdir(parents=True, exist_ok=True)

SETTINGS = """PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_NUM_RETURN=5
PRIMER_PRODUCT_SIZE_RANGE=80-140
PRIMER_OPT_SIZE=20
PRIMER_MIN_SIZE=18
PRIMER_MAX_SIZE=24
PRIMER_OPT_TM=60.0
PRIMER_MIN_TM=59.0
PRIMER_MAX_TM=62.0
PRIMER_PAIR_MAX_DIFF_TM=1.0
PRIMER_MIN_GC=30.0
PRIMER_MAX_GC=65.0
PRIMER_MAX_POLY_X=3
PRIMER_MAX_SELF_ANY_TH=45.0
PRIMER_MAX_SELF_END_TH=35.0
PRIMER_PAIR_MAX_COMPL_ANY_TH=45.0
PRIMER_PAIR_MAX_COMPL_END_TH=35.0
PRIMER_MAX_HAIRPIN_TH=24.0
PRIMER_MAX_NS_ACCEPTED=0
PRIMER_EXPLAIN_FLAG=1
"""

def masked_consensus(path):
    aln = AlignIO.read(path, "fasta")
    out = []
    for i in range(aln.get_alignment_length()):
        col = aln[:, i].upper()
        if "-" in col or len(set(col)) != 1:
            out.append("N")
        else:
            out.append(col[0])
    return "".join(out)

rows = []

for f in sorted(glob.glob("work/candidate_aln/*.aln")):
    name = Path(f).stem
    seq = masked_consensus(f)

    inp = WORK / f"{name}.p3in"
    inp.write_text(f"SEQUENCE_ID={name}\nSEQUENCE_TEMPLATE={seq}\n{SETTINGS}=\n")

    res = subprocess.run(["primer3_core", str(inp)],
                         capture_output=True, text=True)
    out = res.stdout
    (WORK / f"{name}.p3out").write_text(out)

    d = dict(line.split("=", 1) for line in out.splitlines() if "=" in line)
    n = int(d.get("PRIMER_PAIR_NUM_RETURNED", 0))

    for i in range(n):
        rows.append({
            "candidate": name,
            "pair": i,
            "left_seq": d.get(f"PRIMER_LEFT_{i}_SEQUENCE", ""),
            "right_seq": d.get(f"PRIMER_RIGHT_{i}_SEQUENCE", ""),
            "left_tm": d.get(f"PRIMER_LEFT_{i}_TM", ""),
            "right_tm": d.get(f"PRIMER_RIGHT_{i}_TM", ""),
            "left_gc": d.get(f"PRIMER_LEFT_{i}_GC_PERCENT", ""),
            "right_gc": d.get(f"PRIMER_RIGHT_{i}_GC_PERCENT", ""),
            "product_size": d.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE", ""),
            "penalty": d.get(f"PRIMER_PAIR_{i}_PENALTY", ""),
        })

    print(f"{name:<45} {n} pairs")

import pandas as pd
df = pd.DataFrame(rows)
df.to_csv(OUTDIR / "primer_candidates.csv", index=False)
print(f"\n{len(df)} primer pairs written to {OUTDIR}/primer_candidates.csv")
