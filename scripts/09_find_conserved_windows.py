#!/usr/bin/env python3
"""Report conserved regions in each candidate alignment.

A column is conserved if all strains agree and none has a gap. Reports the
longest gap-free fully-conserved blocks and whether the alignment can support
a qPCR amplicon (two 20bp conserved windows separated by 80-140 bp).
"""

from Bio import AlignIO
from pathlib import Path
import glob

PRIMER_LEN = 20
AMP_MIN, AMP_MAX = 80, 140

rows = []

for f in sorted(glob.glob("work/candidate_aln/*.aln")):
    aln = AlignIO.read(f, "fasta")
    n, width = len(aln), aln.get_alignment_length()

    cons = []
    for i in range(width):
        col = aln[:, i]
        cons.append(len(set(col.upper())) == 1 and "-" not in col)

    # Longest run of conserved columns
    best = run = 0
    for c in cons:
        run = run + 1 if c else 0
        best = max(best, run)

    # Can we place two 20bp conserved windows 80-140 apart?
    ok = [i for i in range(width - PRIMER_LEN + 1)
          if all(cons[i:i + PRIMER_LEN])]
    okset = set(ok)
    pairs = sum(1 for i in ok
                for amp in range(AMP_MIN, AMP_MAX + 1)
                if (i + amp - PRIMER_LEN) in okset)

    pct = 100 * sum(cons) / width
    rows.append((Path(f).stem, n, width, round(pct, 1), best, len(ok), pairs))

hdr = f"{'candidate':<45} {'n':>2} {'width':>6} {'%cons':>6} {'longest':>8} {'sites':>6} {'pairs':>6}"
print(hdr)
print("-" * len(hdr))
for r in sorted(rows, key=lambda x: -x[6]):
    print(f"{r[0]:<45} {r[1]:>2} {r[2]:>6} {r[3]:>6} {r[4]:>8} {r[5]:>6} {r[6]:>6}")
