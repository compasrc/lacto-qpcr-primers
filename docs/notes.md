
## BLAST task selection for cross-species screening

The initial specificity screen used blastn's default task (megablast, word
size 28), which is tuned for near-identical sequences. At the 80-90%
nucleotide identity separating these species, megablast frequently failed to
seed alignments, producing false "no hit" results for universally conserved
genes — dnaE, mfd, smc, uvrB, addA and pyrAB all appeared as crispatus-specific.
Switching to -task blastn (word size 11) is required for cross-species
comparison at this divergence. Zero-hit counts dropped from 1203 to 319 for
crispatus and the pronounced asymmetry between species resolved.
