---
name: famscan
description: Guided bioinformatics pipelines. Use for protein family searches, enzyme discovery, genome annotation, drug targets, resistance, toxin and allergen screening, evolution, unknown proteins, lab candidate choice, engineering variants, strain comparison and any BLAST direction.
---

# famscan

1. Work out which use case fits the user's request (run `famscan list` to see them).
   If it's unclear, ask: "What are you trying to do?" and offer the 12 options.
2. Collect the required inputs for that workflow. Inputs can be a protein FASTA,
   a genome (genes are called automatically) or CDS sequences.
3. Choose a search engine by job size:
   - many sequences or whole proteomes: the built-in engine (HMMER), or `blast+`/`diamond`
     against a downloaded database.
   - a few queries against NCBI nr/swissprot: `engine=ncbi`. Never automate the
     NCBI web form; the official API allows about 1 request every 10 s.
4. Run `famscan workflow <name> --out <dir> --set key=value ...` and read `<dir>/report.md`.
5. Summarise for the user: counts at each stage, the final candidates and why they
   rank highly, and anything stopped or skipped (for example, structures that need network access).
6. Point out the sequence-only heuristics (signal peptide, TM, disorder, coiled coil,
   secondary structure) and suggest confirming final picks with SignalP 6, DeepTMHMM,
   IUPred3, PSIPRED or the structure-derived calls.
7. For more than ~50 structures, suggest local ColabFold or ESMFold and `--set pdb_dir=...`.
