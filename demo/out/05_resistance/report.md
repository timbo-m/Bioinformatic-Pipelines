# famscan: Screen for antibiotic resistance

_I want to know which antibiotic-resistance genes a bacterium carries_

- Stages: load -> resistance -> primary -> features
- Parameters: `{'strict_identity': 80, 'input': 'data/bacillus_contig.fna', 'resistance_db': 'data/standin_screen_db.faa'}`
- Run time: 2.6 s; final candidates: 6

## Input

`data/bacillus_contig.fna` read as **genome**: 287 proteins.

## Antibiotic resistance

6 proteins match `data/standin_screen_db.faa`. Classes follow CARD RGI logic: perfect (identical), strict (>= 80 % identity, >= 60 % coverage), loose (weaker homologs, possibly novel variants).

| protein | amr_hit | amr_desc | amr_identity | amr_coverage | amr_evalue | amr_class |
|---|---|---|---|---|---|---|
| 1390.SAMEA104415756.OFHT01000022_121 | standin_3 | demo entry (mutated copy of 1390.SAMEA104415756.OFHT01000022 | 90.5 | 100 | 0 | strict |
| 1390.SAMEA104415756.OFHT01000022_181 | standin_4 | demo entry (mutated copy of 1390.SAMEA104415756.OFHT01000022 | 94 | 100 | 2.0e-204 | strict |
| 1390.SAMEA104415756.OFHT01000022_241 | standin_5 | demo entry (mutated copy of 1390.SAMEA104415756.OFHT01000022 | 89.7 | 100 | 5.3e-154 | strict |
| 1390.SAMEA104415756.OFHT01000022_61 | standin_2 | demo entry (mutated copy of 1390.SAMEA104415756.OFHT01000022 | 96.6 | 100 | 2.8e-41 | strict |
| 1390.SAMEA104415756.OFHT01000022_1 | standin_1 | demo entry (mutated copy of 1390.SAMEA104415756.OFHT01000022 | 94.7 | 86.4 | 8.9e-10 | strict |
| 1390.SAMEA104415756.OFHT01000022_182 | standin_4 | demo entry (mutated copy of 1390.SAMEA104415756.OFHT01000022 | 27.7 | 25 | 2.8e-06 | loose |

## Primary structure

| index | length | mw_kda | pI | instability | gravy |
|---|---|---|---|---|---|
| mean | 310 | 34.8 | 6.65 | 43.9 | -0.28 |
| min | 22 | 2.04 | 4.36 | 27 | -0.74 |
| max | 803 | 90.9 | 10.3 | 72.2 | 0.1 |

## Sequence features (fast predictors)

Signal peptide: 0, >=1 TM helix: 0, coiled coil: 0, >30 % disordered: 0.

| localisation | count |
|---|---|
| cytoplasmic | 6 |

Per-residue H/E/C strings: `secondary_structure.txt`. These are quick heuristics; confirm with SignalP 6, DeepTMHMM, IUPred3 and PSIPRED.

## Final candidates

| id | length | pI | instability | localisation_guess | n_tm | signal_peptide | amr_class | amr_desc |
|---|---|---|---|---|---|---|---|---|
| 1390.SAMEA104415756.OFHT01000022_1 | 22 | 5.27 | 27 | cytoplasmic | 0 | False | strict | demo entry (mutated copy of 1390.SAMEA104415756.OFHT01000022 |
| 1390.SAMEA104415756.OFHT01000022_61 | 59 | 10.3 | 72.2 | cytoplasmic | 0 | False | strict | demo entry (mutated copy of 1390.SAMEA104415756.OFHT01000022 |
| 1390.SAMEA104415756.OFHT01000022_121 | 803 | 8.83 | 39.6 | cytoplasmic | 0 | False | strict | demo entry (mutated copy of 1390.SAMEA104415756.OFHT01000022 |
| 1390.SAMEA104415756.OFHT01000022_181 | 332 | 4.99 | 38.9 | cytoplasmic | 0 | False | strict | demo entry (mutated copy of 1390.SAMEA104415756.OFHT01000022 |
| 1390.SAMEA104415756.OFHT01000022_182 | 378 | 4.36 | 44.4 | cytoplasmic | 0 | False | loose | demo entry (mutated copy of 1390.SAMEA104415756.OFHT01000022 |
| 1390.SAMEA104415756.OFHT01000022_241 | 263 | 6.12 | 41.3 | cytoplasmic | 0 | False | strict | demo entry (mutated copy of 1390.SAMEA104415756.OFHT01000022 |

## Files

- `results.tsv`
- `final_candidates.faa`
- `predicted_proteins.faa`
- `secondary_structure.txt`
