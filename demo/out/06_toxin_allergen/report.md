# famscan: Screen for toxins and allergens

_I need to check proteins for toxin or allergen risk (safety, food, biotech)_

- Stages: load -> primary -> features -> toxin_allergen -> motifs
- Parameters: `{'input': 'data/proteome_938293.faa', 'toxin_db': 'data/standin_screen_db.faa', 'allergen_db': 'data/standin_screen_db.faa'}`
- Run time: 37.9 s; final candidates: 7

## Input

`data/proteome_938293.faa` read as **protein**: 2100 proteins.

## Primary structure

| index | length | mw_kda | pI | instability | gravy |
|---|---|---|---|---|---|
| mean | 324 | 36.5 | 6.32 | 33.9 | -0.24 |
| min | 32 | 2.45 | 4.05 | -1.5 | -1.97 |
| max | 4.56e+03 | 509 | 12 | 105 | 1.72 |

## Sequence features (fast predictors)

Signal peptide: 237, >=1 TM helix: 336, coiled coil: 117, >30 % disordered: 44.

| localisation | count |
|---|---|
| cytoplasmic | 1579 |
| membrane | 309 |
| secreted/periplasmic | 151 |
| lipoprotein | 61 |

Per-residue H/E/C strings: `secondary_structure.txt`. These are quick heuristics; confirm with SignalP 6, DeepTMHMM, IUPred3 and PSIPRED.

## Toxins and allergens

**Toxins:** 4 proteins match `data/standin_screen_db.faa`.

**Allergens:** 2 proteins flagged by Codex criteria (> 35 % identity over >= 80 aa, or an identical 8-residue stretch) against `data/standin_screen_db.faa`.

**Toxin-like by sequence:** 3 small (<= 150 aa), secreted, Cys-rich (>= 4 Cys) proteins - the profile of many peptide toxins.

**7 proteins flagged in total** (listed under Final candidates).

## Motifs

7 motif matches.

| motif | count |
|---|---|
| ATP/GTP-binding P-loop (Walker A) | 3 |
| Walker B (DEAD/DEAH-like) | 3 |
| Peroxisome targeting PTS1 | 1 |

**De novo motifs shared by the set** (x = any residue):

| motif | seq_frac | enrichment | length |
|---|---|---|---|
| NxGE | 0.57 | 8 | 4 |
| EIxE | 0.71 | 5 | 4 |
| GExI | 0.57 | 5.3 | 4 |
| GxDI | 0.57 | 5.3 | 4 |
| IxEI | 0.71 | 4 | 4 |
| IDxI | 0.71 | 4 | 4 |
| IExK | 0.57 | 4 | 4 |
| IDxL | 0.57 | 4 | 4 |
| RSxI | 0.57 | 4 | 4 |
| SxLD | 0.57 | 4 | 4 |

## Final candidates

| id | length | pI | instability | localisation_guess | n_tm | signal_peptide | motifs | allergen_flag |
|---|---|---|---|---|---|---|---|---|
| 938293.PRJEB85.HG003690_132 | 144 | 5.03 | 39.9 | secreted/periplasmic | 0 | True |  | False |
| 938293.PRJEB85.HG003691_69 | 84 | 9.36 | 38.3 | secreted/periplasmic | 1 | True | Peroxisome targeting PTS1 | False |
| 938293.PRJEB85.HG003686_82 | 59 | 9.8 | 86.8 | cytoplasmic | 0 | False |  | False |
| 938293.PRJEB85.HG003686_397 | 103 | 4.9 | 43 | secreted/periplasmic | 0 | True |  | False |
| 938293.PRJEB85.HG003686_776 | 1162 | 5.11 | 33.2 | cytoplasmic | 0 | False | ATP/GTP-binding P-loop (Walker A); Walker B (DEAD/DEAH-like) | False |
| 938293.PRJEB85.HG003686_785 | 756 | 8.83 | 44.2 | cytoplasmic | 0 | False | ATP/GTP-binding P-loop (Walker A); Walker B (DEAD/DEAH-like) | True |
| 938293.PRJEB85.HG003687_64 | 655 | 4.95 | 48.9 | cytoplasmic | 0 | False | ATP/GTP-binding P-loop (Walker A); Walker B (DEAD/DEAH-like) | True |

## Files

- `results.tsv`
- `final_candidates.faa`
- `secondary_structure.txt`
- `motif_hits.tsv`
- `denovo_motifs.tsv`
