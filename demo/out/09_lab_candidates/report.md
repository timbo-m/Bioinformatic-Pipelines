# famscan: Pick candidates for the lab

_I have candidate proteins and need to choose which to clone and express_

- Stages: load -> primary -> features -> filter -> rank -> structure -> rank
- Parameters: `{'length_range': '100-600', 'top_n': 10, 'structure': 'none', 'rules': '{"length":[50,3000]}', 'input': 'out/01_family/final_candidates.faa'}`
- Run time: 0.1 s; final candidates: 7

## Input

`out/01_family/final_candidates.faa` read as **protein**: 7 proteins.

## Primary structure

| index | length | mw_kda | pI | instability | gravy |
|---|---|---|---|---|---|
| mean | 355 | 39.5 | 5.15 | 33.8 | -0.24 |
| min | 311 | 34.2 | 4.71 | 22.8 | -0.5 |
| max | 468 | 50.8 | 5.82 | 41.1 | -0.05 |

## Sequence features (fast predictors)

Signal peptide: 0, >=1 TM helix: 0, coiled coil: 0, >30 % disordered: 0.

| localisation | count |
|---|---|
| cytoplasmic | 7 |

Per-residue H/E/C strings: `secondary_structure.txt`. These are quick heuristics; confirm with SignalP 6, DeepTMHMM, IUPred3 and PSIPRED.

## Filters

`{'length': [50, 3000]}`: 7/7 pass.

## Ranking

Higher is better: stable (instability < 40), predicted soluble, ordered, no TM helices, length 100-600 aa, confident fold.

| id | rank_score | length | pI | instability | wh_soluble | disorder_frac | n_tm | localisation_guess |
|---|---|---|---|---|---|---|---|---|
| 938293.PRJEB85.HG003686_519 | 1.76 | 322 | 4.71 | 36.4 | True | 0 | 0 | cytoplasmic |
| 938293.PRJEB85.HG003689_19 | 1.7 | 327 | 4.98 | 22.8 | True | 0 | 0 | cytoplasmic |
| 938293.PRJEB85.HG003685_432 | 1.57 | 368 | 5.16 | 35 | True | 0.01 | 0 | cytoplasmic |
| 938293.PRJEB85.HG003685_328 | 1.54 | 311 | 4.99 | 23.2 | True | 0 | 0 | cytoplasmic |
| 938293.PRJEB85.HG003686_258 | 1.51 | 328 | 5.2 | 37.3 | True | 0 | 0 | cytoplasmic |
| 938293.PRJEB85.HG003691_27 | 0.49 | 364 | 5.16 | 41 | True | 0.03 | 0 | cytoplasmic |
| 938293.PRJEB85.HG003688_14 | -0.75 | 468 | 5.82 | 41.1 | False | 0.05 | 0 | cytoplasmic |

## Tertiary structure

Skipped (structure=none).

## Ranking

Higher is better: stable (instability < 40), predicted soluble, ordered, no TM helices, length 100-600 aa, confident fold.

| id | rank_score | length | pI | instability | wh_soluble | disorder_frac | n_tm | localisation_guess |
|---|---|---|---|---|---|---|---|---|
| 938293.PRJEB85.HG003686_519 | 1.76 | 322 | 4.71 | 36.4 | True | 0 | 0 | cytoplasmic |
| 938293.PRJEB85.HG003689_19 | 1.7 | 327 | 4.98 | 22.8 | True | 0 | 0 | cytoplasmic |
| 938293.PRJEB85.HG003685_432 | 1.57 | 368 | 5.16 | 35 | True | 0.01 | 0 | cytoplasmic |
| 938293.PRJEB85.HG003685_328 | 1.54 | 311 | 4.99 | 23.2 | True | 0 | 0 | cytoplasmic |
| 938293.PRJEB85.HG003686_258 | 1.51 | 328 | 5.2 | 37.3 | True | 0 | 0 | cytoplasmic |
| 938293.PRJEB85.HG003691_27 | 0.49 | 364 | 5.16 | 41 | True | 0.03 | 0 | cytoplasmic |
| 938293.PRJEB85.HG003688_14 | -0.75 | 468 | 5.82 | 41.1 | False | 0.05 | 0 | cytoplasmic |

## Final candidates

| id | rank_score | length | pI | instability | localisation_guess | n_tm | signal_peptide |
|---|---|---|---|---|---|---|---|
| 938293.PRJEB85.HG003686_519 | 1.76 | 322 | 4.71 | 36.4 | cytoplasmic | 0 | False |
| 938293.PRJEB85.HG003689_19 | 1.7 | 327 | 4.98 | 22.8 | cytoplasmic | 0 | False |
| 938293.PRJEB85.HG003685_432 | 1.57 | 368 | 5.16 | 35 | cytoplasmic | 0 | False |
| 938293.PRJEB85.HG003685_328 | 1.54 | 311 | 4.99 | 23.2 | cytoplasmic | 0 | False |
| 938293.PRJEB85.HG003686_258 | 1.51 | 328 | 5.2 | 37.3 | cytoplasmic | 0 | False |
| 938293.PRJEB85.HG003691_27 | 0.49 | 364 | 5.16 | 41 | cytoplasmic | 0 | False |
| 938293.PRJEB85.HG003688_14 | -0.75 | 468 | 5.82 | 41.1 | cytoplasmic | 0 | False |

## Files

- `results.tsv`
- `final_candidates.faa`
- `secondary_structure.txt`
