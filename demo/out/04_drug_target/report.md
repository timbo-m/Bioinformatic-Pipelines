# famscan: Find drug targets

_I want proteins in a pathogen that humans don't have (potential drug or vaccine targets)_

- Stages: load -> subtract_host -> domains -> primary -> features -> structure -> rank
- Parameters: `{'top_n': 10, 'structure': 'none', 'input': 'data/proteome_938293.faa', 'host': 'data/bacillus_proteins.faa', 'hmm_db': 'data/mini_domain_library.hmm'}`
- Run time: 22.9 s; final candidates: 1755

> Secreted/membrane proteins in the output are vaccine candidates; essential, cytoplasmic enzymes with good folds are classic drug targets.

## Input

`data/proteome_938293.faa` read as **protein**: 2100 proteins.

## Host subtraction

345 proteins have a homolog in `data/bacillus_proteins.faa` and are dropped; 1755 pathogen-specific proteins remain.

## Domain architecture

14/1755 proteins have >= 1 domain.

| domain | count |
|---|---|
| 2-Hacid_dh_C | 11 |
| Thuricin_rSAM_RRE | 1 |
| Ranthipeptide_rSAM_RRE | 1 |
| LuxC | 1 |
| Thurincin_rSAM_RRE | 1 |
| PqqD_RRE | 1 |

## Primary structure

| index | length | mw_kda | pI | instability | gravy |
|---|---|---|---|---|---|
| mean | 311 | 35 | 6.41 | 33.9 | -0.24 |
| min | 32 | 2.45 | 4.05 | -1.5 | -1.97 |
| max | 4.56e+03 | 509 | 12 | 105 | 1.72 |

## Sequence features (fast predictors)

Signal peptide: 223, >=1 TM helix: 287, coiled coil: 84, >30 % disordered: 43.

| localisation | count |
|---|---|
| cytoplasmic | 1292 |
| membrane | 264 |
| secreted/periplasmic | 139 |
| lipoprotein | 60 |

Per-residue H/E/C strings: `secondary_structure.txt`. These are quick heuristics; confirm with SignalP 6, DeepTMHMM, IUPred3 and PSIPRED.

## Tertiary structure

Skipped (structure=none).

## Ranking

Higher is better: stable (instability < 40), predicted soluble, ordered, no TM helices, length 100-600 aa, confident fold.

| id | rank_score | length | pI | instability | wh_soluble | disorder_frac | n_tm | localisation_guess | host_homolog | domains |
|---|---|---|---|---|---|---|---|---|---|---|
| 938293.PRJEB85.HG003685_422 | 2.48 | 188 | 4.56 | 24.7 | True | 0 | 0 | cytoplasmic | False |  |
| 938293.PRJEB85.HG003690_14 | 2.48 | 100 | 4.45 | 31.6 | True | 0 | 0 | cytoplasmic | False |  |
| 938293.PRJEB85.HG003690_6 | 2.48 | 116 | 10.9 | 10.1 | True | 0 | 0 | cytoplasmic | False |  |
| 938293.PRJEB85.HG003687_5 | 2.48 | 149 | 4.59 | 22.4 | True | 0 | 0 | cytoplasmic | False |  |
| 938293.PRJEB85.HG003686_179 | 2.48 | 117 | 11 | 28.4 | True | 0 | 0 | cytoplasmic | False |  |
| 938293.PRJEB85.HG003685_52 | 2.47 | 249 | 4.4 | 31 | True | 0 | 0 | cytoplasmic | False |  |
| 938293.PRJEB85.HG003690_135 | 2.47 | 118 | 4.12 | 36.5 | True | 0 | 0 | cytoplasmic | False |  |
| 938293.PRJEB85.HG003685_321 | 2.46 | 103 | 4.65 | 32.7 | True | 0.01 | 0 | cytoplasmic | False |  |
| 938293.PRJEB85.HG003686_465 | 2.45 | 162 | 4.18 | 35.8 | True | 0 | 0 | cytoplasmic | False |  |
| 938293.PRJEB85.HG003685_259 | 2.45 | 285 | 4.26 | 39.9 | True | 0 | 0 | cytoplasmic | False |  |
| 938293.PRJEB85.HG003686_7 | 2.45 | 168 | 4.44 | 30.9 | True | 0.03 | 0 | cytoplasmic | False |  |
| 938293.PRJEB85.HG003686_231 | 2.45 | 203 | 4.92 | 23.6 | True | 0 | 0 | cytoplasmic | False |  |
| 938293.PRJEB85.HG003686_555 | 2.44 | 233 | 4.37 | 37.8 | True | 0 | 0 | cytoplasmic | False |  |
| 938293.PRJEB85.HG003686_811 | 2.44 | 102 | 4.23 | 35.3 | True | 0 | 0 | cytoplasmic | False |  |
| 938293.PRJEB85.HG003687_118 | 2.44 | 143 | 4.05 | 34.2 | True | 0.03 | 0 | cytoplasmic | False |  |
| 938293.PRJEB85.HG003685_334 | 2.44 | 273 | 4.3 | 24.4 | True | 0 | 0 | cytoplasmic | False |  |
| 938293.PRJEB85.HG003685_70 | 2.43 | 238 | 4.21 | 33.2 | True | 0 | 0 | cytoplasmic | False |  |
| 938293.PRJEB85.HG003687_107 | 2.42 | 153 | 4.45 | 35.8 | True | 0 | 0 | cytoplasmic | False |  |
| 938293.PRJEB85.HG003686_608 | 2.41 | 183 | 4.71 | 38.3 | True | 0.04 | 0 | cytoplasmic | False |  |
| 938293.PRJEB85.HG003686_927 | 2.41 | 282 | 10.2 | 29.5 | True | 0 | 0 | cytoplasmic | False |  |
| 938293.PRJEB85.HG003686_817 | 2.41 | 103 | 4.46 | 12.3 | True | 0 | 0 | cytoplasmic | False |  |
| 938293.PRJEB85.HG003690_205 | 2.41 | 151 | 4.45 | 37 | True | 0.01 | 0 | cytoplasmic | False |  |
| 938293.PRJEB85.HG003686_724 | 2.41 | 173 | 4.65 | 26.6 | True | 0 | 0 | cytoplasmic | False |  |
| 938293.PRJEB85.HG003686_166 | 2.4 | 290 | 6.13 | 23.8 | True | 0 | 0 | cytoplasmic | False |  |
| 938293.PRJEB85.HG003685_246 | 2.4 | 599 | 4.6 | 36.4 | True | 0 | 0 | cytoplasmic | False |  |

## Final candidates

| id | rank_score | length | pI | instability | localisation_guess | n_tm | signal_peptide | domains |
|---|---|---|---|---|---|---|---|---|
| 938293.PRJEB85.HG003685_422 | 2.48 | 188 | 4.56 | 24.7 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003690_14 | 2.48 | 100 | 4.45 | 31.6 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003690_6 | 2.48 | 116 | 10.9 | 10.1 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003687_5 | 2.48 | 149 | 4.59 | 22.4 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003686_179 | 2.48 | 117 | 11 | 28.4 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003685_52 | 2.47 | 249 | 4.4 | 31 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003690_135 | 2.47 | 118 | 4.12 | 36.5 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003685_321 | 2.46 | 103 | 4.65 | 32.7 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003686_465 | 2.45 | 162 | 4.18 | 35.8 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003685_259 | 2.45 | 285 | 4.26 | 39.9 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003686_7 | 2.45 | 168 | 4.44 | 30.9 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003686_231 | 2.45 | 203 | 4.92 | 23.6 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003686_555 | 2.44 | 233 | 4.37 | 37.8 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003686_811 | 2.44 | 102 | 4.23 | 35.3 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003687_118 | 2.44 | 143 | 4.05 | 34.2 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003685_334 | 2.44 | 273 | 4.3 | 24.4 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003685_70 | 2.43 | 238 | 4.21 | 33.2 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003687_107 | 2.42 | 153 | 4.45 | 35.8 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003686_608 | 2.41 | 183 | 4.71 | 38.3 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003686_927 | 2.41 | 282 | 10.2 | 29.5 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003686_817 | 2.41 | 103 | 4.46 | 12.3 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003690_205 | 2.41 | 151 | 4.45 | 37 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003686_724 | 2.41 | 173 | 4.65 | 26.6 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003686_166 | 2.4 | 290 | 6.13 | 23.8 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003685_246 | 2.4 | 599 | 4.6 | 36.4 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003687_199 | 2.4 | 127 | 4.93 | 35.8 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003686_223 | 2.4 | 176 | 4.51 | 28.2 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003687_44 | 2.39 | 147 | 4.49 | 26.7 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003686_284 | 2.39 | 593 | 4.74 | 39.6 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003687_156 | 2.39 | 273 | 4.5 | 26.4 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003685_270 | 2.39 | 209 | 5.01 | 24.8 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003686_371 | 2.38 | 166 | 4.76 | 31.4 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003690_56 | 2.38 | 175 | 4.57 | 22.7 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003685_19 | 2.38 | 226 | 4.67 | 21.5 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003686_528 | 2.38 | 131 | 4.43 | 27 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003685_184 | 2.38 | 157 | 4.59 | 35.4 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003685_292 | 2.37 | 176 | 4.99 | 26.4 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003686_688 | 2.37 | 256 | 4.59 | 31.3 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003686_487 | 2.37 | 408 | 4.65 | 35 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003686_420 | 2.36 | 429 | 4.84 | 36.4 | cytoplasmic | 0 | False |  |

## Files

- `results.tsv`
- `final_candidates.faa`
- `domains.tsv`
- `secondary_structure.txt`
