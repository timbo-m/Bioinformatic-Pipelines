# famscan: Find members of a protein family

_I want every member of a protein family in my genome/proteome, checked by structure_

- Stages: load -> family_search -> domains -> primary -> features -> motifs -> compare_refs -> filter -> structure
- Parameters: `{'rules': '{"length":[50,3000]}', 'evalue': 1e-05, 'top_n': 5, 'structure': 'auto', 'input': 'data/proteome_938293.faa', 'query': 'data/PF02826.hmm', 'hmm_db': 'data/mini_domain_library.hmm'}`
- Run time: 1.2 s; final candidates: 7

## Input

`data/proteome_938293.faa` read as **protein**: 2100 proteins.

## Family search

7 proteins match `data/PF02826.hmm` at E <= 1e-05.

## Domain architecture

7/7 proteins have >= 1 domain.

| domain | count |
|---|---|
| 2-Hacid_dh_C | 8 |
| KR | 2 |

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

## Motifs

6 motif matches.

| motif | count |
|---|---|
| Rossmann NAD(P)-binding GxGxxG | 4 |
| ATP/GTP-binding P-loop (Walker A) | 2 |

**De novo motifs shared by the set** (x = any residue):

| motif | seq_frac | enrichment | length |
|---|---|---|---|
| NxARG | 0.57 | 16 | 5 |
| LxEE | 0.71 | 10 | 4 |
| YExL | 0.71 | 10 | 4 |
| DLxK | 0.86 | 6 | 4 |
| LGxT | 0.71 | 6.7 | 4 |
| ExLE | 0.71 | 6.7 | 4 |
| NxAR | 0.57 | 8 | 4 |
| AxKK | 0.57 | 8 | 4 |
| TAxG | 0.71 | 5 | 4 |
| LxNA | 0.57 | 5.3 | 4 |

## Similarity to references

Compared with 1 reference sequence(s).

## Filters

`{'length': [50, 3000]}`: 7/7 pass.

## Tertiary structure

0/5 structures obtained. Missing ones need network access to AlphaFold DB / ESMFold, or fold them with ColabFold and pass the folder as `pdb_dir`.

## Final candidates

| id | evalue | score | length | pI | instability | localisation_guess | n_tm | signal_peptide | motifs | domains | ref_identity | structure_verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 938293.PRJEB85.HG003685_328 | 5.6e-54 | 179 | 311 | 4.99 | 23.2 | cytoplasmic | 0 | False | ATP/GTP-binding P-loop (Walker A) | 2-Hacid_dh_C(100-279) | 45.1 | no structure |
| 938293.PRJEB85.HG003686_519 | 7.4e-53 | 176 | 322 | 4.71 | 36.4 | cytoplasmic | 0 | False |  | 2-Hacid_dh_C(112-284) | 45.9 | no structure |
| 938293.PRJEB85.HG003689_19 | 6.3e-52 | 172 | 327 | 4.98 | 22.8 | cytoplasmic | 0 | False | ATP/GTP-binding P-loop (Walker A); Rossmann NAD(P)-binding G | 2-Hacid_dh_C(110-295) | 48 | no structure |
| 938293.PRJEB85.HG003686_258 | 4.1e-42 | 140 | 328 | 5.2 | 37.3 | cytoplasmic | 0 | False |  | 2-Hacid_dh_C(113-288) | 39.9 | no structure |
| 938293.PRJEB85.HG003685_432 | 7.6e-09 | 32.1 | 368 | 5.16 | 35 | cytoplasmic | 0 | False | Rossmann NAD(P)-binding GxGxxG | 2-Hacid_dh_C(147-274) | 27.7 | no structure |
| 938293.PRJEB85.HG003691_27 | 3.0e-08 | 30.2 | 364 | 5.16 | 41 | cytoplasmic | 0 | False | Rossmann NAD(P)-binding GxGxxG | 2-Hacid_dh_C(151-256) | 27.7 |  |
| 938293.PRJEB85.HG003688_14 | 4.4e-07 | 26.4 | 468 | 5.82 | 41.1 | cytoplasmic | 0 | False | Rossmann NAD(P)-binding GxGxxG | 2-Hacid_dh_C(118-198) / KR(139-184) / 2-Hacid_dh_C(259-365)  | 32.7 |  |

## Files

- `results.tsv`
- `final_candidates.faa`
- `domains.tsv`
- `secondary_structure.txt`
- `motif_hits.tsv`
- `denovo_motifs.tsv`
