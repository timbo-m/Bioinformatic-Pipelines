# famscan: Discover new enzymes

_I'm looking for new enzymes of a known type (e.g. novel lipases, PETases)_

- Stages: load -> family_search -> domains -> primary -> features -> motifs -> compare_refs -> filter -> cluster -> msa -> tree -> structure -> rank
- Parameters: `{'rules': '{"length":[50,3000]}', 'cluster_identity': 90, 'evalue': 1e-05, 'top_n': 10, 'structure': 'none', 'input': 'data/bacillus_contig.fna', 'query': 'data/Thioesterase.hmm', 'hmm_db': 'data/mini_domain_library.hmm', 'motifs': 'thioesterase_nucleophile=G-x-S-x-G'}`
- Run time: 0.3 s; final candidates: 1

> Use `motifs` for catalytic signatures, e.g. lipase=G-x-S-x-G.

## Input

`data/bacillus_contig.fna` read as **genome**: 287 proteins.

## Family search

1 proteins match `data/Thioesterase.hmm` at E <= 1e-05.

## Domain architecture

1/1 proteins have >= 1 domain.

| domain | count |
|---|---|
| Thioesterase | 1 |

## Primary structure

| index | length | mw_kda | pI | instability | gravy |
|---|---|---|---|---|---|
| mean | 2.48e+03 | 276 | 5.35 | 43.9 | -0.43 |
| min | 2.48e+03 | 276 | 5.35 | 43.9 | -0.43 |
| max | 2.48e+03 | 276 | 5.35 | 43.9 | -0.43 |

## Sequence features (fast predictors)

Signal peptide: 0, >=1 TM helix: 0, coiled coil: 0, >30 % disordered: 0.

| localisation | count |
|---|---|
| cytoplasmic | 1 |

Per-residue H/E/C strings: `secondary_structure.txt`. These are quick heuristics; confirm with SignalP 6, DeepTMHMM, IUPred3 and PSIPRED.

## Motifs

6 motif matches.

| motif | count |
|---|---|
| Ser hydrolase nucleophile elbow GxSxG | 3 |
| thioesterase_nucleophile | 3 |

## Similarity to references

Compared with 1 reference sequence(s).

## Filters

`{'length': [50, 3000]}`: 1/1 pass.

## Clustering

1 sequences -> 1 clusters at 90 % identity; continuing with one representative per cluster.

## Alignment

Skipped: fewer than 3 sequences.

## Tertiary structure

Skipped (structure=none).

## Ranking

Higher is better: stable (instability < 40), predicted soluble, ordered, no TM helices, length 100-600 aa, confident fold.

| id | rank_score | length | pI | instability | wh_soluble | disorder_frac | n_tm | localisation_guess | domains |
|---|---|---|---|---|---|---|---|---|---|
| 1390.SAMEA104415756.OFHT01000022_277 | 0.23 | 2485 | 5.35 | 43.9 | True | 0.08 | 0 | cytoplasmic | Thioesterase(2182-2450) |

## Final candidates

| id | rank_score | evalue | score | length | pI | instability | localisation_guess | n_tm | signal_peptide | motifs | domains | ref_identity |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1390.SAMEA104415756.OFHT01000022_277 | 0.23 | 4.7e-44 | 145 | 2485 | 5.35 | 43.9 | cytoplasmic | 0 | False | Ser hydrolase nucleophile elbow GxSxG x3; thioesterase_nucle | Thioesterase(2182-2450) | 49.5 |

## Files

- `results.tsv`
- `final_candidates.faa`
- `predicted_proteins.faa`
- `domains.tsv`
- `secondary_structure.txt`
- `motif_hits.tsv`
