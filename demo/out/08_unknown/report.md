# famscan: Characterise unknown proteins

_I have proteins with no or weak annotation and want to know what they might be_

- Stages: load -> name_by_homology -> domains -> primary -> features -> motifs -> structure
- Parameters: `{'top_n': 10, 'structure': 'none', 'foldseek_web': False, 'input': 'data/unknown_proteins.faa', 'hmm_db': 'data/mini_domain_library.hmm'}`
- Run time: 0.2 s; final candidates: 6

> Structure search (Foldseek) often finds function when sequence search fails.

## Input

`data/unknown_proteins.faa` read as **protein**: 6 proteins.

## Naming by homology

Skipped: no `annotation_db` (e.g. uniprot_sprot.fasta) given.

## Domain architecture

0/6 proteins have >= 1 domain.

_none_

## Primary structure

| index | length | mw_kda | pI | instability | gravy |
|---|---|---|---|---|---|
| mean | 362 | 40.2 | 5.69 | 36.3 | -0.17 |
| min | 287 | 32.2 | 4.69 | 16.5 | -0.87 |
| max | 511 | 56 | 8.82 | 46 | 0.54 |

## Sequence features (fast predictors)

Signal peptide: 1, >=1 TM helix: 1, coiled coil: 2, >30 % disordered: 1.

| localisation | count |
|---|---|
| cytoplasmic | 4 |
| membrane | 1 |
| secreted/periplasmic | 1 |

Per-residue H/E/C strings: `secondary_structure.txt`. These are quick heuristics; confirm with SignalP 6, DeepTMHMM, IUPred3 and PSIPRED.

## Motifs

4 motif matches.

| motif | count |
|---|---|
| Ser hydrolase nucleophile elbow GxSxG | 1 |
| Leucine zipper | 1 |
| ATP/GTP-binding P-loop (Walker A) | 1 |
| Walker B (DEAD/DEAH-like) | 1 |

**De novo motifs shared by the set** (x = any residue):

| motif | seq_frac | enrichment | length |
|---|---|---|---|
| KMLE | 0.5 | 12 | 4 |
| VxGGV | 0.5 | 12 | 5 |
| VxKAL | 0.5 | 12 | 5 |
| LFxLA | 0.5 | 12 | 5 |
| SxTG | 0.5 | 12 | 4 |
| AFxK | 0.5 | 12 | 4 |
| VxAF | 0.5 | 12 | 4 |
| ELxR | 0.5 | 12 | 4 |
| AxLA | 0.5 | 12 | 4 |
| NKxT | 0.5 | 12 | 4 |

## Tertiary structure

Skipped (structure=none).

## Final candidates

| id | length | pI | instability | localisation_guess | n_tm | signal_peptide | motifs |
|---|---|---|---|---|---|---|---|
| 938293.PRJEB85.HG003689_7 | 302 | 5.14 | 40.4 | membrane | 1 | False | Ser hydrolase nucleophile elbow GxSxG |
| 938293.PRJEB85.HG003690_76 | 287 | 4.69 | 46 | secreted/periplasmic | 0 | True |  |
| 938293.PRJEB85.HG003690_176 | 320 | 8.82 | 36.8 | cytoplasmic | 0 | False | Leucine zipper |
| 938293.PRJEB85.HG003691_10 | 511 | 5.1 | 41.1 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003691_110 | 314 | 5.11 | 16.5 | cytoplasmic | 0 | False |  |
| 938293.PRJEB85.HG003685_100 | 439 | 5.28 | 37.1 | cytoplasmic | 0 | False | ATP/GTP-binding P-loop (Walker A); Walker B (DEAD/DEAH-like) |

## Files

- `results.tsv`
- `final_candidates.faa`
- `domains.tsv`
- `secondary_structure.txt`
- `motif_hits.tsv`
- `denovo_motifs.tsv`
