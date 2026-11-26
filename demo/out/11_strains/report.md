# famscan: Compare strains

_I have several strains/isolates and want shared and unique genes_

- Stages: load_many -> strain_compare -> domains -> features
- Parameters: `{'proteomes': ['data/bacillus_contig.fna', 'data/bacillus_contig_variant.fna'], 'hmm_db': 'data/mini_domain_library.hmm'}`
- Run time: 31.7 s; final candidates: 9

## Inputs

| proteome | kind | proteins |
|---|---|---|
| bacillus_contig | genome | 287 |
| bacillus_contig_variant | genome | 299 |

## Strain comparison

Reference strain `bacillus_contig`; orthologs by reciprocal best hit (E <= 1e-10).

| category | genes |
|---|---|
| core (all strains) | 278 |
| accessory | 0 |
| only in bacillus_contig | 9 |
| only in bacillus_contig_variant (no hit in bacillus_contig) | 1 |

## Domain architecture

1/9 proteins have >= 1 domain.

| domain | count |
|---|---|
| KR | 1 |

## Sequence features (fast predictors)

Signal peptide: 0, >=1 TM helix: 2, coiled coil: 0, >30 % disordered: 0.

| localisation | count |
|---|---|
| cytoplasmic | 7 |
| membrane | 2 |

Per-residue H/E/C strings: `secondary_structure.txt`. These are quick heuristics; confirm with SignalP 6, DeepTMHMM, IUPred3 and PSIPRED.

## Final candidates

| id | length | localisation_guess | n_tm | signal_peptide | domains |
|---|---|---|---|---|---|
| bacillus_contig/1390.SAMEA104415756.OFHT01000022_1 | 22 | cytoplasmic | 0 | False |  |
| bacillus_contig/1390.SAMEA104415756.OFHT01000022_5 | 1590 | cytoplasmic | 0 | False | KR(275-388) |
| bacillus_contig/1390.SAMEA104415756.OFHT01000022_63 | 160 | cytoplasmic | 0 | False |  |
| bacillus_contig/1390.SAMEA104415756.OFHT01000022_64 | 306 | cytoplasmic | 0 | False |  |
| bacillus_contig/1390.SAMEA104415756.OFHT01000022_65 | 558 | cytoplasmic | 0 | False |  |
| bacillus_contig/1390.SAMEA104415756.OFHT01000022_66 | 154 | cytoplasmic | 0 | False |  |
| bacillus_contig/1390.SAMEA104415756.OFHT01000022_67 | 311 | cytoplasmic | 0 | False |  |
| bacillus_contig/1390.SAMEA104415756.OFHT01000022_68 | 117 | membrane | 1 | False |  |
| bacillus_contig/1390.SAMEA104415756.OFHT01000022_69 | 714 | membrane | 1 | False |  |

## Files

- `results.tsv`
- `final_candidates.faa`
- `presence_absence.tsv`
- `strain_unique_genes.json`
- `domains.tsv`
- `secondary_structure.txt`
