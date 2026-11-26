# famscan: Annotate a new genome

_I have a newly sequenced genome and want to know what its proteins are_

- Stages: load -> name_by_homology -> domains -> primary -> features -> motifs
- Parameters: `{'input': 'data/bacillus_contig.fna', 'hmm_db': 'data/mini_domain_library.hmm', 'annotation_db': 'data/reference_proteins.faa'}`
- Run time: 10.2 s; final candidates: 287

> For a fully curated annotation also run Bakta or Prokka; eukaryotes need BRAKER first.

## Input

`data/bacillus_contig.fna` read as **genome**: 287 proteins.

## Naming by homology

20/287 proteins have a hit in `data/reference_proteins.faa`.

## Domain architecture

15/287 proteins have >= 1 domain.

| domain | count |
|---|---|
| KR | 21 |
| 2-Hacid_dh_C | 3 |
| Thioesterase | 1 |

## Primary structure

| index | length | mw_kda | pI | instability | gravy |
|---|---|---|---|---|---|
| mean | 424 | 47.1 | 6.5 | 36.8 | -0.21 |
| min | 22 | 2.04 | 4.05 | 6.9 | -1.09 |
| max | 5.43e+03 | 600 | 11.5 | 72.2 | 1.16 |

## Sequence features (fast predictors)

Signal peptide: 26, >=1 TM helix: 36, coiled coil: 16, >30 % disordered: 4.

| localisation | count |
|---|---|
| cytoplasmic | 231 |
| membrane | 35 |
| secreted/periplasmic | 17 |
| lipoprotein | 4 |

Per-residue H/E/C strings: `secondary_structure.txt`. These are quick heuristics; confirm with SignalP 6, DeepTMHMM, IUPred3 and PSIPRED.

## Motifs

136 motif matches.

| motif | count |
|---|---|
| Ser hydrolase nucleophile elbow GxSxG | 48 |
| ATP/GTP-binding P-loop (Walker A) | 38 |
| Rossmann NAD(P)-binding GxGxxG | 15 |
| Leucine zipper | 13 |
| Cupin-like metal-binding HxH...E | 10 |
| Walker B (DEAD/DEAH-like) | 9 |
| Peroxisome targeting PTS1 | 3 |

**De novo motifs shared by the set** (x = any residue):

_none_

## Final candidates

| id | length | pI | instability | localisation_guess | n_tm | signal_peptide | motifs | domains | product |
|---|---|---|---|---|---|---|---|---|---|
| 1390.SAMEA104415756.OFHT01000022_1 | 22 | 5.27 | 27 | cytoplasmic | 0 | False |  |  |  |
| 1390.SAMEA104415756.OFHT01000022_2 | 214 | 9.24 | 20.1 | lipoprotein | 0 | True |  |  |  |
| 1390.SAMEA104415756.OFHT01000022_3 | 768 | 5.73 | 43 | cytoplasmic | 0 | False | Ser hydrolase nucleophile elbow GxSxG x2 |  | Narbonolide/10-deoxymethynolide synthase PikA1, modules 1 an |
| 1390.SAMEA104415756.OFHT01000022_4 | 4084 | 5.77 | 39.5 | cytoplasmic | 0 | False | ATP/GTP-binding P-loop (Walker A) x4; Walker B (DEAD/DEAH-li | KR(1160-1290) / KR(2680-2804) | Narbonolide/10-deoxymethynolide synthase PikA1, modules 1 an |
| 1390.SAMEA104415756.OFHT01000022_5 | 1590 | 5.56 | 41.7 | cytoplasmic | 0 | False | Leucine zipper; Ser hydrolase nucleophile elbow GxSxG |  | Narbonolide/10-deoxymethynolide synthase PikA1, modules 1 an |
| 1390.SAMEA104415756.OFHT01000022_6 | 2922 | 5.46 | 40.4 | cytoplasmic | 0 | False | Rossmann NAD(P)-binding GxGxxG | KR(296-449) / KR(1766-1898) | Narbonolide/10-deoxymethynolide synthase PikA1, modules 1 an |
| 1390.SAMEA104415756.OFHT01000022_7 | 2333 | 5.91 | 38.9 | cytoplasmic | 0 | False | Ser hydrolase nucleophile elbow GxSxG | KR(269-384) / KR(1463-1600) | Narbonolide/10-deoxymethynolide synthase PikA1, modules 1 an |
| 1390.SAMEA104415756.OFHT01000022_8 | 1903 | 5.51 | 37.5 | cytoplasmic | 0 | False | ATP/GTP-binding P-loop (Walker A); Ser hydrolase nucleophile | KR(752-907) | Narbonolide/10-deoxymethynolide synthase PikA1, modules 1 an |
| 1390.SAMEA104415756.OFHT01000022_9 | 2460 | 5.96 | 43.6 | cytoplasmic | 0 | False | Cupin-like metal-binding HxH...E | KR(273-462) / KR(1603-1717) | Narbonolide/10-deoxymethynolide synthase PikA1, modules 1 an |
| 1390.SAMEA104415756.OFHT01000022_10 | 1284 | 5.69 | 39.8 | cytoplasmic | 0 | False | Ser hydrolase nucleophile elbow GxSxG; Cupin-like metal-bind | KR(708-879) | Narbonolide/10-deoxymethynolide synthase PikA1, modules 1 an |
| 1390.SAMEA104415756.OFHT01000022_11 | 363 | 5.22 | 38.7 | cytoplasmic | 0 | False |  |  |  |
| 1390.SAMEA104415756.OFHT01000022_12 | 371 | 6.08 | 44.6 | cytoplasmic | 0 | False | Ser hydrolase nucleophile elbow GxSxG |  |  |
| 1390.SAMEA104415756.OFHT01000022_13 | 325 | 4.76 | 42.9 | cytoplasmic | 0 | False |  |  |  |
| 1390.SAMEA104415756.OFHT01000022_14 | 442 | 5.18 | 35 | cytoplasmic | 0 | False |  |  |  |
| 1390.SAMEA104415756.OFHT01000022_15 | 470 | 4.91 | 28.3 | cytoplasmic | 0 | False | Rossmann NAD(P)-binding GxGxxG; Cupin-like metal-binding HxH |  |  |
| 1390.SAMEA104415756.OFHT01000022_16 | 124 | 7.1 | 42.6 | lipoprotein | 0 | True |  |  |  |
| 1390.SAMEA104415756.OFHT01000022_17 | 263 | 9.3 | 40.1 | lipoprotein | 0 | True |  |  |  |
| 1390.SAMEA104415756.OFHT01000022_18 | 490 | 5.05 | 35.5 | cytoplasmic | 0 | False | Leucine zipper; Ser hydrolase nucleophile elbow GxSxG |  |  |
| 1390.SAMEA104415756.OFHT01000022_19 | 88 | 4.97 | 36.6 | cytoplasmic | 0 | False |  |  |  |
| 1390.SAMEA104415756.OFHT01000022_20 | 210 | 6.39 | 40.7 | cytoplasmic | 0 | False |  |  |  |
| 1390.SAMEA104415756.OFHT01000022_21 | 62 | 8.04 | 14.5 | cytoplasmic | 0 | False |  |  |  |
| 1390.SAMEA104415756.OFHT01000022_22 | 265 | 5.1 | 31.9 | cytoplasmic | 0 | False |  |  |  |
| 1390.SAMEA104415756.OFHT01000022_23 | 150 | 4.93 | 33.6 | cytoplasmic | 0 | False |  |  |  |
| 1390.SAMEA104415756.OFHT01000022_24 | 521 | 8.48 | 25.6 | secreted/periplasmic | 0 | True | ATP/GTP-binding P-loop (Walker A) x2 |  |  |
| 1390.SAMEA104415756.OFHT01000022_25 | 97 | 5.91 | 52.3 | membrane | 1 | False |  |  |  |
| 1390.SAMEA104415756.OFHT01000022_26 | 205 | 4.49 | 44.8 | secreted/periplasmic | 0 | True |  |  |  |
| 1390.SAMEA104415756.OFHT01000022_27 | 62 | 10 | 10.1 | secreted/periplasmic | 0 | True |  |  |  |
| 1390.SAMEA104415756.OFHT01000022_28 | 612 | 5.16 | 38.8 | cytoplasmic | 0 | False | ATP/GTP-binding P-loop (Walker A) |  |  |
| 1390.SAMEA104415756.OFHT01000022_29 | 106 | 9.61 | 28.1 | membrane | 1 | False |  |  |  |
| 1390.SAMEA104415756.OFHT01000022_30 | 69 | 9.78 | 38 | cytoplasmic | 0 | False |  |  |  |
| 1390.SAMEA104415756.OFHT01000022_31 | 207 | 5.03 | 31.7 | lipoprotein | 0 | True |  |  |  |
| 1390.SAMEA104415756.OFHT01000022_32 | 442 | 5.75 | 31.1 | cytoplasmic | 0 | False | ATP/GTP-binding P-loop (Walker A); Walker B (DEAD/DEAH-like) |  |  |
| 1390.SAMEA104415756.OFHT01000022_33 | 158 | 5.02 | 34 | cytoplasmic | 0 | False |  |  |  |
| 1390.SAMEA104415756.OFHT01000022_34 | 309 | 6.25 | 32.8 | cytoplasmic | 0 | False | Ser hydrolase nucleophile elbow GxSxG |  |  |
| 1390.SAMEA104415756.OFHT01000022_35 | 93 | 4.61 | 18.9 | cytoplasmic | 0 | False |  |  |  |
| 1390.SAMEA104415756.OFHT01000022_36 | 403 | 9.43 | 36.9 | membrane | 4 | True |  |  |  |
| 1390.SAMEA104415756.OFHT01000022_37 | 1148 | 5.42 | 36.7 | cytoplasmic | 0 | False |  |  |  |
| 1390.SAMEA104415756.OFHT01000022_38 | 308 | 9.99 | 26.6 | membrane | 3 | True |  |  |  |
| 1390.SAMEA104415756.OFHT01000022_39 | 305 | 9.37 | 31.3 | membrane | 2 | False |  |  |  |
| 1390.SAMEA104415756.OFHT01000022_40 | 356 | 8.83 | 27.8 | membrane | 2 | True |  |  |  |

## Files

- `results.tsv`
- `final_candidates.faa`
- `predicted_proteins.faa`
- `domains.tsv`
- `secondary_structure.txt`
- `motif_hits.tsv`
- `denovo_motifs.tsv`
