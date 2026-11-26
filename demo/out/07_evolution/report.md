# famscan: Study a family's evolution across species

_I want to compare a protein family across several species (copy number, conservation, tree)_

- Stages: load_many -> family_search -> primary -> features -> motifs -> msa -> tree
- Parameters: `{'evalue': 1e-05, 'tree_method': 'nj', 'query': 'data/PF02826.hmm', 'proteomes': ['data/proteome_938293.faa', 'data/bacillus_contig.fna', 'data/cp040672_cds.fna']}`
- Run time: 0.6 s; final candidates: 9

## Inputs

| proteome | kind | proteins |
|---|---|---|
| proteome_938293 | protein | 2100 |
| bacillus_contig | genome | 287 |
| cp040672_cds | cds | 100 |

## Family search

9 proteins match `data/PF02826.hmm` at E <= 1e-05.

Copies per proteome:

| proteome | copies |
|---|---|
| bacillus_contig | 1 |
| cp040672_cds | 1 |
| proteome_938293 | 7 |

## Primary structure

| index | length | mw_kda | pI | instability | gravy |
|---|---|---|---|---|---|
| mean | 351 | 38.6 | 5.26 | 32.6 | -0.16 |
| min | 296 | 31.8 | 4.71 | 22.8 | -0.5 |
| max | 468 | 50.8 | 5.9 | 41.1 | 0.16 |

## Sequence features (fast predictors)

Signal peptide: 0, >=1 TM helix: 0, coiled coil: 0, >30 % disordered: 0.

| localisation | count |
|---|---|
| cytoplasmic | 9 |

Per-residue H/E/C strings: `secondary_structure.txt`. These are quick heuristics; confirm with SignalP 6, DeepTMHMM, IUPred3 and PSIPRED.

## Motifs

10 motif matches.

| motif | count |
|---|---|
| Rossmann NAD(P)-binding GxGxxG | 6 |
| ATP/GTP-binding P-loop (Walker A) | 2 |
| Peroxisome targeting PTS1 | 1 |
| Ser hydrolase nucleophile elbow GxSxG | 1 |

**De novo motifs shared by the set** (x = any residue):

| motif | seq_frac | enrichment | length |
|---|---|---|---|
| LGxT | 0.56 | 20 | 4 |
| YExL | 0.56 | 20 | 4 |
| KxKE | 0.56 | 10 | 4 |
| DLxK | 0.78 | 7 | 4 |
| TAxG | 0.67 | 6 | 4 |
| ExLE | 0.67 | 4.8 | 4 |
| GAxV | 0.67 | 4.8 | 4 |
| LxEE | 0.67 | 4 | 4 |
| AKxK | 0.56 | 4 | 4 |
| SxAE | 0.56 | 4 | 4 |

## Alignment and conservation

9 sequences, 569 columns. 5 residues of `proteome_938293|938293.PRJEB85.HG003685_328` are in columns with conservation >= 0.85 (candidate active/binding-site residues):

G143, G145, G148, M218, G230

## Phylogeny

Newick: `tree.nwk` (open in iTOL or FigTree). For a maximum-likelihood tree run IQ-TREE on `alignment.afa`.

```
      ____________________ proteome_938293|938293.PRJEB85.HG0036...
    ,|
    ||________________ proteome_938293|938293.PRJEB85.HG0036...
  __|
 |  |  _______________ proteome_938293|938293.PRJEB85.HG0036...
 |  |_|
 |    |_________________ proteome_938293|938293.PRJEB85.HG0036...
 |
_|              _________ cp040672_cds|lcl|NZ_CP040672.1_cds_WP...
 |  ___________|
 |,|           |_________ proteome_938293|938293.PRJEB85.HG0036...
 |||
 ,||_____________________ proteome_938293|938293.PRJEB85.HG0036...
 ||
 ||_______________________ proteome_938293|938293.PRJEB85.HG0036...
 |
 |_______________________ bacillus_contig|1390.SAMEA104415756.O...

```

## Final candidates

| id | evalue | score | length | pI | instability | localisation_guess | n_tm | signal_peptide | motifs |
|---|---|---|---|---|---|---|---|---|---|
| proteome_938293/938293.PRJEB85.HG003685_328 | 6.6e-54 | 179 | 311 | 4.99 | 23.2 | cytoplasmic | 0 | False | ATP/GTP-binding P-loop (Walker A) |
| proteome_938293/938293.PRJEB85.HG003686_519 | 8.8e-53 | 176 | 322 | 4.71 | 36.4 | cytoplasmic | 0 | False |  |
| proteome_938293/938293.PRJEB85.HG003689_19 | 7.5e-52 | 172 | 327 | 4.98 | 22.8 | cytoplasmic | 0 | False | ATP/GTP-binding P-loop (Walker A); Rossmann NAD(P)-binding G |
| proteome_938293/938293.PRJEB85.HG003686_258 | 4.8e-42 | 140 | 328 | 5.2 | 37.3 | cytoplasmic | 0 | False |  |
| cp040672_cds/lcl/NZ_CP040672.1_cds_WP_003152011.1_73 | 3.3e-09 | 33.5 | 376 | 5.44 | 30 | cytoplasmic | 0 | False | Peroxisome targeting PTS1; Rossmann NAD(P)-binding GxGxxG; S |
| proteome_938293/938293.PRJEB85.HG003685_432 | 9.0e-09 | 32.1 | 368 | 5.16 | 35 | cytoplasmic | 0 | False | Rossmann NAD(P)-binding GxGxxG |
| proteome_938293/938293.PRJEB85.HG003691_27 | 3.6e-08 | 30.2 | 364 | 5.16 | 41 | cytoplasmic | 0 | False | Rossmann NAD(P)-binding GxGxxG |
| bacillus_contig/1390.SAMEA104415756.OFHT01000022_223 | 4.7e-08 | 29.8 | 296 | 5.9 | 26.9 | cytoplasmic | 0 | False | Rossmann NAD(P)-binding GxGxxG |
| proteome_938293/938293.PRJEB85.HG003688_14 | 5.2e-07 | 26.4 | 468 | 5.82 | 41.1 | cytoplasmic | 0 | False | Rossmann NAD(P)-binding GxGxxG |

## Files

- `results.tsv`
- `final_candidates.faa`
- `secondary_structure.txt`
- `motif_hits.tsv`
- `denovo_motifs.tsv`
- `alignment.afa`
- `conservation.tsv`
- `conserved_residues.tsv`
- `tree.nwk`
