# famscan: Mine variants for protein engineering

_I want to improve a protein (stability, activity) using natural sequence variation_

- Stages: load -> family_search -> cluster -> msa -> variants -> tree
- Parameters: `{'cluster_identity': 95, 'min_variant_freq': 0.15, 'evalue': 1e-05, 'reference': 'data/luxC_reference.faa', 'input': 'data/luxC_homologs.faa', 'query': 'data/luxC_reference.faa', 'reference_id': 'sp|P19841|LUXC_PHOPO'}`
- Run time: 0.5 s; final candidates: 12

> `input` is the sequence database to mine homologs from (e.g. UniRef50 subset).

## Input

`data/luxC_homologs.faa` read as **protein**: 12 proteins.

## Family search

12 proteins match `data/luxC_reference.faa` at E <= 1e-05.

## Clustering

12 sequences -> 12 clusters at 95 % identity; continuing with one representative per cluster.

## Alignment and conservation

12 sequences, 560 columns. 106 residues of `sp|P19841|LUXC_PHOPO` are in columns with conservation >= 0.85 (candidate active/binding-site residues):

M11, G21, I68, N70, F71, L72, Q77, W79, E82, Y84, R87, L94, A105, L116, S118, L122, L130, G131, D137, W139, P141, Y146, A149, G153, H157, L158, L159, A160, G161, N162, G167, S170, I171, L176, T177, N179, I182, K184, S186, S187, P190, T192, S199, F200, I201, D223, D235, V237, W240, G241, G242, I246, W248, A249, P254, D258, K261, G263, K265, S267, A281, G283, A285, D287, I288, C289, D292, Q293, C296, S298, Q300, Y303, G306, F312, L320, L327, P328, F334, E346, V354

## Variant mining

266 positions of `sp|P19841|LUXC_PHOPO` vary among 11 homologs. 23 back-to-consensus mutations suggested (a common way to gain stability):

| suggested | consensus_freq | conservation | alternatives |
|---|---|---|---|
| R24Q | 0.75 | 0.538 | Q:0.75,F:0.25 |
| Y39N | 0.67 | 0.53 | N:0.67,K:0.22 |
| K78R | 0.73 | 0.788 | R:0.73 |
| N83E | 0.73 | 0.71 | E:0.73,T:0.18 |
| L88R | 0.8 | 0.737 | R:0.8 |
| S117C | 0.9 | 0.733 | C:0.9 |
| I172L | 0.73 | 0.725 | L:0.73 |
| T203V | 0.73 | 0.672 | V:0.73 |
| S210T | 0.8 | 0.737 | T:0.8 |
| I226L | 0.73 | 0.725 | L:0.73 |
| S239A | 0.82 | 0.721 | A:0.82 |
| E251K | 0.73 | 0.725 | K:0.73 |
| N256H | 0.64 | 0.588 | H:0.64,D:0.18 |
| S280A | 0.82 | 0.721 | A:0.82 |
| S339A | 0.73 | 0.633 | A:0.73 |
| Q348L | 0.64 | 0.641 | L:0.64 |
| K351G | 0.82 | 0.812 | G:0.82 |
| A381V | 0.73 | 0.788 | V:0.73 |
| N398K | 0.73 | 0.725 | K:0.73 |
| T405S | 0.73 | 0.633 | S:0.73 |
| V406I | 0.64 | 0.588 | I:0.64 |
| S464T | 0.9 | 0.733 | T:0.9 |
| K471E | 0.7 | 0.654 | E:0.7 |

## Phylogeny

Newick: `tree.nwk` (open in iTOL or FigTree). For a maximum-likelihood tree run IQ-TREE on `alignment.afa`.

```
  __________________________________ tr|A0A2H5XB72|A0A2H5XB72_9BACT
 |
_|_________________________________ tr|A0A644X819|A0A644X819_9ZZZZ
 |
 |                           __ sp|P23113|LUXC_PHOLU
 |                  ________|
 |                ,|        |__ sp|Q7N577|LUXC_PHOLL
 |                ||
 |                ||____________ sp|P08639|LUXC_VIBHA
 |________________|
                  | ___________ tr|S3DGB3|S3DGB3_9GAMM
                  ||
                  ||     ______ tr|B6ESM7|B6ESM7_ALISL
                  ||   _|
                   |  | |     ___ sp|P12748|LUXC_ALIFS
                   |  | |____|
                   |__|      |___ tr|Q5DZ03|Q5DZ03_ALIF1
                      |
                      |     ___ sp|Q03324|LUXC1_PHOLE
                      |    |
                      |____|____ sp|P29236|LUXC2_PHOLE
                           |
                           |_____ sp|P19841|LUXC_PHOPO

```

## Final candidates

| id | evalue | score | length |
|---|---|---|---|
| sp/P19841/LUXC_PHOPO | 0 | 1.12e+03 | 488 |
| sp/Q03324/LUXC1_PHOLE | 5.1e-285 | 935 | 478 |
| sp/P29236/LUXC2_PHOLE | 6.6e-274 | 898 | 478 |
| tr/B6ESM7/B6ESM7_ALISL | 2.5e-240 | 788 | 476 |
| sp/P12748/LUXC_ALIFS | 1.2e-226 | 743 | 479 |
| tr/Q5DZ03/Q5DZ03_ALIF1 | 6.4e-224 | 734 | 479 |
| tr/S3DGB3/S3DGB3_9GAMM | 2.5e-217 | 712 | 479 |
| sp/P23113/LUXC_PHOLU | 2.1e-211 | 692 | 480 |
| sp/Q7N577/LUXC_PHOLL | 2.3e-207 | 679 | 480 |
| sp/P08639/LUXC_VIBHA | 8.9e-203 | 664 | 477 |
| tr/A0A644X819/A0A644X819_9ZZZZ | 2.8e-17 | 52.3 | 427 |
| tr/A0A2H5XB72/A0A2H5XB72_9BACT | 6.5e-13 | 37.9 | 502 |

## Files

- `results.tsv`
- `final_candidates.faa`
- `alignment.afa`
- `conservation.tsv`
- `conserved_residues.tsv`
- `variants.tsv`
- `tree.nwk`
