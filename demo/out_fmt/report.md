# famscan: Run a BLAST-style search (any direction)

_I just want a homology search: protein or DNA, translated or not_

- Stages: homology
- Parameters: `{'program': 'blastp', 'engine': 'builtin', 'timeout': 900, 'evalue': 1e-05, 'query': 'data/luxC_reference.faa', 'target': 'data/luxC_homologs.faa', 'formats': 'md'}`
- Run time: 0.2 s; final candidates: 12

> builtin needs no install and no network. ncbi/ebi/mmseqs are web services, fine for a few queries; for big jobs run `famscan fetchdb swissprot` and search the file locally.

## Homology search

12 hits (blastp, engine builtin).

| query | target | evalue | score | method | domain_from | domain_to |
|---|---|---|---|---|---|---|
| sp/P19841/LUXC_PHOPO | sp/P19841/LUXC_PHOPO | 0 | 1.12e+03 | blastp(phmmer) | 1 | 488 |
| sp/P19841/LUXC_PHOPO | sp/Q03324/LUXC1_PHOLE | 5.1e-285 | 935 | blastp(phmmer) | 1 | 478 |
| sp/P19841/LUXC_PHOPO | sp/P29236/LUXC2_PHOLE | 6.6e-274 | 899 | blastp(phmmer) | 1 | 478 |
| sp/P19841/LUXC_PHOPO | tr/B6ESM7/B6ESM7_ALISL | 2.5e-240 | 788 | blastp(phmmer) | 1 | 476 |
| sp/P19841/LUXC_PHOPO | sp/P12748/LUXC_ALIFS | 1.2e-226 | 743 | blastp(phmmer) | 1 | 479 |
| sp/P19841/LUXC_PHOPO | tr/Q5DZ03/Q5DZ03_ALIF1 | 6.4e-224 | 734 | blastp(phmmer) | 1 | 479 |
| sp/P19841/LUXC_PHOPO | tr/S3DGB3/S3DGB3_9GAMM | 2.5e-217 | 712 | blastp(phmmer) | 1 | 479 |
| sp/P19841/LUXC_PHOPO | sp/P23113/LUXC_PHOLU | 2.1e-211 | 692 | blastp(phmmer) | 1 | 480 |
| sp/P19841/LUXC_PHOPO | sp/Q7N577/LUXC_PHOLL | 2.3e-207 | 679 | blastp(phmmer) | 1 | 480 |
| sp/P19841/LUXC_PHOPO | sp/P08639/LUXC_VIBHA | 8.9e-203 | 664 | blastp(phmmer) | 1 | 477 |
| sp/P19841/LUXC_PHOPO | tr/A0A644X819/A0A644X819_9ZZZZ | 2.8e-17 | 52.3 | blastp(phmmer) | 48 | 427 |
| sp/P19841/LUXC_PHOPO | tr/A0A2H5XB72/A0A2H5XB72_9BACT | 6.5e-13 | 37.9 | blastp(phmmer) | 84 | 431 |

## Final candidates

| id | evalue | score |
|---|---|---|
| sp/P19841/LUXC_PHOPO->sp/P19841/LUXC_PHOPO | 0 | 1.12e+03 |
| sp/P19841/LUXC_PHOPO->sp/Q03324/LUXC1_PHOLE | 5.1e-285 | 935 |
| sp/P19841/LUXC_PHOPO->sp/P29236/LUXC2_PHOLE | 6.6e-274 | 899 |
| sp/P19841/LUXC_PHOPO->tr/B6ESM7/B6ESM7_ALISL | 2.5e-240 | 788 |
| sp/P19841/LUXC_PHOPO->sp/P12748/LUXC_ALIFS | 1.2e-226 | 743 |
| sp/P19841/LUXC_PHOPO->tr/Q5DZ03/Q5DZ03_ALIF1 | 6.4e-224 | 734 |
| sp/P19841/LUXC_PHOPO->tr/S3DGB3/S3DGB3_9GAMM | 2.5e-217 | 712 |
| sp/P19841/LUXC_PHOPO->sp/P23113/LUXC_PHOLU | 2.1e-211 | 692 |
| sp/P19841/LUXC_PHOPO->sp/Q7N577/LUXC_PHOLL | 2.3e-207 | 679 |
| sp/P19841/LUXC_PHOPO->sp/P08639/LUXC_VIBHA | 8.9e-203 | 664 |
| sp/P19841/LUXC_PHOPO->tr/A0A644X819/A0A644X819_9ZZZZ | 2.8e-17 | 52.3 |
| sp/P19841/LUXC_PHOPO->tr/A0A2H5XB72/A0A2H5XB72_9BACT | 6.5e-13 | 37.9 |

## Files

- `results.tsv`
- `final_candidates.faa`
- `homology_hits.tsv`
