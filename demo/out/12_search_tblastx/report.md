# famscan: Run a BLAST-style search (any direction)

_I just want a homology search: protein or DNA, translated or not_

- Stages: homology
- Parameters: `{'program': 'tblastx', 'engine': 'builtin', 'evalue': 1e-05, 'query': 'data/bmyD_gene.fna', 'target': 'data/bacillus_contig.fna'}`
- Run time: 0.4 s; final candidates: 1

> engine=ncbi uses NCBI's URL API: fine for a few queries; download databases for big jobs.

## Homology search

1 hits (tblastx, engine builtin).

| query | target | evalue | score | method | domain_from | domain_to | frame | nt_start | nt_end |
|---|---|---|---|---|---|---|---|---|---|
| ENA/CDG29762/CDG29762.1 | 1390.SAMEA104415756.OFHT01000022 | 3.3e-69 | 233 | tblastx(phmmer) | 7 | 316 | f+2 | 185567 | 186496 |

## Final candidates

| id | evalue | score |
|---|---|---|
| ENA/CDG29762/CDG29762.1->1390.SAMEA104415756.OFHT01000022 | 3.3e-69 | 233 |

## Files

- `results.tsv`
- `final_candidates.faa`
- `homology_hits.tsv`
