# famscan: Run a BLAST-style search (any direction)

_I just want a homology search: protein or DNA, translated or not_

- Stages: homology
- Parameters: `{'program': 'blastx', 'engine': 'builtin', 'evalue': 1e-05, 'query': 'data/bmyD_gene.fna', 'target': 'data/bacillus_proteins.faa'}`
- Run time: 0.1 s; final candidates: 6

> engine=ncbi uses NCBI's URL API: fine for a few queries; download databases for big jobs.

## Homology search

6 hits (blastx, engine builtin).

| query | target | evalue | score | method | domain_from | domain_to |
|---|---|---|---|---|---|---|
| ENA/CDG29762/CDG29762.1 | 1390.SAMEA104415756.OFHT01000022_139 | 1.2e-70 | 233 | blastx(phmmer) | 1 | 310 |
| ENA/CDG29762/CDG29762.1 | 1390.SAMEA104415756.OFHT01000022_3 | 9.7e-56 | 184 | blastx(phmmer) | 1 | 322 |
| ENA/CDG29762/CDG29762.1 | 1390.SAMEA104415756.OFHT01000022_264 | 2.7e-55 | 182 | blastx(phmmer) | 1 | 288 |
| ENA/CDG29762/CDG29762.1 | 1390.SAMEA104415756.OFHT01000022_266 | 7.5e-38 | 125 | blastx(phmmer) | 1 | 217 |
| ENA/CDG29762/CDG29762.1 | 1390.SAMEA104415756.OFHT01000022_265 | 5.0e-23 | 76.1 | blastx(phmmer) | 1 | 300 |
| ENA/CDG29762/CDG29762.1 | 1390.SAMEA104415756.OFHT01000022_267 | 5.6e-06 | 20 | blastx(phmmer) | 1 | 78 |

## Final candidates

| id | evalue | score |
|---|---|---|
| ENA/CDG29762/CDG29762.1->1390.SAMEA104415756.OFHT01000022_13 | 1.2e-70 | 233 |
| ENA/CDG29762/CDG29762.1->1390.SAMEA104415756.OFHT01000022_3 | 9.7e-56 | 184 |
| ENA/CDG29762/CDG29762.1->1390.SAMEA104415756.OFHT01000022_26 | 2.7e-55 | 182 |
| ENA/CDG29762/CDG29762.1->1390.SAMEA104415756.OFHT01000022_26 | 7.5e-38 | 125 |
| ENA/CDG29762/CDG29762.1->1390.SAMEA104415756.OFHT01000022_26 | 5.0e-23 | 76.1 |
| ENA/CDG29762/CDG29762.1->1390.SAMEA104415756.OFHT01000022_26 | 5.6e-06 | 20 |

## Files

- `results.tsv`
- `final_candidates.faa`
- `homology_hits.tsv`
