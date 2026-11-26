# famscan: Run a BLAST-style search (any direction)

_I just want a homology search: protein or DNA, translated or not_

- Stages: homology
- Parameters: `{'program': 'blastn', 'engine': 'builtin', 'evalue': 1e-05, 'query': 'data/bmyD_gene.fna', 'target': 'data/bacillus_contig.fna'}`
- Run time: 1.5 s; final candidates: 1

> engine=ncbi uses NCBI's URL API: fine for a few queries; download databases for big jobs.

## Homology search

4 hits (blastn, engine builtin).

| query | target | evalue | score | method | domain_from | domain_to | strand |
|---|---|---|---|---|---|---|---|
| ENA/CDG29762/CDG29762.1 | 1390.SAMEA104415756.OFHT01000022 | 5.4e-29 | 96.3 | blastn(nhmmer) | 185567 | 186450 | + |
| ENA/CDG29762/CDG29762.1 | 1390.SAMEA104415756.OFHT01000022 | 7.8e-21 | 69.4 | blastn(nhmmer) | 315274 | 316130 | + |
| ENA/CDG29762/CDG29762.1 | 1390.SAMEA104415756.OFHT01000022 | 1.8e-17 | 58.3 | blastn(nhmmer) | 313291 | 314140 | + |
| ENA/CDG29762/CDG29762.1 | 1390.SAMEA104415756.OFHT01000022 | 1.4e-15 | 52 | blastn(nhmmer) | 1392 | 2225 | + |

## Final candidates

| id | evalue | score |
|---|---|---|
| ENA/CDG29762/CDG29762.1->1390.SAMEA104415756.OFHT01000022 | 5.4e-29 | 96.3 |
| ENA/CDG29762/CDG29762.1->1390.SAMEA104415756.OFHT01000022 | 7.8e-21 | 69.4 |
| ENA/CDG29762/CDG29762.1->1390.SAMEA104415756.OFHT01000022 | 1.8e-17 | 58.3 |
| ENA/CDG29762/CDG29762.1->1390.SAMEA104415756.OFHT01000022 | 1.4e-15 | 52 |

## Files

- `results.tsv`
- `final_candidates.faa`
- `homology_hits.tsv`
