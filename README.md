# famscan: a guided bioinformatics pipeline

Tell it what you're trying to do; it builds the full pipeline for that job,
from raw genome or proteome to ranked candidates with structure checks.

```
$ famscan start

famscan - what are you trying to do?

   1. Find members of a protein family
   2. Discover new enzymes
   3. Annotate a new genome
   4. Find drug targets                    (in the pathogen, not in the host)
   5. Screen for antibiotic resistance
   6. Screen for toxins and allergens
   7. Study a family's evolution across species
   8. Characterise unknown proteins
   9. Pick candidates for the lab
  10. Mine variants for protein engineering
  11. Compare strains
  12. Run a BLAST-style search (any direction)
```

It then shows the stages it will run, asks only for the inputs that workflow
needs, prints the equivalent one-line command, and runs.

## Install

```bash
pip install -e .            # biopython, pyhmmer, pyrodigal, pyfamsa, pandas, tmtools
pip install -e ".[mcp]"     # + MCP server so Claude can drive it
# optional external tools it will use when present:
conda install -c bioconda blast diamond foldseek meme
```

Everything except structure prediction and NCBI queries runs **offline**.

## Inputs

- **Protein FASTA**: used as is.
- **Genome (DNA FASTA/GenBank, contigs > 20 kb)**: genes are called with Prodigal (pyrodigal).
  For eukaryotes, give the predicted proteome (BRAKER/AUGUSTUS) instead.
- **CDS / gene FASTA (short DNA)**: translated directly.

## Stages (mix-and-match building blocks)

| Level | Stage | What it does |
|---|---|---|
| Search | `family_search` | HMMER profile search (Pfam accession, `.hmm`, or example sequences); DNA targets are 6-frame translated |
| | `homology` | BLAST-style search in **all five directions**, see below |
| | `name_by_homology` | names proteins from their best hit in e.g. Swiss-Prot |
| Primary structure | `primary` | ProtParam (MW, pI, instability, GRAVY, aromaticity, extinction coefficient, charge), composition, Wilkinson-Harrison solubility |
| Secondary structure | `features` | Chou-Fasman H/E/C prediction, structure-derived H/E/C when a model exists |
| Small-scale structure | `features`, `motifs`, `structure` | super-secondary motifs (helix-turn-helix, beta-hairpin, beta-alpha-beta, alpha-hairpin), coiled coils, TM helices, signal peptides (Sec, Tat, lipobox), disorder, low complexity, disulfide bonds from 3D |
| Motifs | `motifs` | PROSITE patterns (built-in set, or the full `prosite.dat`), your own regex/PROSITE motifs, de novo shared motifs (or MEME/STREME if installed) |
| Domains | `domains` | domain architecture against Pfam-A.hmm or your HMM library |
| Alignment | `compare_refs` | pairwise BLOSUM62 alignment to known members (or the HMM consensus) |
| | `msa` | multiple alignment (FAMSA), per-column conservation, conserved residues (candidate active/binding sites) |
| | `cluster` | redundancy removal at a chosen % identity (CD-HIT style) |
| Evolution | `tree` | neighbour-joining/UPGMA tree, Newick output |
| | `strain_compare` | reciprocal best hits, core/accessory/unique genes, presence-absence matrix |
| Screening | `subtract_host`, `resistance`, `toxin_allergen` | host-homolog removal and essentiality; CARD-style perfect/strict/loose classes; Codex allergen rules (> 35 % identity over 80 aa, identical 8-mers) |
| Engineering | `variants` | natural variation per position, back-to-consensus mutation suggestions |
| Tertiary structure | `structure` | AlphaFold DB, ESMFold API or your own PDBs; pLDDT, TM-score to a reference, disulfides, Foldseek |
| Decision | `filter`, `rank` | your property rules; composite lab-readiness score |

Sequence-only predictors (signal peptide, TM, disorder, coiled coil, secondary
structure) are fast first-pass heuristics for triaging many hits, labelled as
such in the report. Confirm final candidates with SignalP 6, DeepTMHMM,
IUPred3 and PSIPRED, or use the structure-derived calls.

## Search directions

| Program | Query | Target | Built-in engine (no install) |
|---|---|---|---|
| blastp | protein | protein | phmmer |
| blastn | DNA | DNA | nhmmer (both strands) |
| blastx | DNA (translated) | protein | 6-frame + phmmer |
| tblastn | protein | DNA (translated) | 6-frame + phmmer, hits mapped to genome coordinates |
| tblastx | DNA (translated) | DNA (translated) | 6-frame both + phmmer |
| hmmsearch | profile HMM | protein or DNA | hmmsearch |

```bash
famscan search -q gene.fna -t genome.fna -p tblastx
famscan search -q protein.faa -t nr -p blastp -e ncbi        # NCBI's servers, small batches
famscan search -q proteins.faa -t uniref50.dmnd -e diamond   # big local jobs
```

## Should BLAST go through the NCBI website in bulk?

No. Driving the web form with many jobs is slow (each job queues for minutes),
breaks NCBI's usage guidelines (about 1 request every 10 s, ~100 searches a day,
heavy users get blocked), and breaks when the page changes.

| Batch size | Use |
|---|---|
| 1-20 queries | `-e ncbi`: NCBI's official URL API, rate limited politely |
| many queries or whole proteomes | download the database once and search locally (built-in engine, BLAST+ or DIAMOND): thousands of sequences a minute, no limits |

BLAST+ is US-government public-domain software and NCBI provides its
databases for download precisely for large local jobs. UniProt data is CC BY
4.0 (credit it); Pfam/InterPro is CC0. (General information, not legal advice.)

Structure prediction has the same split: the ESMFold and AlphaFold DB APIs are
fine for tens of proteins; for hundreds run ColabFold or ESMFold on a GPU and
pass the folder as `pdb_dir`.

## Non-interactive use

```bash
famscan list                                       # workflows and their inputs
famscan workflow enzyme --out out \
  --set input=genome.fna --set query=PF00561 \
  --set "motifs=nucleophile=G-x-S-x-G" --set hmm_db=Pfam-A.hmm --set top_n=10
famscan workflow drug_target --set input=pathogen.faa --set host=human.faa --set essential_db=deg.faa
famscan workflow evolution --set query=PF02826 --set "proteomes=sp1.faa sp2.faa sp3.fna"
```

Every run writes `report.md`, `results.tsv` (one row per protein, every
computed column), `final_candidates.faa`, plus stage files (alignment,
tree, domains, motifs, variants, structures).

Want Excel or PDF too? Add `--format`:

```
famscan workflow family --set input=proteome.faa --set query=PF02826 --format md,xlsx,pdf   # or --format all
famscan export famscan_out --format xlsx,pdf                                                # convert an earlier run
```

`report.xlsx` has a Summary sheet (parameters + report text) and one sheet per
table (`results`, motifs, domains, ...) with filters and frozen headers;
`report.pdf` is the rendered report (wide tables show their first 12 columns;
every column is in the Excel file). The interactive menu (`famscan start`)
asks which formats you want.

## Databases you may want

| For | Download |
|---|---|
| domains | Pfam-A.hmm: https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/ |
| naming | Swiss-Prot: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz |
| host proteome | human UP000005640 from UniProt proteomes |
| resistance | CARD: https://card.mcmaster.ca/download |
| essential genes | DEG: http://origin.tubic.org/deg/ |
| motifs | prosite.dat: https://ftp.expasy.org/databases/prosite/ |

## Use from Claude

- **Claude Code skill**: copy `skill/` to `~/.claude/skills/famscan/`.
- **MCP server**: `claude mcp add famscan -- python -m famscan.mcp_server`, which exposes
  `list_workflows`, `run`, `blast_search`, `protein_profile` and `blosum_compare`.

## Demo

`bash demo/run_demos.sh` runs all 12 workflows and all 5 search directions on
the bundled data (`demo/out/*/report.md`, about 2 minutes). The data comes from the pyhmmer
test set (MIT): a 2,100-protein bacterial proteome (taxon 938293), a 391 kb
*Bacillus* contig (taxon 1390), 12 LuxC homologs and Pfam PF02826. The screening
demos use **stand-in databases** (a 5-HMM mini library in place of Pfam, the
*Bacillus* proteins in place of the human proteome, mutated copies of a few proteins
in place of CARD/toxin/allergen sets), so they show how the pipeline works rather
than real findings.

## Tests

`python -m pytest tests/` (offline, about 5 s).
