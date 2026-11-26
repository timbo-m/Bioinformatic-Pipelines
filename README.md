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

## Try it in two minutes

```bash
git clone https://github.com/timbo-m/Bioinformatic-Pipelines
cd Bioinformatic-Pipelines
python -m venv .venv && source .venv/bin/activate     # Windows: .venv\Scripts\activate
pip install -e .
famscan demo                                          # three worked examples, ~30 s, no network
famscan start                                         # then your own data
```

`famscan demo` runs on the bundled data and needs no internet, no databases and
no other programs, so a tester can see real output before deciding whether to
set anything up. `famscan list` shows every workflow and the inputs it needs.

Needs Python 3.10 or newer. Tested on Linux and Windows 11.

<details>
<summary>If something goes wrong (notes from real installs)</summary>

- **Windows and `tmtools`**: it only has Windows builds up to Python 3.12, so use
  3.12 there, or skip it (`pip install -e . --no-deps` then install the rest) if
  you don't need TM-score comparisons.
- **Windows Smart App Control** may block some compiled packages or the
  `famscan.exe` shim. Run `python -m famscan.cli ...` instead of `famscan ...`,
  and pin `pyhmmer==0.12.2` if its wheel is blocked. Don't turn the setting off.
- **`famscan: command not found`**: the virtual environment isn't active, or use
  `python -m famscan.cli`.
- **A web service is unreachable**: that's expected on a locked-down network. The
  run continues without that step and says so; see "When a web service is down".

Please report anything else as a GitHub issue, with the command you ran and the
last few lines of output.
</details>

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

## When a web service is down

Nothing is lost if ESMFold, AlphaFold DB, NCBI, the EBI or MMseqs2 is unreachable.
Only the stage that needed the service is affected; every other stage keeps its
results. `famscan start` tells you which service failed and asks whether to
continue without that step. Scripts choose up front:

```bash
--set on_service_down=ask    # explain and ask (the default in famscan start)
--set on_service_down=skip   # carry on without that step (the default for scripts)
--set on_service_down=stop   # stop the run
```

Either way the report gains a section naming the service, the stage, and exactly
what is missing from the results, so a gap is never silent.

## Avoiding web services entirely

```bash
famscan fetchdb --list                 # what's available, what you already have
famscan fetchdb swissprot              # ~90 MB, 570k curated proteins
famscan fetchdb proteome:UP000005640   # one proteome (human here)
famscan fetchdb pfam                   # Pfam-A.hmm for domain architecture
```

Downloads land in `~/.famscan/db` (or `$FAMSCAN_DB`, or `--dir`), and a DIAMOND
or BLAST+ index is built when those tools are installed. After that, searches
are local: thousands of sequences a minute, no queue, no rate limit, no network,
and the same answer every time. `famscan start` offers the download when an
input needs a database you don't have, and reuses one you already downloaded.

## Should BLAST go through the NCBI website in bulk?

No. Driving the web form with many jobs is slow (each job queues for minutes),
breaks NCBI's usage guidelines (about 1 request every 10 s, ~100 searches a day,
heavy users get blocked), and breaks when the page changes.

| Batch size | Use |
|---|---|
| 1-20 queries | a web service: `-e ncbi` (official URL API), `-e ebi` (separate queue, usually quicker; the EBI requires a contact email, which
famscan asks for and can save with `famscan config --email you@example.org` - it
is kept on your computer only, never committed, and written to reports as
`email: (set)` so sharing a report or spreadsheet does not share your address) or `-e mmseqs` (seconds per query) |
| many queries or whole proteomes | `famscan fetchdb swissprot` once, then search locally with the built-in engine, BLAST+ or DIAMOND |

Every web service call has a time limit (`--set timeout=900` by default) and
reports its job id if it runs out, so a job stuck in a queue can never hang a
run. NCBI queues can reach 30 minutes at busy times; one local search of
Swiss-Prot takes seconds.

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

## Output formats

Every run writes `report.md`, `results.tsv` (one row per protein, every computed
column), `final_candidates.faa`, plus stage files (alignment, tree, domains,
motifs, variants, structures). Add `--format` for more:

```bash
famscan workflow family --format all   ...   # md + xlsx + pdf
famscan workflow family --format excel ...   # md + xlsx
famscan search -q a.faa -t b.faa --format pdf ...
famscan export out/my_run --format xlsx,pdf      # convert a finished run
```

- **xlsx**: a Summary sheet, a Results sheet (one row per protein, frozen header
  and filters on every column) and one sheet per stage table.
- **pdf**: the report typeset in landscape with headings and tables.

`famscan start` asks which formats you want; Markdown is always written.

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

`python -m pytest tests/` (24 tests, offline, about 6 s).

## Status

Verified end to end on Linux and Windows 11: all 12 workflows, all 5 search
directions, Swiss-Prot download and local search, ESMFold, AlphaFold DB,
Foldseek, TM-score and the EBI engine. The MMseqs2 web server is down at
their end (ColabFold hardware move), so that engine is tested only against a
simulated outage. Fast sequence-only predictors are first-pass filters, not
substitutes for SignalP 6, DeepTMHMM, IUPred3 or PSIPRED.

Contributions and bug reports welcome: please include the command and the last
few lines of output.
