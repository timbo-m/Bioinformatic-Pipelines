"""Download and prepare local sequence databases, so no web service is needed.

    famscan fetchdb swissprot            ~90 MB  570k curated proteins
    famscan fetchdb pfam                 ~1.6 GB Pfam-A.hmm (domain architecture)
    famscan fetchdb uniref50             ~12 GB  clustered, broad coverage
    famscan fetchdb proteome:UP000005640 ~15 MB  one proteome (human here)

Downloaded once, a local search is thousands of sequences a minute with no
queue, no rate limit and no network. Files land in ~/.famscan/db by default
(override with --dir or $FAMSCAN_DB).
"""
from __future__ import annotations

import gzip
import os
import shutil
import subprocess
import time
from pathlib import Path

import requests

UNIPROT = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete"

CATALOG: dict[str, dict] = {
    "swissprot": {
        "url": f"{UNIPROT}/uniprot_sprot.fasta.gz",
        "filename": "uniprot_sprot.fasta",
        "size": "~90 MB compressed, ~290 MB unpacked",
        "about": "UniProtKB/Swiss-Prot: ~570,000 manually curated proteins. "
                 "The best default for naming proteins and for blastp.",
        "kind": "protein",
    },
    "trembl": {
        "url": f"{UNIPROT}/uniprot_trembl.fasta.gz",
        "filename": "uniprot_trembl.fasta",
        "size": "~60 GB compressed - only with a reason",
        "about": "UniProtKB/TrEMBL: all unreviewed entries. Huge; prefer uniref50.",
        "kind": "protein",
    },
    "uniref50": {
        "url": "https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz",
        "filename": "uniref50.fasta",
        "size": "~12 GB compressed",
        "about": "UniRef50: sequences clustered at 50 % identity. Broad coverage "
                 "at a fraction of the size of nr.",
        "kind": "protein",
    },
    "pfam": {
        "url": "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz",
        "filename": "Pfam-A.hmm",
        "size": "~1.6 GB unpacked",
        "about": "Pfam-A profile HMMs: pass as hmm_db for domain architecture.",
        "kind": "hmm",
    },
    "prosite": {
        "url": "https://ftp.expasy.org/databases/prosite/prosite.dat",
        "filename": "prosite.dat",
        "size": "~30 MB",
        "about": "Full PROSITE pattern library: pass as prosite_dat for motif scanning.",
        "kind": "patterns",
    },
}

PROTEOME_URL = ("https://rest.uniprot.org/uniprotkb/stream?format=fasta&compressed=true"
                "&query=%28proteome%3A{acc}%29")


def db_dir(override: str | Path | None = None) -> Path:
    d = Path(override or os.environ.get("FAMSCAN_DB") or Path.home() / ".famscan" / "db")
    d.mkdir(parents=True, exist_ok=True)
    return d


def describe() -> str:
    lines = [f"{k:<12} {v['size']:<34} {v['about']}" for k, v in CATALOG.items()]
    lines.append(f"{'proteome:ACC':<12} {'~5-30 MB':<34} "
                 "One UniProt proteome, e.g. proteome:UP000005640 (human), UP000000625 (E. coli K-12).")
    return "\n".join(lines)


def _download(url: str, dest: Path, log=print) -> Path:
    """Stream to <dest>.part, then move into place. Resumes are not attempted;
    an interrupted download leaves only the .part file."""
    part = dest.with_suffix(dest.suffix + ".part")
    headers = {"User-Agent": "famscan (bioinformatics pipeline; local database download)"}
    with requests.get(url, stream=True, timeout=120, headers=headers) as r:
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        done = last = 0
        t0 = time.time()
        with open(part, "wb") as fh:
            for chunk in r.iter_content(1 << 20):
                fh.write(chunk)
                done += len(chunk)
                if done - last > 50 << 20:
                    last = done
                    pct = f" ({100 * done / total:.0f} %)" if total else ""
                    log(f"      {done / 1e6:.0f} MB{pct} at {done / 1e6 / max(time.time() - t0, 1):.1f} MB/s")
    part.replace(dest)
    return dest


def fetch(name: str, directory: str | Path | None = None, force: bool = False,
          make_blast_db: bool = False, log=print) -> Path:
    """Download one catalog entry (or proteome:ACC) and unpack it. Returns the path."""
    out = db_dir(directory)
    if name.startswith("proteome:"):
        acc = name.split(":", 1)[1].strip().upper()
        entry = {"url": PROTEOME_URL.format(acc=acc), "filename": f"{acc}.fasta", "kind": "protein"}
    elif name in CATALOG:
        entry = CATALOG[name]
    else:
        raise ValueError(f"unknown database '{name}'. Available:\n{describe()}")

    final = out / entry["filename"]
    if final.exists() and not force:
        log(f"      already have {final} ({final.stat().st_size / 1e6:.0f} MB); --force to redownload")
    else:
        gz = entry["url"].endswith(".gz")
        tmp = out / (entry["filename"] + (".gz" if gz else ".tmp"))
        log(f"      downloading {entry['url']}")
        _download(entry["url"], tmp, log=log)
        if gz:
            log("      unpacking")
            with gzip.open(tmp, "rb") as fi, open(final, "wb") as fo:
                shutil.copyfileobj(fi, fo, 1 << 20)
            tmp.unlink()
        else:
            tmp.replace(final)
        log(f"      ready: {final} ({final.stat().st_size / 1e6:.0f} MB)")

    if make_blast_db and entry.get("kind") == "protein":
        prepare(final, log=log)
    return final


def prepare(fasta: str | Path, log=print) -> dict:
    """Build the optional index files that speed up repeat searches:
    a DIAMOND .dmnd and/or a BLAST+ database, when those tools are installed.
    The built-in HMMER engine needs none of this."""
    fasta = Path(fasta)
    made = {}
    if shutil.which("diamond"):
        dmnd = fasta.with_suffix(".dmnd")
        if not dmnd.exists():
            log("      building DIAMOND index")
            subprocess.run(["diamond", "makedb", "--in", str(fasta), "-d", str(dmnd.with_suffix(""))],
                           check=True, capture_output=True)
        made["diamond"] = str(dmnd)
    if shutil.which("makeblastdb") and not Path(str(fasta) + ".pin").exists():
        log("      building BLAST+ database")
        subprocess.run(["makeblastdb", "-in", str(fasta), "-dbtype", "prot"],
                       check=True, capture_output=True)
        made["blast+"] = str(fasta)
    if not made:
        log("      no DIAMOND/BLAST+ installed - the built-in HMMER engine needs no index")
    return made


def installed(directory: str | Path | None = None) -> list[dict]:
    out = db_dir(directory)
    rows = []
    for p in sorted(out.iterdir()):
        if p.suffix in (".part", ".gz") or p.is_dir():
            continue
        rows.append({"name": p.name, "path": str(p), "size_mb": round(p.stat().st_size / 1e6, 1)})
    return rows
