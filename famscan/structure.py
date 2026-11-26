"""Stage 3: tertiary-structure confirmation for the top hits.

Prediction / retrieval sources:
  * AlphaFold DB lookup by UniProt accession (instant, if the hit has one)
  * ESMFold public API (api.esmatlas.com) - single sequence, <= 400 aa,
    seconds per protein, no MSA needed
  * a folder of PDBs you produced yourself (ColabFold / AlphaFold3 /
    local ESMFold on a GPU) - pass `--pdb-dir`

Comparison:
  * TM-score against a reference structure with TM-align (`tmtools`, offline)
  * Foldseek easy-search against a local DB if `foldseek` is installed, or
    the public Foldseek web server for a quick check of a few structures.
"""
from __future__ import annotations

import re
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

import numpy as np
import requests
from Bio.PDB import PDBParser
from Bio.Data.IUPACData import protein_letters_3to1

ESMFOLD_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"
AFDB_URL = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"
FOLDSEEK_URL = "https://search.foldseek.com/api"
_UNIPROT = re.compile(r"\b([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})\b")


def uniprot_accession(identifier: str) -> str | None:
    """Pull a UniProt accession out of IDs like sp|P0A9K9|ADHE_ECOLI."""
    m = _UNIPROT.search(identifier.replace("|", " "))
    return m.group(1) if m else None


def fetch_afdb(acc: str, dest: Path) -> Path | None:
    r = requests.get(AFDB_URL.format(acc=acc), timeout=30)
    if r.status_code != 200 or not r.json():
        return None
    pdb_url = r.json()[0]["pdbUrl"]
    dest.write_text(requests.get(pdb_url, timeout=60).text)
    return dest


def esmfold(seq: str, dest: Path, retries: int = 3) -> Path:
    if len(seq) > 400:
        raise ValueError("ESMFold API limit is 400 aa; fold locally/ColabFold instead")
    for attempt in range(retries):
        r = requests.post(ESMFOLD_URL, data=seq, timeout=180)
        if r.status_code == 200 and r.text.startswith(("HEADER", "ATOM", "PARENT", "MODEL")):
            dest.write_text(r.text)
            return dest
        time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"ESMFold API failed: HTTP {r.status_code} {r.text[:200]}")


def plddt(pdb: Path) -> float:
    """Mean pLDDT from CA B-factors (0-100; ESMFold API's 0-1 scale is rescaled)."""
    vals = [float(l[60:66]) for l in pdb.read_text().splitlines()
            if l.startswith("ATOM") and l[12:16].strip() == "CA"]
    if not vals:
        return float("nan")
    m = float(np.mean(vals))
    return round(m * 100 if m <= 1.0 else m, 1)


def _ca_coords(pdb: Path) -> tuple[np.ndarray, str]:
    s = PDBParser(QUIET=True).get_structure("x", str(pdb))
    chain = next(s[0].get_chains())
    coords, seq = [], []
    for res in chain:
        if res.id[0] == " " and "CA" in res:
            coords.append(res["CA"].coord)
            seq.append(protein_letters_3to1.get(res.get_resname().capitalize(), "X"))
    return np.array(coords, dtype=float), "".join(seq)


def tm_score(model: Path, reference: Path) -> dict:
    """TM-align via tmtools. TM > 0.5 ~ same fold; > 0.8 ~ near-identical."""
    from tmtools import tm_align

    c1, s1 = _ca_coords(model)
    c2, s2 = _ca_coords(reference)
    r = tm_align(c1, c2, s1, s2)
    return {"tm_score": round(r.tm_norm_chain2, 3), "rmsd": round(r.rmsd, 2)}


def foldseek_local(pdbs: list[Path], target_db: str, out: Path, threads: int = 4) -> Path:
    if not shutil.which("foldseek"):
        raise RuntimeError("foldseek not installed (conda install -c bioconda foldseek)")
    with tempfile.TemporaryDirectory() as tmp:
        qdir = Path(tmp) / "q"
        qdir.mkdir()
        for p in pdbs:
            shutil.copy(p, qdir / p.name)
        subprocess.run(["foldseek", "easy-search", str(qdir), target_db, str(out),
                        str(Path(tmp) / "tmp"), "--threads", str(threads),
                        "--format-output", "query,target,evalue,bits,alntmscore,prob,theader"],
                       check=True)
    return out


def foldseek_web(pdb: Path, databases=("afdb50", "pdb100"), timeout: int = 600) -> list[dict]:
    """One query to the public Foldseek server. Be gentle: a few at a time."""
    files = {"q": (pdb.name, pdb.read_bytes())}
    data = [("mode", "3diaa")] + [("database[]", d) for d in databases]
    t = requests.post(f"{FOLDSEEK_URL}/ticket", files=files, data=data, timeout=60).json()
    tid, t0 = t["id"], time.time()
    while time.time() - t0 < timeout:
        st = requests.get(f"{FOLDSEEK_URL}/ticket/{tid}", timeout=30).json()["status"]
        if st == "COMPLETE":
            break
        if st == "ERROR":
            raise RuntimeError("Foldseek job failed")
        time.sleep(5)
    res = requests.get(f"{FOLDSEEK_URL}/result/{tid}/0", timeout=60).json()
    out = []
    for db in res.get("results", []):
        for aln in (db.get("alignments") or [[]])[0][:5]:
            out.append({"db": db["db"], "target": aln["target"], "evalue": aln["eval"],
                        "prob": aln.get("prob"), "tm": aln.get("alntmscore")})
    return out


def get_structure(hit_id: str, seq: str, outdir: Path, pdb_dir: Path | None = None,
                  method: str = "auto", log=print) -> tuple[Path | None, str]:
    """Return (pdb_path, source) for a hit, trying the cheapest source first.

    source is "provided", "alphafold-db", "esmfold-api", "unavailable" (this
    protein could not be folded) or "unreachable" (the servers are down, so
    there is no point trying the rest)."""
    unreachable = 0
    outdir.mkdir(parents=True, exist_ok=True)
    safe = re.sub(r"[^A-Za-z0-9_.-]", "_", hit_id)
    if pdb_dir:
        for ext in (".pdb", ".cif"):
            p = Path(pdb_dir) / f"{safe}{ext}"
            if p.exists():
                return p, "provided"
    if method in ("auto", "afdb"):
        acc = uniprot_accession(hit_id)
        if acc:
            try:
                p = fetch_afdb(acc, outdir / f"{safe}.afdb.pdb")
                if p:
                    return p, "alphafold-db"
            except requests.RequestException as e:
                unreachable += isinstance(e, (requests.ConnectionError, requests.Timeout))
                log(f"[structure] AFDB lookup failed for {hit_id}: {str(e)[:120]}")
    if method in ("auto", "esmfold"):
        try:
            return esmfold(seq, outdir / f"{safe}.esmfold.pdb"), "esmfold-api"
        except (requests.RequestException, RuntimeError, ValueError) as e:
            unreachable += isinstance(e, (requests.ConnectionError, requests.Timeout))
            log(f"[structure] ESMFold failed for {hit_id}: {str(e)[:120]}")
    return None, "unreachable" if unreachable else "unavailable"
