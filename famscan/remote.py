"""Remote search services, as alternatives to NCBI's BLAST URL API.

    engine=ebi      EBI Job Dispatcher (same NCBI BLAST programs, separate queue,
                    usually minutes; requires an email address, which EBI uses to
                    contact you about misuse - it is sent only when you pass it)
    engine=mmseqs   MMseqs2 web server: seconds per query, profile-sensitive
    engine=ncbi     NCBI BLAST URL API (see search.remote_blast)

None of these are needed if you download a database: see famscan.databases.
Every call here has a wall-clock limit, so a stuck job cannot hang a run.
"""
from __future__ import annotations

import time
from pathlib import Path

import requests

EBI = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast"
MMSEQS = "https://search.mmseqs.com/api"
UA = {"User-Agent": "famscan (bioinformatics pipeline)"}

EBI_DATABASES = ("uniprotkb", "uniprotkb_swissprot", "uniprotkb_reference_proteomes",
                 "uniref50", "uniref90", "pdb", "afdb")


class SearchTimeout(RuntimeError):
    """A remote job did not finish inside the time limit. Carries the job id so
    the result can be collected later."""

    def __init__(self, service: str, job_id: str, waited: float):
        self.service, self.job_id, self.waited = service, job_id, waited
        super().__init__(f"{service} job {job_id} still running after {waited / 60:.0f} min. "
                         f"It is not lost - collect it later, or use a local database "
                         f"(famscan fetchdb swissprot) to avoid the queue.")


def _poll(status, done, failed, timeout: float, interval: float, service: str, job_id: str, log):
    t0 = time.time()
    waited = 0.0
    while time.time() - t0 < timeout:
        s = status()
        if s in done:
            return s
        if s in failed:
            raise RuntimeError(f"{service} job {job_id} failed with status {s}")
        time.sleep(min(interval, max(timeout - (time.time() - t0), 1)))
        interval = min(interval * 1.3, 30)
        if time.time() - t0 - waited > 60:
            waited = time.time() - t0
            log(f"      {service}: still {s} after {waited / 60:.0f} min")
    raise SearchTimeout(service, job_id, time.time() - t0)


# ------------------------------------------------------------------- EBI
def ebi_blast(sequence: str, program: str = "blastp", database: str = "uniprotkb_swissprot",
              email: str | None = None, evalue: float = 1e-5, hitlist_size: int = 100,
              timeout: float = 900, log=print) -> list[dict]:
    """One EBI Job Dispatcher BLAST job. `email` is required by the EBI and is
    sent to them with the job; famscan never fills it in for you."""
    if not email or "@" not in email:
        raise ValueError("the EBI requires a contact email for its BLAST service: pass email=...")
    stype = "dna" if program in ("blastn", "tblastx") else "protein"
    data = {"email": email, "program": program, "stype": stype, "database": database,
            "sequence": sequence, "exp": str(evalue), "alignments": str(hitlist_size),
            "scores": str(hitlist_size)}
    job_id = requests.post(f"{EBI}/run", data=data, timeout=60, headers=UA).text.strip()
    log(f"      ebi: submitted {job_id}")
    _poll(lambda: requests.get(f"{EBI}/status/{job_id}", timeout=30, headers=UA).text.strip(),
          {"FINISHED"}, {"FAILURE", "ERROR", "NOT_FOUND"}, timeout, 5, "ebi", job_id, log)
    out = requests.get(f"{EBI}/result/{job_id}/json", timeout=120, headers=UA).json()
    rows = []
    for hit in out.get("hits", []):
        hsp = (hit.get("hit_hsps") or [{}])[0]
        rows.append({"target": hit.get("hit_acc", ""), "description": (hit.get("hit_desc") or "")[:200],
                     "evalue": float(hsp.get("hsp_expect", 1)), "score": float(hsp.get("hsp_bit_score", 0)),
                     "pident": float(hsp.get("hsp_identity", 0)),
                     "domain_from": int(hsp.get("hsp_hit_from", 0) or 0),
                     "domain_to": int(hsp.get("hsp_hit_to", 0) or 0),
                     "sequence": (hsp.get("hsp_hseq") or "").replace("-", ""),
                     "method": f"ebi-{program}:{database}"})
    return rows


# ---------------------------------------------------------------- MMseqs2
def mmseqs_search(fasta_path: str | Path, databases=("UniProtKB/Swiss-Prot",),
                  timeout: float = 600, log=print) -> list[dict]:
    """MMseqs2 web server: fast, profile-sensitive protein search."""
    files = {"q": (Path(fasta_path).name, Path(fasta_path).read_bytes())}
    data = [("mode", "accept")] + [("database[]", d) for d in databases]
    t = requests.post(f"{MMSEQS}/ticket", files=files, data=data, timeout=60, headers=UA).json()
    job_id = t.get("id", "")
    if t.get("status") == "ERROR" or not job_id:
        raise RuntimeError(f"MMseqs2 server rejected the job: {t}")
    log(f"      mmseqs: submitted {job_id}")
    _poll(lambda: requests.get(f"{MMSEQS}/ticket/{job_id}", timeout=30, headers=UA).json().get("status"),
          {"COMPLETE"}, {"ERROR"}, timeout, 3, "mmseqs", job_id, log)
    res = requests.get(f"{MMSEQS}/result/{job_id}/0", timeout=120, headers=UA).json()
    rows = []
    for db in res.get("results", []):
        for aln in (db.get("alignments") or []):
            for a in (aln if isinstance(aln, list) else [aln]):
                rows.append({"target": a.get("target", ""), "query": a.get("query", ""),
                             "description": (a.get("theader") or "")[:200],
                             "evalue": float(a.get("eval", 1)), "score": float(a.get("score", 0)),
                             "pident": round(100 * float(a.get("seqId", 0)), 1),
                             "domain_from": int(a.get("dbStartPos", 0) or 0),
                             "domain_to": int(a.get("dbEndPos", 0) or 0),
                             "sequence": (a.get("tSeq") or "").replace("-", ""),
                             "method": f"mmseqs:{db.get('db', '')}"})
    return rows


# ------------------------------------------------------------------ NCBI
NCBI = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"


def ncbi_blast(sequence: str, program: str = "blastp", database: str = "nr", evalue: float = 1e-5,
               hitlist_size: int = 100, entrez_query: str | None = None, email: str | None = None,
               timeout: float = 900, log=print) -> list[dict]:
    """One NCBI BLAST job through the official URL API, with a wall-clock limit.

    Biopython's qblast polls forever, so a job stuck in NCBI's queue hangs the
    run; this raises SearchTimeout instead and reports the RID, which can be
    collected later at blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=<rid>.
    """
    import io

    from Bio.Blast import NCBIXML

    params = {"CMD": "Put", "PROGRAM": program, "DATABASE": database, "QUERY": sequence,
              "EXPECT": evalue, "HITLIST_SIZE": hitlist_size}
    if entrez_query:
        params["ENTREZ_QUERY"] = entrez_query
    if email:
        params["EMAIL"] = email
    r = requests.post(NCBI, data=params, timeout=60, headers=UA)
    r.raise_for_status()
    rid = next((line.split("=", 1)[1].strip() for line in r.text.splitlines() if line.strip().startswith("RID =")), "")
    if not rid:
        raise RuntimeError("NCBI did not return a request id (RID)")
    log(f"      ncbi: submitted {rid}")

    def status() -> str:
        s = requests.get(NCBI, params={"CMD": "Get", "FORMAT_OBJECT": "SearchInfo", "RID": rid},
                         timeout=30, headers=UA).text
        if "Status=WAITING" in s:
            return "WAITING"
        if "Status=FAILED" in s:
            return "FAILED"
        if "Status=UNKNOWN" in s:
            return "UNKNOWN"
        return "READY" if "Status=READY" in s else "WAITING"

    _poll(status, {"READY"}, {"FAILED", "UNKNOWN"}, timeout, 10, "ncbi", rid, log)
    xml = requests.get(NCBI, params={"CMD": "Get", "FORMAT_TYPE": "XML", "RID": rid},
                       timeout=180, headers=UA).text
    rec = NCBIXML.read(io.StringIO(xml))
    rows = []
    for aln in rec.alignments:
        hsp = aln.hsps[0]
        rows.append({"target": aln.accession, "description": aln.hit_def[:200],
                     "evalue": hsp.expect, "score": hsp.bits,
                     "pident": round(100 * hsp.identities / hsp.align_length, 1),
                     "domain_from": hsp.sbjct_start, "domain_to": hsp.sbjct_end,
                     "sequence": hsp.sbjct.replace("-", "") if program in ("blastp", "blastx") else "",
                     "method": f"ncbi-{program}:{database}"})
    return rows
