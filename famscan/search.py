"""Homology search in every direction BLAST supports, plus profile (HMM) search.

    program   query        target        built-in engine (pyhmmer, no install)
    -------   ----------   -----------   --------------------------------------
    blastp    protein      protein       phmmer
    blastn    nucleotide   nucleotide    nhmmer (both strands)
    blastx    nucleotide   protein       6-frame translate query  -> phmmer
    tblastn   protein      nucleotide    6-frame translate target -> phmmer
    tblastx   nucleotide   nucleotide    translate both           -> phmmer
    hmmsearch profile HMM  protein/nt    hmmsearch (nt targets are 6-frame translated)

Engines:
  * "builtin" - pyhmmer, in-process and multithreaded. Default; no installs.
  * "blast+"  - local NCBI BLAST+ (or DIAMOND for blastp/blastx with a .dmnd db).
  * "ncbi"    - NCBI's official URL API via Biopython. Small batches only: NCBI
                asks for <= 1 request per 10 s and ~100 searches/day. Never
                automate the web form instead.
"""
from __future__ import annotations

import gzip
import io
import re
import shutil
import subprocess
import tempfile
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Iterable

import pyhmmer

from .seqio import alphabet_of, read_seqs, six_frame, write_fasta

PROGRAMS = {
    #           query       target      translate_query translate_target
    "blastp":  ("protein", "protein", False, False),
    "blastn":  ("dna", "dna", False, False),
    "blastx":  ("dna", "protein", True, False),
    "tblastn": ("protein", "dna", False, True),
    "tblastx": ("dna", "dna", True, True),
}


@dataclass
class Hit:
    target: str
    query: str
    evalue: float
    score: float
    method: str
    description: str = ""
    domain_from: int | None = None
    domain_to: int | None = None
    sequence: str = ""
    extra: dict = field(default_factory=dict)

    def row(self) -> dict:
        d = asdict(self)
        d.pop("extra")
        d.pop("sequence")
        d.update(self.extra)
        return d


def read_fasta(path) -> dict[str, str]:  # backwards compatible name
    return read_seqs(path)


def _s(x) -> str:
    return x.decode() if isinstance(x, (bytes, bytearray)) else str(x)


def _block(seqs: dict[str, str], kind: str):
    alpha = pyhmmer.easel.Alphabet.dna() if kind == "dna" else pyhmmer.easel.Alphabet.amino()
    clean = re.compile(r"[^A-Z]") if kind == "protein" else re.compile(r"[^ACGTURYSWKMBDHVN]")
    ts = [pyhmmer.easel.TextSequence(name=k, sequence=clean.sub("X" if kind == "protein" else "N", v))
          for k, v in seqs.items() if v]
    return pyhmmer.easel.TextSequenceBlock(ts).digitize(alpha), alpha


def guess_program(query_kind: str, target_kind: str, translated: bool = False) -> str:
    if query_kind == "protein":
        return "blastp" if target_kind == "protein" else "tblastn"
    if target_kind == "protein":
        return "blastx"
    return "tblastx" if translated else "blastn"


# ----------------------------------------------------------------- public API
ENGINES = ("builtin", "blast+", "diamond", "ncbi", "ebi", "mmseqs")


def remote_search(query: str, program: str, engine: str, database: str, evalue: float,
                  email: str | None = None, timeout: float = 900, pause: float = 11.0,
                  max_hits: int = 100, log=print) -> list[Hit]:
    """One job per query against a web service (ncbi, ebi or mmseqs). Every
    service has a time limit, so a queued job cannot hang the run. For more than
    about 20 queries download a database instead: `famscan fetchdb swissprot`."""
    from . import remote as rs

    queries = read_seqs(query)
    if len(queries) > 20:
        log(f"      note: {len(queries)} queries through {engine}; a local database "
            f"(famscan fetchdb swissprot) would be far faster and has no limits")
    hits: list[Hit] = []
    if engine == "mmseqs":  # takes the whole FASTA in one job
        for r in rs.mmseqs_search(query, timeout=timeout, log=log):
            hits.append(_from_row(r, r.get("query") or next(iter(queries))))
        return hits
    fn = rs.ebi_blast if engine == "ebi" else rs.ncbi_blast
    kw = {"email": email} if engine == "ebi" or email else {}
    for i, (qid, seq) in enumerate(queries.items()):
        if i:
            time.sleep(pause)  # NCBI asks for <= 1 request / 10 s
        log(f"      {engine}: {program} {qid} ({i + 1}/{len(queries)}) vs {database}")
        for r in fn(seq, program=program, database=database, evalue=evalue,
                    hitlist_size=max_hits, timeout=timeout, log=log, **kw):
            hits.append(_from_row(r, qid))
    return hits


def _from_row(r: dict, qid: str) -> Hit:
    return Hit(target=r["target"], query=qid, evalue=r["evalue"], score=r["score"],
               method=r["method"], description=r.get("description", ""),
               domain_from=r.get("domain_from"), domain_to=r.get("domain_to"),
               sequence=r.get("sequence", ""), extra={"pident": r.get("pident")})


def homology_search(query: str, target: str, program: str = "auto", engine: str = "builtin",
                    evalue: float = 1e-5, cpus: int = 0, max_hits: int = 500,
                    workdir: Path | None = None, email: str | None = None,
                    timeout: float = 900, log=print) -> list[Hit]:
    """Run any BLAST-style search. `query`/`target` are FASTA paths (target may
    also be a BLAST db name for engine='ncbi' or 'blast+')."""
    if program == "hmmsearch" or _is_hmm(query) or re.fullmatch(r"PF\d{5}(\.\d+)?", query or ""):
        return hmmer_search(query, target, evalue=evalue, cpus=cpus, workdir=workdir)
    qkind = alphabet_of(query)
    if program == "auto":
        tkind = alphabet_of(target) if Path(target).exists() else "protein"
        program = guess_program(qkind, tkind)
    if program not in PROGRAMS:
        raise ValueError(f"unknown program {program}; choose {list(PROGRAMS)} or hmmsearch")
    if engine in ("ncbi", "ebi", "mmseqs"):
        return remote_search(query, program, engine, target, evalue, email=email,
                             timeout=timeout, max_hits=min(max_hits, 100), log=log)
    if engine in ("blast+", "diamond"):
        return local_blast(query, target, program=program, evalue=evalue, max_targets=max_hits,
                           tool="diamond" if engine == "diamond" else "auto")
    return builtin_search(query, target, program, evalue=evalue, cpus=cpus)


def builtin_search(query: str, target: str, program: str, evalue: float = 1e-5,
                   cpus: int = 0) -> list[Hit]:
    qk, tk, tq, tt = PROGRAMS[program]
    q, t = read_seqs(query), read_seqs(target)
    if program == "blastn":
        qb, alpha = _block(q, "dna")
        tb, _ = _block(t, "dna")
        hits = []
        for top in pyhmmer.nhmmer(qb, tb, cpus=cpus, E=evalue):
            qn = _s(top.query.name)
            for h in top:
                if not h.included:
                    continue
                d = h.best_domain
                a, b = d.env_from, d.env_to
                hits.append(Hit(target=_s(h.name), query=qn, evalue=h.evalue, score=h.score,
                                method="blastn(nhmmer)", description=_s(h.description or ""),
                                domain_from=min(a, b), domain_to=max(a, b),
                                extra={"strand": "+" if a <= b else "-"}))
        return hits
    qp = six_frame(q) if tq else q
    tp = six_frame(t) if tt else t
    qb, _ = _block(qp, "protein")
    tb, _ = _block(tp, "protein")
    hits = []
    for top in pyhmmer.phmmer(qb, tb, cpus=cpus, E=evalue):
        qn = _s(top.query.name)
        for h in top:
            if not h.included:
                continue
            qparts = qn.rsplit("|", 2) if tq else [qn]
            hits.append(_mk_hit(h, qparts[0], f"{program}(phmmer)", tp, translated_target=tt,
                                extra={"query_frame": qparts[1], "query_nt": qparts[2]} if tq else {}))
    return dedupe(hits, key=lambda h: (h.query, h.target))


def _mk_hit(h, qname: str, method: str, tseqs: dict, translated_target: bool, extra=None) -> Hit:
    name = _s(h.name)
    d = h.best_domain
    extra = dict(extra or {})
    target = name
    if translated_target:  # map peptide coords back to the nucleotide record
        sid, frame, span = name.rsplit("|", 2)
        s0, e0 = map(int, span.split("-"))
        fr = int(frame[1:])
        if fr > 0:
            nt_a, nt_b = s0 + 3 * (d.env_from - 1), s0 + 3 * d.env_to - 1
        else:
            nt_a, nt_b = e0 - 3 * d.env_to + 1, e0 - 3 * (d.env_from - 1)
        extra.update(frame=frame, nt_start=nt_a, nt_end=nt_b)
        target = sid
    return Hit(target=target, query=qname, evalue=h.evalue, score=h.score, method=method,
               description=_s(h.description or ""), domain_from=d.env_from, domain_to=d.env_to,
               sequence=tseqs.get(name, ""), extra=extra)


# --------------------------------------------------------------------- HMMER
def _is_hmm(path) -> bool:
    try:
        return Path(path).exists() and Path(path).read_text(errors="ignore")[:10].startswith("HMMER")
    except (OSError, TypeError):
        return False


def fetch_pfam_hmm(accession: str, dest: Path) -> Path:
    """Download a Pfam HMM from InterPro (needs network)."""
    import requests

    acc = accession.split(".")[0].upper()
    r = requests.get(f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{acc}?annotation=hmm", timeout=60)
    r.raise_for_status()
    data = gzip.decompress(r.content) if r.content[:2] == b"\x1f\x8b" else r.content
    dest.write_bytes(data)
    return dest


def load_hmms(query: str, workdir: Path | None = None) -> list:
    qpath = Path(query)
    if not qpath.exists() and query.upper().startswith("PF"):
        workdir = Path(workdir or tempfile.mkdtemp())
        workdir.mkdir(parents=True, exist_ok=True)
        qpath = fetch_pfam_hmm(query, workdir / f"{query}.hmm")
    with pyhmmer.plan7.HMMFile(str(qpath)) as hf:
        return list(hf)


def hmmer_search(query: str, target: str | dict, evalue: float = 1e-5, cpus: int = 0,
                 workdir: Path | None = None) -> list[Hit]:
    """Profile search. `query`: HMM file, Pfam accession or protein FASTA (phmmer).
    `target`: protein or nucleotide FASTA path, or a {id: protein} dict."""
    seqs = read_seqs(target) if not isinstance(target, dict) else target
    translated = alphabet_of(seqs) == "dna"
    tseqs = six_frame(seqs) if translated else seqs
    tb, alpha = _block(tseqs, "protein")
    if _is_hmm(query) or not Path(query).exists():
        queries, method = load_hmms(query, workdir), "hmmsearch"
        results = pyhmmer.hmmsearch(queries, tb, cpus=cpus, E=evalue)
    else:
        qb, _ = _block(read_seqs(query), "protein")
        method = "phmmer"
        results = pyhmmer.phmmer(qb, tb, cpus=cpus, E=evalue)
    hits = []
    for top in results:
        qn = _s(top.query.name)
        for h in top:
            if h.included:
                hits.append(_mk_hit(h, qn, method + ("(6-frame)" if translated else ""), tseqs, translated))
    return hits


def hmmscan_domains(seqs: dict[str, str], hmm_db: str, evalue: float = 1e-3, cpus: int = 0) -> list[dict]:
    """Domain architecture: scan proteins against an HMM library (e.g. Pfam-A.hmm).
    Returns one row per domain, sorted by position."""
    qb, _ = _block(seqs, "protein")
    hmms = load_hmms(hmm_db)
    rows = []
    for top in pyhmmer.hmmsearch(hmms, qb, cpus=cpus, E=evalue):
        acc = _s(top.query.accession or "") if getattr(top.query, "accession", None) else ""
        for h in top:
            if not h.included:
                continue
            for d in h.domains:
                if d.included:
                    rows.append({"protein": _s(h.name), "domain": _s(top.query.name), "acc": acc,
                                 "start": d.env_from, "end": d.env_to, "evalue": d.i_evalue,
                                 "score": round(d.score, 1)})
    return sorted(rows, key=lambda r: (r["protein"], r["start"]))


# ------------------------------------------------------- local BLAST / DIAMOND
_FMT = "qseqid sseqid evalue bitscore pident length sstart send stitle"


def local_blast(query_fasta: str, db: str, program: str = "blastp", evalue: float = 1e-5,
                threads: int = 4, tool: str = "auto", max_targets: int = 500) -> list[Hit]:
    """Run a BLAST+ program (all 5 directions) or DIAMOND against a local db.
    A FASTA target is turned into a BLAST db automatically."""
    use_diamond = tool == "diamond" or (tool == "auto" and db.endswith(".dmnd"))
    exe = "diamond" if use_diamond else program
    if not shutil.which(exe):
        raise RuntimeError(f"{exe} not installed (conda install -c bioconda blast diamond)")
    with tempfile.TemporaryDirectory() as tmp:
        if not use_diamond and Path(db).is_file() and not Path(db + ".pin").exists() \
                and not Path(db + ".nin").exists():
            dbtype = "nucl" if alphabet_of(db) == "dna" else "prot"
            subprocess.run(["makeblastdb", "-in", db, "-dbtype", dbtype, "-out", f"{tmp}/db"],
                           check=True, capture_output=True)
            db = f"{tmp}/db"
        out = f"{tmp}/out.tsv"
        if use_diamond:
            cmd = ["diamond", program, "-q", query_fasta, "-d", db, "-o", out, "-e", str(evalue),
                   "-p", str(threads), "-k", str(max_targets), "--sensitive", "-f", "6", *_FMT.split()]
        else:
            cmd = [program, "-query", query_fasta, "-db", db, "-out", out, "-evalue", str(evalue),
                   "-num_threads", str(threads), "-max_target_seqs", str(max_targets),
                   "-outfmt", f"6 {_FMT}"]
        subprocess.run(cmd, check=True)
        return _parse_tab(Path(out).read_text(), exe if use_diamond else program)


def _parse_tab(text: str, method: str) -> list[Hit]:
    hits = []
    for line in text.splitlines():
        f = line.split("\t")
        if len(f) < 8:
            continue
        hits.append(Hit(target=f[1], query=f[0], evalue=float(f[2]), score=float(f[3]),
                        method=method, description=f[8] if len(f) > 8 else "",
                        domain_from=int(f[6]), domain_to=int(f[7]),
                        extra={"pident": float(f[4]), "aln_len": int(f[5])}))
    return hits


# ------------------------------------------------------------ remote NCBI BLAST
def remote_blast(queries: dict[str, str], program: str = "blastp", database: str = "nr",
                 evalue: float = 1e-5, hitlist_size: int = 100, entrez_query: str | None = None,
                 email: str | None = None, pause: float = 11.0, log=print) -> list[Hit]:
    """Submit queries one at a time to NCBI BLAST via its URL API (any of the
    5 programs). For more than ~20 queries use a local database instead."""
    from Bio.Blast import NCBIWWW, NCBIXML

    if email:
        NCBIWWW.email = email
    hits: list[Hit] = []
    for i, (qid, seq) in enumerate(queries.items()):
        if i:
            time.sleep(pause)
        log(f"[ncbi] {program} {qid} ({i + 1}/{len(queries)}) vs {database}")
        kw = dict(expect=evalue, hitlist_size=hitlist_size)
        if entrez_query:
            kw["entrez_query"] = entrez_query
        rec = NCBIXML.read(io.StringIO(NCBIWWW.qblast(program, database, seq, **kw).read()))
        for aln in rec.alignments:
            hsp = aln.hsps[0]
            hits.append(Hit(target=aln.accession, query=qid, evalue=hsp.expect, score=hsp.bits,
                            method=f"ncbi-{program}:{database}", description=aln.hit_def[:200],
                            domain_from=hsp.sbjct_start, domain_to=hsp.sbjct_end,
                            sequence=hsp.sbjct.replace("-", "") if program in ("blastp", "blastx") else "",
                            extra={"pident": round(100 * hsp.identities / hsp.align_length, 1)}))
    return hits


def dedupe(hits: Iterable[Hit], key=lambda h: h.target) -> list[Hit]:
    """Keep the best (lowest E-value) hit per key (default: per target)."""
    best: dict = {}
    for h in hits:
        k = key(h)
        if k not in best or h.evalue < best[k].evalue:
            best[k] = h
    return sorted(best.values(), key=lambda h: h.evalue)


def best_hit_per_query(hits: Iterable[Hit]) -> dict[str, Hit]:
    out: dict[str, Hit] = {}
    for h in hits:
        if h.query not in out or h.evalue < out[h.query].evalue:
            out[h.query] = h
    return out
