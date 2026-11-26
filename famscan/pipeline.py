"""End-to-end pipeline: search -> property/BLOSUM filter -> structure -> report."""
from __future__ import annotations

import json
import time
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd
import pyhmmer

from . import analyze, search, structure

DEFAULT_RULES = {
    "length": [80, 2000],
    "ref_coverage": [30, None],
    "instability": [None, 60],
}


@dataclass
class Config:
    proteome: str | None = None          # FASTA to search (local mode)
    query: str | None = None             # HMM file, Pfam acc, or query FASTA
    blast_db: str | None = None          # local BLAST/DIAMOND db (instead of proteome)
    remote_db: str | None = None         # e.g. "swissprot" -> NCBI URL API
    references: str | None = None        # FASTA of known family members
    reference_pdb: str | None = None     # known structure for TM-score
    evalue: float = 1e-5
    rules: dict = field(default_factory=lambda: dict(DEFAULT_RULES))
    top_n: int = 10
    structure_method: str = "auto"       # auto | esmfold | afdb | none
    pdb_dir: str | None = None
    foldseek_db: str | None = None
    foldseek_web: bool = False
    plddt_min: float = 70.0
    tm_min: float = 0.5
    cpus: int = 0
    out: str = "famscan_out"


def run(cfg: Config, log=print) -> pd.DataFrame:
    out = Path(cfg.out)
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    timings = {}

    # ---------------------------------------------------------------- search
    log("[1/4] searching")
    if cfg.remote_db:
        q = search.read_fasta(cfg.query)
        hits = search.remote_blast(q, database=cfg.remote_db, evalue=cfg.evalue, log=log)
    elif cfg.blast_db:
        hits = search.local_blast(cfg.query, cfg.blast_db, evalue=cfg.evalue)
        if cfg.proteome:
            seqs = search.read_fasta(cfg.proteome)
            for h in hits:
                h.sequence = seqs.get(h.target, "")
    else:
        hits = search.hmmer_search(cfg.query, cfg.proteome, evalue=cfg.evalue,
                                   cpus=cfg.cpus, workdir=out)
    hits = search.dedupe(hits)
    timings["search_s"] = round(time.time() - t0, 2)
    log(f"      {len(hits)} unique hits")
    if not hits:
        pd.DataFrame().to_csv(out / "results.tsv", sep="\t")
        return pd.DataFrame()

    # ------------------------------------------------------------- analyse
    t1 = time.time()
    log("[2/4] ProtParam + BLOSUM62 scoring")
    refs = search.read_fasta(cfg.references) if cfg.references else _hmm_consensus(cfg.query)
    rows = []
    for h in hits:
        row = h.row()
        if h.sequence:
            row.update(analyze.protparam(h.sequence))
            if refs:
                row.update(analyze.best_reference_match(h.sequence, refs))
        ok, why = analyze.passes(row, cfg.rules)
        row["pass_filter"], row["fail_reason"] = ok, why
        rows.append(row)
    df = pd.DataFrame(rows)
    timings["analyse_s"] = round(time.time() - t1, 2)
    log(f"      {int(df.pass_filter.sum())}/{len(df)} pass filters {cfg.rules}")

    # ----------------------------------------------------------- structure
    t2 = time.time()
    seqmap = {h.target: h.sequence for h in hits}
    top = df[df.pass_filter].sort_values("score", ascending=False).head(cfg.top_n)
    for col in ("structure_source", "plddt", "tm_score", "rmsd", "pdb", "foldseek_top"):
        df[col] = None
    if cfg.structure_method != "none" and len(top):
        log(f"[3/4] structures for top {len(top)} hits ({cfg.structure_method})")
        pdbs = []
        for idx, r in top.iterrows():
            pdb, src = structure.get_structure(r.target, seqmap[r.target], out / "structures",
                                               Path(cfg.pdb_dir) if cfg.pdb_dir else None,
                                               cfg.structure_method, log=log)
            df.at[idx, "structure_source"] = src
            if pdb is None:
                continue
            pdbs.append(pdb)
            df.at[idx, "pdb"] = str(pdb)
            df.at[idx, "plddt"] = structure.plddt(pdb)
            if cfg.reference_pdb:
                tm = structure.tm_score(pdb, Path(cfg.reference_pdb))
                df.at[idx, "tm_score"], df.at[idx, "rmsd"] = tm["tm_score"], tm["rmsd"]
            if cfg.foldseek_web:
                try:
                    fs = structure.foldseek_web(pdb)
                    df.at[idx, "foldseek_top"] = fs[0]["target"] if fs else None
                except Exception as e:  # network / server errors are non-fatal
                    log(f"      foldseek web failed for {r.target}: {e}")
        if cfg.foldseek_db and pdbs:
            structure.foldseek_local(pdbs, cfg.foldseek_db, out / "foldseek.m8")
    else:
        log("[3/4] structure stage skipped")
    timings["structure_s"] = round(time.time() - t2, 2)

    # --------------------------------------------------------------- verdict
    df["verdict"] = df.apply(lambda r: _verdict(r, cfg), axis=1)
    log("[4/4] writing report")
    df.to_csv(out / "results.tsv", sep="\t", index=False)
    _write_fasta(df[df.pass_filter], seqmap, out / "filtered_hits.faa")
    _write_report(df, cfg, timings, out / "report.md")
    (out / "run.json").write_text(json.dumps({"config": cfg.__dict__, "timings": timings}, indent=2, default=str))
    log(f"done in {time.time() - t0:.1f}s -> {out}/report.md")
    return df


def _hmm_consensus(query: str | None) -> dict[str, str]:
    """If no reference FASTA is given, use the HMM consensus as reference."""
    if not query or not Path(query).exists():
        return {}
    try:
        with pyhmmer.plan7.HMMFile(query) as hf:
            return {f"{search._s(h.name)}_consensus": h.consensus.upper() for h in hf}
    except Exception:
        return {}


def _verdict(r, cfg: Config) -> str:
    if not r.get("pass_filter"):
        return "filtered"
    if r.get("plddt") is None or pd.isna(r.get("plddt")):
        return "candidate (no structure yet)"
    if r["plddt"] < cfg.plddt_min:
        return "low-confidence fold"
    tm = r.get("tm_score")
    if tm is not None and not pd.isna(tm):
        return "confirmed" if tm >= cfg.tm_min else "different fold"
    return "well-folded (no reference to compare)"


def _write_fasta(df: pd.DataFrame, seqs: dict, path: Path) -> None:
    with open(path, "w") as fh:
        for t in df.target:
            if seqs.get(t):
                fh.write(f">{t}\n{seqs[t]}\n")


def _write_report(df: pd.DataFrame, cfg: Config, timings: dict, path: Path) -> None:
    cols = [c for c in ["target", "evalue", "score", "length", "mw_kda", "pI", "instability",
                        "gravy", "ref_identity", "ref_coverage", "plddt", "tm_score", "verdict"]
            if c in df.columns]
    view = df.sort_values("score", ascending=False)[cols].head(30).copy()
    view["evalue"] = view["evalue"].map(lambda x: f"{x:.1e}")
    view["score"] = view["score"].round(1)
    lines = [
        "# famscan report", "",
        f"- Query: `{cfg.query}`",
        f"- Searched: `{cfg.proteome or cfg.blast_db or cfg.remote_db}`",
        f"- Hits: {len(df)}; passing filters: {int(df.pass_filter.sum())}",
        f"- Filters: `{cfg.rules}`",
        f"- Timings (s): {timings}", "",
        "## Verdicts", "",
        *[f"- {k}: {v}" for k, v in df.verdict.value_counts().items()], "",
        "## Top hits", "",
        _md_table(view), "",
        "Full table: `results.tsv`. Sequences passing filters: `filtered_hits.faa`.",
    ]
    path.write_text("\n".join(lines))


def _md_table(df: pd.DataFrame) -> str:
    head = "| " + " | ".join(df.columns) + " |"
    sep = "|" + "---|" * len(df.columns)
    body = ["| " + " | ".join("" if pd.isna(v) else str(v) for v in row) + " |"
            for row in df.itertuples(index=False)]
    return "\n".join([head, sep, *body])
