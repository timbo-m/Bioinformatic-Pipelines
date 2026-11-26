"""famscan command line.

  famscan start                         interactive: pick what you're doing
  famscan list                          show workflows and their inputs
  famscan workflow NAME --set k=v ...   run a workflow without questions
  famscan search -q Q -t T [-p blastx]  one homology search, any direction
  famscan run ...                       classic family pipeline (v0.1 options)
  famscan export DIR                    add Excel/PDF to a finished run

Add --format md,xlsx,pdf (or all) to any run for Excel and PDF reports.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .pipeline import Config, DEFAULT_RULES, run


def _coerce(v: str):
    for t in (int, float):
        try:
            return t(v)
        except ValueError:
            pass
    if v.lower() in ("true", "false"):
        return v.lower() == "true"
    return v


def main(argv=None):
    p = argparse.ArgumentParser(prog="famscan", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("start", help="interactive menu")
    s.add_argument("--out", default="famscan_out")
    sub.add_parser("list", help="list workflows")
    w = sub.add_parser("workflow", help="run a workflow non-interactively")
    w.add_argument("name")
    w.add_argument("--set", action="append", default=[], metavar="KEY=VALUE")
    w.add_argument("--out", default="famscan_out")
    w.add_argument("--format", default="md", help="md, xlsx, pdf, comma separated, or all")
    h = sub.add_parser("search", help="homology search in any BLAST direction")
    h.add_argument("-q", "--query", required=True)
    h.add_argument("-t", "--target", required=True)
    h.add_argument("-p", "--program", default="auto",
                   choices=["auto", "blastp", "blastn", "blastx", "tblastn", "tblastx", "hmmsearch"])
    h.add_argument("-e", "--engine", default="builtin", choices=["builtin", "blast+", "diamond", "ncbi"])
    h.add_argument("--evalue", type=float, default=1e-5)
    h.add_argument("--out", default="famscan_out")
    h.add_argument("--format", default="md", help="md, xlsx, pdf, comma separated, or all")

    r = sub.add_parser("run", help="classic family pipeline")
    r.add_argument("--query", required=True, help="HMM file, Pfam accession or query FASTA")
    src = r.add_mutually_exclusive_group(required=True)
    src.add_argument("--proteome")
    src.add_argument("--blast-db")
    src.add_argument("--remote-db")
    r.add_argument("--seqs")
    r.add_argument("--references")
    r.add_argument("--reference-pdb")
    r.add_argument("--evalue", type=float, default=1e-5)
    r.add_argument("--rules", help=f"JSON, default {json.dumps(DEFAULT_RULES)}")
    r.add_argument("--top", type=int, default=10)
    r.add_argument("--structure", default="auto", choices=["auto", "esmfold", "afdb", "none"])
    r.add_argument("--pdb-dir")
    r.add_argument("--foldseek-db")
    r.add_argument("--foldseek-web", action="store_true")
    r.add_argument("--plddt-min", type=float, default=70)
    r.add_argument("--tm-min", type=float, default=0.5)
    r.add_argument("--cpus", type=int, default=0)
    r.add_argument("--out", default="famscan_out")
    r.add_argument("--format", default="md", help="md, xlsx, pdf, comma separated, or all")

    x = sub.add_parser("export", help="write report.xlsx / report.pdf for a finished run")
    x.add_argument("dir", help="output folder of an earlier run (must contain report.md)")
    x.add_argument("--format", default="xlsx,pdf", help="md, xlsx, pdf, comma separated, or all")
    a = p.parse_args(argv)
    if getattr(a, "format", None):
        from .export import parse_formats
        try:
            parse_formats(a.format)
        except ValueError as e:
            sys.exit(f"famscan: {e}")

    if a.cmd == "start":
        from .interactive import main as interactive
        interactive(a.out)
    elif a.cmd == "list":
        from .workflows import WORKFLOWS
        for k, wf in WORKFLOWS.items():
            print(f"{k:15s} {wf.title}\n{'':15s} {wf.question}\n{'':15s} stages: {' -> '.join(wf.stages[:-1])}")
            for prm in wf.params:
                print(f"{'':17s}{prm.name}{'*' if prm.required else ''}: {prm.prompt}")
            print()
    elif a.cmd == "workflow":
        from .workflows import WORKFLOWS, run_workflow
        if a.name not in WORKFLOWS:
            sys.exit(f"unknown workflow; choose from {', '.join(WORKFLOWS)}")
        params = {}
        for kv in a.set:
            k, v = kv.split("=", 1)
            params[k] = v.split() if k == "proteomes" else _coerce(v)
        params.setdefault("formats", a.format)
        run_workflow(a.name, params, a.out)
    elif a.cmd == "export":
        from .export import export
        if not (Path(a.dir) / "report.md").exists():
            sys.exit(f"famscan: {a.dir}/report.md not found - give the --out folder of an earlier run")
        export(a.dir, a.format)
    elif a.cmd == "search":
        from .workflows import run_workflow
        run_workflow("search", {"query": a.query, "target": a.target, "program": a.program,
                                "engine": a.engine, "evalue": a.evalue, "formats": a.format}, a.out)
    else:
        cfg = Config(query=a.query, proteome=a.proteome or a.seqs, blast_db=a.blast_db,
                     remote_db=a.remote_db, references=a.references, reference_pdb=a.reference_pdb,
                     evalue=a.evalue, rules=json.loads(a.rules) if a.rules else dict(DEFAULT_RULES),
                     top_n=a.top, structure_method=a.structure, pdb_dir=a.pdb_dir,
                     foldseek_db=a.foldseek_db, foldseek_web=a.foldseek_web,
                     plddt_min=a.plddt_min, tm_min=a.tm_min, cpus=a.cpus, out=a.out)
        run(cfg)
        from .export import export
        export(a.out, a.format)


if __name__ == "__main__":
    main()
