"""famscan command line.

  famscan demo                          worked examples on bundled data (start here)
  famscan start                         interactive: pick what you're doing
  famscan list                          show workflows and their inputs
  famscan workflow NAME --set k=v ...   run a workflow without questions
  famscan search -q Q -t T [-p blastx]  one homology search, any direction
  famscan run ...                       classic family pipeline (v0.1 options)
  famscan export DIR                    add Excel/PDF to a finished run
  famscan fetchdb swissprot             download a database to search locally
  famscan config --email you@org         save your EBI contact email on this computer

Add --format md,xlsx,pdf (or all) to any run for Excel and PDF reports.
"""
from __future__ import annotations

import argparse
import json
import os
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


def _terminal_ask():
    """A yes/no prompt when someone is at the keyboard, else None (scripts never block).
    Lets on_service_down=ask work from the command line, not only in `famscan start`."""
    if not (sys.stdin and sys.stdin.isatty()):
        return None

    def ask(question: str) -> bool:
        try:
            return input(f"{question} [Y/n]: ").strip().lower() in ("", "y", "yes")
        except (EOFError, KeyboardInterrupt):
            return False
    return ask


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
    h.add_argument("-e", "--engine", default="builtin",
                   choices=["builtin", "blast+", "diamond", "ncbi", "ebi", "mmseqs"])
    h.add_argument("--email", help="contact email, required by engine ebi only "
                                   "(or save it once: famscan config --email you@example.org)")
    h.add_argument("--on-service-down", default="ask", choices=["ask", "skip", "stop"],
                   help="if the service is down: ask, carry on (skip), or stop")
    h.add_argument("--timeout", type=int, default=900, help="give up on a web service after N seconds")
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

    d = sub.add_parser("demo", help="run a worked example on the bundled data (offline, ~30 s)")
    d.add_argument("--out", default="demo_out")
    d.add_argument("--format", default="md", help="md, xlsx, pdf, comma separated, or all")

    c = sub.add_parser("config", help="save settings on this computer (e.g. your EBI contact email)")
    c.add_argument("--email", help="contact email for the EBI service")
    c.add_argument("--show", action="store_true", help="print the saved settings and where they live")

    f = sub.add_parser("fetchdb", help="download a database so searches run locally")
    f.add_argument("name", nargs="?", help="swissprot | uniref50 | pfam | prosite | trembl | proteome:ACC")
    f.add_argument("--dir", help="where to put it (default ~/.famscan/db or $FAMSCAN_DB)")
    f.add_argument("--force", action="store_true", help="redownload even if present")
    f.add_argument("--list", action="store_true", help="show what is available and what you already have")

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
    from .remote import JobRejected, SearchTimeout, ServiceUnavailable
    try:
        _dispatch(a)
    except (SearchTimeout, ServiceUnavailable, JobRejected) as e:  # web-service trouble: one clear line
        sys.exit(f"famscan: {e}")


def _demo_data() -> Path | None:
    """The bundled demo/data folder, when famscan is run from a clone."""
    for base in (Path(__file__).resolve().parents[1], Path.cwd(), Path.cwd().parent):
        d = base / "demo" / "data"
        if (d / "PF02826.hmm").exists():
            return d
    return None


def _how_invoked(windows: bool | None = None) -> str:
    """How the user started famscan, so printed hints are commands they can paste.
    On Windows `famscan` is often unavailable (no PATH, or blocked by Smart App
    Control) and activation is blocked too, so prefer the famscan.cmd launcher,
    then the running interpreter, and only then the bare `famscan` command."""
    windows = os.name == "nt" if windows is None else windows
    launcher = Path(__file__).resolve().parent.parent / "famscan.cmd"
    if windows and launcher.exists():
        return r".\famscan" if Path.cwd() == launcher.parent else f'& "{launcher}"'
    if Path(sys.argv[0]).stem in ("famscan", "famscan-script"):
        return "famscan"
    exe = Path(sys.executable)
    try:
        exe = exe.relative_to(Path.cwd())
    except ValueError:
        pass
    return f"{exe} -m famscan.cli"


def _demo(out: str, formats: str) -> None:
    """A worked example, so a new user sees real output in half a minute."""
    from .workflows import run_workflow

    data = _demo_data()
    if data is None:
        sys.exit("famscan: the demo data lives in the repository. Get it with:\n"
                 "  git clone https://github.com/timbo-m/Bioinformatic-Pipelines\n"
                 "  cd Bioinformatic-Pipelines && famscan demo")
    runs = [
        ("family", {"input": str(data / "proteome_938293.faa"), "query": str(data / "PF02826.hmm"),
                    "hmm_db": str(data / "mini_domain_library.hmm"), "structure": "none"},
         "Find every member of a protein family (Pfam PF02826) in a 2,100-protein bacterial proteome"),
        ("annotate", {"input": str(data / "bacillus_contig.fna"),
                      "annotation_db": str(data / "reference_proteins.faa"),
                      "hmm_db": str(data / "mini_domain_library.hmm")},
         "Annotate a 391 kb genome contig from scratch: call the genes, then name and describe them"),
        ("engineering", {"reference": str(data / "luxC_reference.faa"),
                         "input": str(data / "luxC_homologs.faa")},
         "Mine 12 natural homologs for mutations that might stabilise your protein"),
    ]
    print("famscan demo - three worked examples on the bundled data. No network needed.\n")
    for i, (key, params, what) in enumerate(runs, 1):
        print(f"[{i}/{len(runs)}] {what}")
        params["formats"] = formats
        folder = Path(out) / f"{i}_{key}"
        ctx = run_workflow(key, params, folder, log=lambda m: print(f"   {m.strip()}"))
        print(f"      -> {folder / 'report.md'} ({len(ctx.focus)} candidates)\n")
    me = _how_invoked()
    print(f"Reports are in {Path(out)}{os.sep}. Now try your own data:\n"
          f"  {me} start          pick what you are trying to do\n"
          f"  {me} list           every workflow and its inputs")


def _dispatch(a):
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
        run_workflow(a.name, params, a.out, ask=_terminal_ask())
    elif a.cmd == "demo":
        _demo(a.out, a.format)
    elif a.cmd == "config":
        from . import config
        if a.email and not config.valid_email(a.email):
            sys.exit(f"famscan: '{a.email}' does not look like an email address")
        if a.email:
            print(f"saved in {config.save(email=a.email)}")
        if a.show or not a.email:
            data = config.load()
            print(f"{config.config_path()}:")
            print("\n".join(f"  {k} = {v}" for k, v in data.items()) if data else "  (empty)")
            print("\nThis file stays on this computer; it is not part of the repository.")
    elif a.cmd == "fetchdb":
        from .databases import describe, fetch, installed
        if a.list or not a.name:
            print("Available:\n" + describe())
            have = installed(a.dir)
            print("\nAlready downloaded:" if have else "\nNothing downloaded yet.")
            for row in have:
                print(f"  {row['path']}  ({row['size_mb']} MB)")
            if not a.name:
                return
        path = fetch(a.name, a.dir, force=a.force, make_blast_db=True)
        print(f"\nUse it with:\n  famscan search -q query.faa -t {path}\n"
              f"  famscan workflow annotate --set input=genome.fna --set annotation_db={path}")
    elif a.cmd == "export":
        from .export import export
        if not (Path(a.dir) / "report.md").exists():
            sys.exit(f"famscan: {a.dir}/report.md not found - give the --out folder of an earlier run")
        export(a.dir, a.format)
    elif a.cmd == "search":
        from .workflows import run_workflow
        run_workflow("search", {"query": a.query, "target": a.target, "program": a.program,
                                "engine": a.engine, "evalue": a.evalue, "email": a.email,
                                "timeout": a.timeout, "on_service_down": a.on_service_down,
                                "formats": a.format}, a.out, ask=_terminal_ask())
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
