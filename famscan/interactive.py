"""`famscan start`: ask what the user is doing and route to the right workflow."""
from __future__ import annotations

import re
import shlex
from pathlib import Path

from . import config
from .databases import CATALOG, db_dir, fetch, installed
from .workflows import STAGES, WORKFLOWS, Param, run_workflow

# Which downloadable database each input can be satisfied by.
DB_FOR_PARAM = {
    "annotation_db": "swissprot",
    "hmm_db": "pfam",
    "prosite_dat": "prosite",
    "references": "swissprot",
    "target": "swissprot",
    "host": "proteome:UP000005640",
    "input": "swissprot",
}


def _local_match(name: str) -> str | None:
    """A copy of this database already downloaded?"""
    want = CATALOG.get(name, {}).get("filename") or (name.split(":", 1)[-1] + ".fasta")
    for row in installed():
        if row["name"] == want:
            return row["path"]
    return None


def _offer_download(p: Param) -> str | None:
    """When an input needs a database the user doesn't have, offer to fetch it."""
    name = DB_FOR_PARAM.get(p.name)
    if not name:
        return None
    have = _local_match(name)
    if have:
        if _ask(f"    use the {name} you already have at {have}? (y/n)", "y").lower().startswith("y"):
            return have
        return None
    entry = CATALOG.get(name, {})
    size = entry.get("size", "~5-30 MB")
    about = entry.get("about", f"UniProt proteome {name.split(':')[-1]}")
    print(f"    {name}: {about}")
    print(f"    Downloading it once ({size}) means local searches: no queue, no limits, works offline.")
    if not _ask(f"    download {name} to {db_dir()} now? (y/n)", "n").lower().startswith("y"):
        return None
    try:
        return str(fetch(name, make_blast_db=True))
    except Exception as e:  # noqa: BLE001 - network/disk problems are not fatal here
        print(f"    download failed: {str(e)[:150]}")
        return None


class _Cancelled(Exception):
    """stdin closed or the user pressed Ctrl-C."""


def _ask(prompt: str, default=None) -> str:
    sfx = f" [{default}]" if default not in (None, "") else ""
    try:
        v = input(f"{prompt}{sfx}: ").strip()
    except (EOFError, KeyboardInterrupt):
        raise _Cancelled from None
    return v if v else ("" if default is None else str(default))


def _confirm(question: str) -> bool:
    """Yes/no prompt used mid-run when a web service turns out to be down."""
    try:
        return _ask(f"{question} (y/n)", "y").lower().startswith("y")
    except _Cancelled:
        return False


def _ask_email(saved: str | None) -> str | None:
    """The EBI requires a contact address with every job."""
    print("    The EBI requires a contact email with each job; they use it only to reach you about "
          "misuse. It is sent to the EBI and saved on this computer only.")
    while True:
        v = _ask("    Your email", saved or "")
        if not v:
            print("    needed for the EBI service - pick another engine to skip it")
            return None
        if config.valid_email(v):
            if v != saved and _ask(f"    remember {v} for next time? (y/n)", "y").lower().startswith("y"):
                print(f"    saved in {config.save(email=v)}")
            return v
        print("    that does not look like an email address")


def _ask_param(p: Param):
    while True:
        extra = f" ({'/'.join(p.choices)})" if p.choices else ""
        tag = "" if p.required else " (Enter to skip)" if p.default is None else ""
        v = _ask(f"  {p.prompt}{extra}{tag}", p.default)
        if not v:
            got = _offer_download(p) if p.kind == "file" else None
            if got:
                return got
            if p.required:
                print("    required")
                continue
            return None
        if p.kind == "file" and not Path(v).expanduser().exists() and not re.fullmatch(r"PF\d{5}(\.\d+)?", v.upper()):
            print(f"    not found: {v}")
            got = _offer_download(p)
            if got:
                return got
            continue
        if p.kind == "files":
            parts = shlex.split(v)
            bad = [x for x in parts if not Path(x).expanduser().exists()]
            if bad:
                print(f"    not found: {' '.join(bad)}")
                continue
            return parts
        if p.kind == "choice" and p.choices and v not in p.choices:
            print(f"    choose one of {', '.join(p.choices)}")
            continue
        if p.kind in ("int", "float"):
            try:
                return int(v) if p.kind == "int" else float(v)
            except ValueError:
                print(f"    {p.kind} expected, got {v!r}")
                continue
        if p.kind == "bool":
            return v.lower() in ("y", "yes", "true", "1")
        return str(Path(v).expanduser()) if p.kind == "file" and Path(v).expanduser().exists() else v


def main(out_default: str = "famscan_out"):
    try:
        _menu(out_default)
    except _Cancelled:
        print("\ncancelled - nothing was run. Use `famscan workflow ...` to run without questions.")


def _menu(out_default: str = "famscan_out"):
    keys = list(WORKFLOWS)
    print("\nfamscan - what are you trying to do?\n")
    for i, k in enumerate(keys, 1):
        w = WORKFLOWS[k]
        print(f"  {i:2d}. {w.title:<42} {w.question}")
    while True:
        c = _ask("\nPick a number")
        if c.isdigit() and 1 <= int(c) <= len(keys):
            wf = WORKFLOWS[keys[int(c) - 1]]
            break
        print("  enter a number from the list")
    print(f"\n{wf.title}\nPipeline: " + " -> ".join(wf.stages[:-1]) + "\n")
    for s in wf.stages[:-1]:
        print(f"   - {s}: {STAGES[s][1]}")
    if wf.notes:
        print(f"\n  Note: {wf.notes}")
    print("\nInputs (files can be protein or DNA; genomes get gene calling automatically):")
    params = {}
    for p in wf.params:
        if p.name == "email":
            if params.get("engine") != "ebi":
                continue
            v = _ask_email(config.email())
        else:
            v = _ask_param(p)
        if v is not None:
            params[p.name] = v
    if params.get("engine") == "ebi" and not params.get("email"):
        print("  no email given, so switching the engine back to builtin (no network needed)")
        params["engine"] = "builtin"
    out = _ask("Output folder", out_default)
    while True:
        fmt = _ask("Report formats: md, xlsx (Excel), pdf - comma separated, or all", "md")
        try:
            from .export import parse_formats
            params["formats"] = ",".join(parse_formats(fmt))
            break
        except ValueError as e:
            print(f"    {e}")
    cmd = ["famscan", "workflow", wf.key, "--out", out] + [
        x for k, v in params.items() for x in ("--set", f"{k}={' '.join(v) if isinstance(v, list) else v}")]
    print("\nTo rerun without questions:\n  " + " ".join(shlex.quote(c) for c in cmd) + "\n")
    if _ask("Run now? (y/n)", "y").lower().startswith("y"):
        run_workflow(wf.key, params, out, ask=_confirm)
