"""`famscan start`: ask what the user is doing and route to the right workflow."""
from __future__ import annotations

import re
import shlex
from pathlib import Path

from .workflows import STAGES, WORKFLOWS, Param, run_workflow


def _ask(prompt: str, default=None) -> str:
    sfx = f" [{default}]" if default not in (None, "") else ""
    try:
        v = input(f"{prompt}{sfx}: ").strip()
    except EOFError:
        v = ""
    return v if v else ("" if default is None else str(default))


def _ask_param(p: Param):
    while True:
        extra = f" ({'/'.join(p.choices)})" if p.choices else ""
        tag = "" if p.required else " (Enter to skip)" if p.default is None else ""
        v = _ask(f"  {p.prompt}{extra}{tag}", p.default)
        if not v:
            if p.required:
                print("    required")
                continue
            return None
        if p.kind == "file" and not Path(v).expanduser().exists() and not re.fullmatch(r"PF\d{5}(\.\d+)?", v.upper()):
            print(f"    not found: {v}")
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
        if p.kind == "int":
            return int(v)
        if p.kind == "float":
            return float(v)
        if p.kind == "bool":
            return v.lower() in ("y", "yes", "true", "1")
        return str(Path(v).expanduser()) if p.kind == "file" and Path(v).expanduser().exists() else v


def main(out_default: str = "famscan_out"):
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
        v = _ask_param(p)
        if v is not None:
            params[p.name] = v
    out = _ask("Output folder", out_default)
    cmd = ["famscan", "workflow", wf.key, "--out", out] + [
        x for k, v in params.items() for x in ("--set", f"{k}={' '.join(v) if isinstance(v, list) else v}")]
    print("\nTo rerun without questions:\n  " + " ".join(shlex.quote(c) for c in cmd) + "\n")
    if _ask("Run now? (y/n)", "y").lower().startswith("y"):
        run_workflow(wf.key, params, out)
