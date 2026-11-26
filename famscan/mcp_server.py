"""MCP server so Claude (Desktop / Code) can drive famscan as tools.

Register in Claude Code:   claude mcp add famscan -- python -m famscan.mcp_server
"""
from __future__ import annotations

import json
from pathlib import Path

from mcp.server.fastmcp import FastMCP

from . import analyze, features, motifs, search
from .workflows import WORKFLOWS, run_workflow

mcp = FastMCP("famscan")


@mcp.tool()
def list_workflows() -> str:
    """List the use-case workflows (enzyme discovery, genome annotation, drug targets,
    resistance, toxins/allergens, evolution, unknown proteins, lab candidates,
    engineering, strains, family search, BLAST-style search) and their inputs."""
    out = []
    for k, w in WORKFLOWS.items():
        ps = ", ".join(f"{p.name}{'*' if p.required else ''}" for p in w.params)
        out.append(f"{k}: {w.title} - {w.question}\n  inputs: {ps}\n  stages: {' -> '.join(w.stages[:-1])}")
    return "\n".join(out) + "\n(* = required)"


@mcp.tool()
def run(workflow: str, params_json: str, out: str = "famscan_out", formats: str = "md") -> str:
    """Run a workflow. params_json maps input names to values, e.g.
    {"input": "genome.fna", "query": "PF02826", "top_n": 5}. formats: "md" (default),
    or add "xlsx" / "pdf" (e.g. "md,xlsx,pdf" or "all") to also write report.xlsx /
    report.pdf in the output folder. Returns the Markdown report."""
    logs: list[str] = []
    run_workflow(workflow, json.loads(params_json), out, log=logs.append, formats=formats)
    return "\n".join(logs) + "\n\n" + (Path(out) / "report.md").read_text()


@mcp.tool()
def blast_search(query: str, target: str, program: str = "auto", engine: str = "builtin",
                 evalue: float = 1e-5) -> list[dict]:
    """Homology search in any direction (blastp, blastn, blastx, tblastn, tblastx, hmmsearch).
    engine: builtin (offline), blast+, diamond, or ncbi (official API, small batches only)."""
    return [h.row() for h in search.homology_search(query, target, program, engine, evalue)][:200]


@mcp.tool()
def protein_profile(sequence: str) -> dict:
    """ProtParam properties plus predicted secondary structure, TM helices, signal
    peptide, disorder, coiled coils and PROSITE motifs for one protein."""
    r = analyze.protparam(sequence)
    r.update(features.profile(sequence))
    r["motifs"] = motifs.scan(sequence)
    return r


@mcp.tool()
def blosum_compare(sequence: str, reference: str, mode: str = "local") -> dict:
    """BLOSUM62 alignment score, % identity, % similarity and reference coverage."""
    return analyze.blosum_align(sequence, reference, mode=mode).__dict__


if __name__ == "__main__":
    mcp.run()
