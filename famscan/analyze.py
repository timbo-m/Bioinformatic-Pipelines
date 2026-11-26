"""Stage 2: fast local filtering - ProtParam properties and BLOSUM scoring.

Everything here is pure Python / Biopython and runs in milliseconds per
sequence, so it is applied to *every* hit before the expensive structure
stage.
"""
from __future__ import annotations

from dataclasses import dataclass

from Bio import Align
from Bio.Align import substitution_matrices
from Bio.SeqUtils.ProtParam import ProteinAnalysis

_STD = set("ACDEFGHIKLMNPQRSTVWY")


def protparam(seq: str) -> dict:
    """ExPASy ProtParam-style properties."""
    clean = "".join(c for c in seq.upper() if c in _STD)
    if not clean:
        return {}
    pa = ProteinAnalysis(clean)
    helix, turn, sheet = pa.secondary_structure_fraction()
    ext_red, ext_ox = pa.molar_extinction_coefficient()
    return {
        "length": len(seq),
        "mw_kda": round(pa.molecular_weight() / 1000, 2),
        "pI": round(pa.isoelectric_point(), 2),
        "instability": round(pa.instability_index(), 1),
        "stable": pa.instability_index() < 40,
        "gravy": round(pa.gravy(), 3),
        "aromaticity": round(pa.aromaticity(), 3),
        "charge_pH7": round(pa.charge_at_pH(7.0), 1),
        "frac_helix": round(helix, 3),
        "frac_sheet": round(sheet, 3),
        "ext_coeff": ext_ox,
    }


@dataclass
class AlnResult:
    score: float
    identity: float      # % identical over aligned columns
    similarity: float    # % positive BLOSUM pairs over aligned columns
    coverage: float      # % of reference covered
    aln_len: int = 0     # aligned columns (no gaps)


def _aligner(matrix: str = "BLOSUM62", mode: str = "local") -> Align.PairwiseAligner:
    a = Align.PairwiseAligner()
    a.substitution_matrix = substitution_matrices.load(matrix)
    a.open_gap_score = -11
    a.extend_gap_score = -1
    a.mode = mode
    return a


def blosum_align(query: str, reference: str, matrix: str = "BLOSUM62",
                 mode: str = "local") -> AlnResult:
    aligner = _aligner(matrix, mode)
    aln = aligner.align(query, reference)[0]
    m = aligner.substitution_matrix
    ident = pos = cols = 0
    ref_cov = 0
    for (qs, qe), (rs, re_) in zip(*aln.aligned):
        ref_cov += re_ - rs
        for a, b in zip(query[qs:qe], reference[rs:re_]):
            cols += 1
            ident += a == b
            try:
                pos += m[a, b] > 0
            except (KeyError, IndexError):
                pass
    cols = cols or 1
    return AlnResult(score=float(aln.score), identity=round(100 * ident / cols, 1),
                     similarity=round(100 * pos / cols, 1),
                     coverage=round(100 * ref_cov / max(len(reference), 1), 1), aln_len=cols)


def best_reference_match(seq: str, references: dict[str, str], **kw) -> dict:
    """Align `seq` to every reference; return the best by score."""
    best_id, best = None, None
    for rid, rseq in references.items():
        r = blosum_align(seq, rseq, **kw)
        if best is None or r.score > best.score:
            best_id, best = rid, r
    if best is None:
        return {}
    return {"ref": best_id, "blosum_score": best.score, "ref_identity": best.identity,
            "ref_similarity": best.similarity, "ref_coverage": best.coverage}


def passes(row: dict, rules: dict) -> tuple[bool, str]:
    """Apply simple min/max rules, e.g. {"length": [150, 600], "pI": [None, 8]}."""
    for key, (lo, hi) in rules.items():
        v = row.get(key)
        if v is None:
            continue
        if lo is not None and v < lo:
            return False, f"{key}={v}<{lo}"
        if hi is not None and v > hi:
            return False, f"{key}={v}>{hi}"
    return True, ""
