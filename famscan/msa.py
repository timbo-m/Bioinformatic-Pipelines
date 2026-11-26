"""Alignment-level analyses: multiple sequence alignment, conservation,
phylogenetic trees, redundancy clustering and variant mining for engineering."""
from __future__ import annotations

import io
import math
from collections import Counter
from pathlib import Path

from Bio import Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .analyze import blosum_align


# ------------------------------------------------------------------- MSA
def align(seqs: dict[str, str]) -> dict[str, str]:
    """Multiple sequence alignment with FAMSA (pyfamsa). Returns {id: gapped}."""
    import pyfamsa

    recs = [pyfamsa.Sequence(k.encode(), v.encode()) for k, v in seqs.items()]
    aln = pyfamsa.Aligner(guide_tree="upgma").align(recs)
    return {s.id.decode(): s.sequence.decode() for s in aln}


def write_aln(aln: dict[str, str], path: Path) -> Path:
    with open(path, "w") as fh:
        for k, v in aln.items():
            fh.write(f">{k}\n{v}\n")
    return path


def conservation(aln: dict[str, str]) -> list[dict]:
    """Per-column Shannon-entropy conservation (1 = invariant, 0 = random)
    plus consensus residue and gap fraction."""
    rows = list(aln.values())
    ncol = len(rows[0]) if rows else 0
    out = []
    for j in range(ncol):
        col = [r[j] for r in rows]
        res = [c for c in col if c != "-"]
        gap = 1 - len(res) / len(col)
        if not res:
            out.append({"col": j + 1, "consensus": "-", "conservation": 0.0, "gap_frac": 1.0})
            continue
        cnt = Counter(res)
        h = -sum(v / len(res) * math.log2(v / len(res)) for v in cnt.values())
        cons = (1 - h / math.log2(20)) * (1 - gap)
        top, freq = cnt.most_common(1)[0]
        out.append({"col": j + 1, "consensus": top, "consensus_freq": round(freq / len(res), 2),
                    "conservation": round(cons, 3), "gap_frac": round(gap, 2)})
    return out


def conserved_positions(aln: dict[str, str], ref_id: str, min_cons: float = 0.8) -> list[dict]:
    """Highly conserved columns mapped to residue numbers of `ref_id` -
    candidates for active/binding sites."""
    cons = conservation(aln)
    ref = aln[ref_id]
    pos, out = 0, []
    for c, a in zip(cons, ref):
        if a != "-":
            pos += 1
            if c["conservation"] >= min_cons:
                out.append({"residue": f"{a}{pos}", "conservation": c["conservation"]})
    return out


# ------------------------------------------------------------------ trees
def tree(aln: dict[str, str], method: str = "nj", model: str = "blosum62"):
    """Distance tree (NJ or UPGMA) with Biopython. For ML trees use IQ-TREE
    or FastTree on the written alignment."""
    msa = MultipleSeqAlignment([SeqRecord(Seq(v), id=k[:60]) for k, v in aln.items()])
    dm = DistanceCalculator(model).get_distance(msa)
    t = DistanceTreeConstructor().nj(dm) if method == "nj" else DistanceTreeConstructor().upgma(dm)
    t.root_at_midpoint()
    return t


def tree_to_newick(t) -> str:
    buf = io.StringIO()
    Phylo.write(t, buf, "newick")
    return buf.getvalue().strip()


def tree_ascii(t) -> str:
    buf = io.StringIO()
    Phylo.draw_ascii(t, file=buf, column_width=70)
    return buf.getvalue()


# ------------------------------------------------------------- clustering
def _kmers(s: str, k: int = 3) -> set:
    return {s[i:i + k] for i in range(len(s) - k + 1)}


def cluster(seqs: dict[str, str], min_identity: float = 90.0, min_cov: float = 80.0,
            k: int = 3) -> dict[str, list[str]]:
    """Greedy CD-HIT-style clustering: longest first; a sequence joins the first
    representative it matches at >= min_identity over >= min_cov of its length.
    k-mer Jaccard prefilter keeps it fast. For >10k sequences use MMseqs2."""
    order = sorted(seqs, key=lambda x: -len(seqs[x]))
    reps: dict[str, list[str]] = {}
    kmer = {i: _kmers(seqs[i], k) for i in order}
    for sid in order:
        placed = False
        for rep in reps:
            a, b = kmer[sid], kmer[rep]
            if len(a & b) / max(len(a | b), 1) < (min_identity / 100) ** 6 * 0.5:
                continue
            r = blosum_align(seqs[rep], seqs[sid], mode="local")  # coverage of the shorter seq
            if r.identity >= min_identity and r.coverage >= min_cov:
                reps[rep].append(sid)
                placed = True
                break
        if not placed:
            reps[sid] = [sid]
    return reps


# --------------------------------------------------------- variant mining
def variant_table(aln: dict[str, str], ref_id: str, min_freq: float = 0.2) -> list[dict]:
    """Natural variation at each position of the reference: which substitutions
    occur in homologs and how often. Positions where the reference differs
    from a strong consensus are 'back-to-consensus' candidates, a classic
    stabilising-mutation strategy."""
    ref = aln[ref_id]
    others = [v for k, v in aln.items() if k != ref_id]
    cons = conservation(aln)
    rows, pos = [], 0
    for j, a in enumerate(ref):
        if a == "-":
            continue
        pos += 1
        col = [o[j] for o in others if o[j] != "-"]
        if not col:
            continue
        cnt = Counter(col)
        top, n = cnt.most_common(1)[0]
        alts = {aa: round(c / len(col), 2) for aa, c in cnt.most_common() if aa != a and c / len(col) >= min_freq}
        if not alts:
            continue
        rows.append({"position": pos, "ref": a, "alternatives": ",".join(f"{k}:{v}" for k, v in alts.items()),
                     "consensus": top, "consensus_freq": round(n / len(col), 2),
                     "conservation": cons[j]["conservation"],
                     "back_to_consensus": top != a and n / len(col) >= 0.6,
                     "suggested": f"{a}{pos}{top}" if top != a and n / len(col) >= 0.6 else ""})
    return rows
