"""Sequence input: read FASTA/GenBank, detect DNA vs protein, translate, call genes."""
from __future__ import annotations

from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

_DNA = set("ACGTUN")


def read_seqs(path: str | Path) -> dict[str, str]:
    """Read FASTA or GenBank into {id: sequence} (upper case, no stop '*')."""
    path = str(path)
    fmt = "genbank" if path.endswith((".gb", ".gbk", ".gbff")) else "fasta"
    return {r.id: str(r.seq).upper().rstrip("*") for r in SeqIO.parse(path, fmt)}


def read_descriptions(path: str | Path) -> dict[str, str]:
    """{id: header text after the id} - used to name proteins by homology."""
    path = str(path)
    fmt = "genbank" if path.endswith((".gb", ".gbk", ".gbff")) else "fasta"
    return {r.id: r.description[len(r.id):].strip() for r in SeqIO.parse(path, fmt)}


def write_fasta(seqs: dict[str, str], path: str | Path) -> Path:
    path = Path(path)
    with open(path, "w") as fh:
        for k, v in seqs.items():
            fh.write(f">{k}\n")
            for i in range(0, len(v), 60):
                fh.write(v[i:i + 60] + "\n")
    return path


def alphabet_of(seqs: dict[str, str] | str) -> str:
    """'dna' or 'protein' from the first ~10 kb of sequence."""
    if isinstance(seqs, (str, Path)) and Path(seqs).exists():
        seqs = read_seqs(seqs)
    sample = "".join(list(seqs.values()))[:10000] if isinstance(seqs, dict) else seqs[:10000]
    if not sample:
        return "protein"
    frac = sum(c in _DNA for c in sample) / len(sample)
    return "dna" if frac > 0.9 else "protein"


def six_frame(seqs: dict[str, str], min_len: int = 30) -> dict[str, str]:
    """Translate every frame and split at stops: ORF-like peptides >= min_len aa.

    IDs encode the source: <id>|f<+1..-3>|<nt_start>-<nt_end> (1-based, forward coords).
    """
    out = {}
    for sid, nt in seqs.items():
        nt = nt.replace("U", "T")
        n = len(nt)
        for strand, s in ((1, nt), (-1, str(Seq(nt).reverse_complement()))):
            for frame in range(3):
                sub = s[frame:]
                sub = sub[: len(sub) - len(sub) % 3]
                prot = str(Seq(sub).translate())
                pos = 0
                for pep in prot.split("*"):
                    if len(pep) >= min_len:
                        a = frame + 3 * pos
                        b = a + 3 * len(pep)
                        start, end = (a + 1, b) if strand == 1 else (n - b + 1, n - a)
                        out[f"{sid}|f{strand * (frame + 1):+d}|{start}-{end}"] = pep
                    pos += len(pep) + 1
    return out


def call_genes(seqs: dict[str, str], meta: bool | None = None) -> tuple[dict[str, str], dict[str, str], list[dict]]:
    """Prokaryotic gene calling with pyrodigal (Prodigal). Returns (proteins, cds_nt, table).

    For a single genome >= 100 kb the model is trained on it; otherwise metagenomic
    mode is used. Eukaryotic genomes need a splice-aware tool (AUGUSTUS, BRAKER) -
    pass their predicted proteome instead.
    """
    import pyrodigal

    total = sum(len(s) for s in seqs.values())
    meta = (total < 100_000) if meta is None else meta
    finder = pyrodigal.GeneFinder(meta=meta)
    if not meta:
        finder.train(*[s.encode() for s in seqs.values()])
    prots, cds, table = {}, {}, []
    for sid, nt in seqs.items():
        for i, g in enumerate(finder.find_genes(nt.encode()), 1):
            gid = f"{sid}_{i}"
            prots[gid] = g.translate().rstrip("*")
            cds[gid] = g.sequence()
            table.append({"gene": gid, "contig": sid, "start": g.begin, "end": g.end,
                          "strand": g.strand, "partial": g.partial_begin or g.partial_end,
                          "gc": round(g.gc_cont, 3), "rbs": g.rbs_motif})
    return prots, cds, table


def to_proteins(path: str | Path, log=print) -> tuple[dict[str, str], str]:
    """Load any input as proteins: protein FASTA as is; DNA via gene calling
    (long contigs) or direct translation (short CDS-like records)."""
    seqs = read_seqs(path)
    if alphabet_of(seqs) == "protein":
        return seqs, "protein"
    if max(len(s) for s in seqs.values()) > 20_000:
        prots, _, _ = call_genes(seqs)
        log(f"      called {len(prots)} genes with Prodigal")
        return prots, "genome"
    prots = {k: str(Seq(v[: len(v) - len(v) % 3]).translate()).rstrip("*").split("*")[0]
             for k, v in seqs.items()}
    log(f"      translated {len(prots)} coding sequences")
    return prots, "cds"
