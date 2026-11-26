"""Motifs: PROSITE-pattern scanning, user regex motifs and de novo motif discovery.

* A small built-in set of well-known PROSITE patterns works offline.
* Download the full library (https://ftp.expasy.org/databases/prosite/prosite.dat)
  and pass it as `prosite_dat` to scan all ~1,300 patterns.
* De novo discovery is a fast enriched-k-mer search. If MEME/STREME is installed
  (`streme` on PATH) `denovo_streme` runs the real thing.
"""
from __future__ import annotations

import itertools
import random
import re
import shutil
import subprocess
import tempfile
from collections import Counter
from pathlib import Path

# id, name, PROSITE pattern. Very common "post-translational" patterns
# (PS00001-PS00008) match by chance in most proteins; they are flagged `frequent`.
BUILTIN = [
    ("PS00017", "ATP/GTP-binding P-loop (Walker A)", "[AG]-x(4)-G-K-[ST]", False),
    ("PS00028", "Zinc finger C2H2", "C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H", False),
    ("PS00018", "EF-hand calcium-binding", "D-{W}-[DNS]-{ILVFYW}-[DENSTG]-[DNQGHRK]-{GP}-[LIVMC]-[DENQSTAGC]-x(2)-[DE]-[LIVMFYW]", False),
    ("PS00029", "Leucine zipper", "L-x(6)-L-x(6)-L-x(6)-L", False),
    ("PS00016", "Cell attachment RGD", "R-G-D", True),
    ("PS00014", "ER retention signal", "[KRHQSA]-[DENQ]-E-L>", False),
    ("PS00342", "Peroxisome targeting PTS1", "[STAGCN]-[RKH]-[LIVMAFY]>", False),
    ("ROSSMANN", "Rossmann NAD(P)-binding GxGxxG", "[IVLA]-x-[IVLA]-x-G-x-G-x(2)-G", False),
    ("WALKER_B", "Walker B (DEAD/DEAH-like)", "[LIVMF](4)-D-E", False),
    ("SER_HYDROLASE", "Ser hydrolase nucleophile elbow GxSxG", "G-x-S-x-G", False),
    ("CUPIN_HXH", "Cupin-like metal-binding HxH...E", "H-x-H-x(3,4)-E", False),
    ("PS00001", "N-glycosylation site", "N-{P}-[ST]-{P}", True),
    ("PS00004", "cAMP/cGMP kinase phospho site", "[RK](2)-x-[ST]", True),
    ("PS00005", "PKC phosphorylation site", "[ST]-x-[RK]", True),
    ("PS00006", "CK2 phosphorylation site", "[ST]-x(2)-[DE]", True),
    ("PS00008", "N-myristoylation site", "G-{EDRKHPFYW}-x(2)-[STAGCN]-{P}", True),
    ("PS00009", "Amidation site", "x-G-[RK]-[RK]", True),
]


def prosite_to_regex(pattern: str) -> str:
    """Convert PROSITE syntax (e.g. 'C-x(2,4)-[DE]-{P}>') to a Python regex."""
    p = pattern.strip().rstrip(".")
    start = p.startswith("<")
    end = p.endswith(">")
    p = p.strip("<>")
    out = []
    for el in p.split("-"):
        m = re.fullmatch(r"(.+?)(?:\((\d+)(?:,(\d+))?\))?", el)
        core, a, b = m.group(1), m.group(2), m.group(3)
        core = core.replace(">", "")
        if core == "x":
            r = "."
        elif core.startswith("["):
            r = core
        elif core.startswith("{"):
            r = "[^" + core[1:-1] + "]"
        else:
            r = core
        if a:
            r += "{" + a + ("," + b if b else "") + "}"
        out.append(r)
    rx = "".join(out)
    return ("^" if start else "") + rx + ("$" if end else "")


def load_prosite_dat(path: str | Path) -> list[tuple[str, str, str, bool]]:
    """Parse PATTERN entries from prosite.dat."""
    entries, cur = [], {}
    for line in Path(path).read_text(errors="ignore").splitlines():
        tag, val = line[:2], line[5:]
        if tag == "ID":
            cur = {"type": val.split(";")[1].strip().rstrip(".")}
        elif tag == "AC":
            cur["ac"] = val.strip().rstrip(";")
        elif tag == "DE":
            cur["de"] = val.strip()
        elif tag == "PA":
            cur["pa"] = cur.get("pa", "") + val.strip()
        elif tag == "CC" and "/SKIP-FLAG=TRUE" in val:
            cur["skip"] = True
        elif tag == "//":
            if cur.get("type") == "PATTERN" and cur.get("pa"):
                entries.append((cur["ac"], cur.get("de", ""), cur["pa"], cur.get("skip", False)))
            cur = {}
    return entries


def scan(seq: str, patterns=None, user_regex: dict[str, str] | None = None,
         include_frequent: bool = False) -> list[dict]:
    """All motif matches in one protein (overlapping matches included)."""
    patterns = BUILTIN if patterns is None else patterns
    hits = []
    for ac, name, pat, frequent in patterns:
        if frequent and not include_frequent:
            continue
        try:
            rx = re.compile("(?=(" + prosite_to_regex(pat) + "))")
        except re.error:
            continue
        for m in rx.finditer(seq):
            hits.append({"id": ac, "motif": name, "start": m.start() + 1,
                         "end": m.start() + len(m.group(1)), "match": m.group(1)})
    for name, rx in (user_regex or {}).items():
        for m in re.finditer("(?=(" + rx + "))", seq):
            hits.append({"id": "user", "motif": name, "start": m.start() + 1,
                         "end": m.start() + len(m.group(1)), "match": m.group(1)})
    return hits


def summarise(hits: list[dict]) -> str:
    c = Counter(h["motif"] for h in hits)
    return "; ".join(f"{k} x{v}" if v > 1 else k for k, v in c.items())


# --------------------------------------------------------------- de novo
def denovo_kmers(seqs: dict[str, str], k_range=(4, 5, 6), min_frac: float = 0.5,
                 shuffles: int = 3, top: int = 15, seed: int = 0) -> list[dict]:
    """Find k-mers present in >= min_frac of sequences, scored by enrichment
    over shuffled copies (preserves composition). Allows one wildcard position."""
    rnd = random.Random(seed)
    n = len(seqs)
    if n < 3:
        return []

    def occurrence(sset):
        occ = Counter()
        for s in sset:
            seen = set()
            for k in k_range:
                for i in range(len(s) - k + 1):
                    w = s[i:i + k]
                    seen.add(w)
                    for j in range(1, k - 1):
                        seen.add(w[:j] + "x" + w[j + 1:])
            occ.update(seen)
        return occ

    fg = occurrence(seqs.values())
    bg = Counter()
    for _ in range(shuffles):
        sh = []
        for s in seqs.values():
            l = list(s)
            rnd.shuffle(l)
            sh.append("".join(l))
        bg.update(occurrence(sh))
    rows = []
    for w, c in fg.items():
        if c / n < min_frac or w.count("x") > 1:
            continue
        exp = (bg.get(w, 0) + 1) / (shuffles + 1)
        rows.append({"motif": w, "seq_frac": round(c / n, 2), "enrichment": round(c / exp, 1),
                     "length": len(w)})
    rows.sort(key=lambda r: (r["enrichment"] * r["seq_frac"] * (r["length"] - r["motif"].count("x"))), reverse=True)
    # drop motifs contained in a better one
    keep = []
    for r in rows:
        core = r["motif"].replace("x", ".")
        if not any(re.search(core, k["motif"].replace("x", "")) for k in keep):
            keep.append(r)
        if len(keep) >= top:
            break
    return keep


def denovo_streme(seqs: dict[str, str], outdir: Path, nmotifs: int = 5) -> Path:
    if not shutil.which("streme"):
        raise RuntimeError("STREME (MEME suite) not installed: conda install -c bioconda meme")
    with tempfile.NamedTemporaryFile("w", suffix=".faa", delete=False) as fh:
        for k, v in seqs.items():
            fh.write(f">{k}\n{v}\n")
    subprocess.run(["streme", "--p", fh.name, "--protein", "--nmotifs", str(nmotifs),
                    "--oc", str(outdir)], check=True, capture_output=True)
    return outdir / "streme.txt"
