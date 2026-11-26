"""Primary, secondary and small-scale structural features from sequence (and
optionally from a 3D model).

All predictors here are fast, dependency-free *first-pass* methods meant for
triaging hundreds of hits. They are clearly labelled in the output. For
publication-grade calls use the dedicated tools named in each docstring.
"""
from __future__ import annotations

import math
import re
from collections import Counter
from pathlib import Path

import numpy as np

KD = dict(I=4.5, V=4.2, L=3.8, F=2.8, C=2.5, M=1.9, A=1.8, G=-0.4, T=-0.7, S=-0.8, W=-0.9,
          Y=-1.3, P=-1.6, H=-3.2, E=-3.5, Q=-3.5, D=-3.5, N=-3.5, K=-3.9, R=-4.5)
# Chou-Fasman propensities (P-alpha, P-beta, P-turn) x100
CF = dict(A=(142, 83, 66), R=(98, 93, 95), N=(67, 89, 156), D=(101, 54, 146), C=(70, 119, 119),
          Q=(111, 110, 98), E=(151, 37, 74), G=(57, 75, 156), H=(100, 87, 95), I=(108, 160, 47),
          L=(121, 130, 59), K=(114, 74, 101), M=(145, 105, 60), F=(113, 138, 60), P=(57, 55, 152),
          S=(77, 75, 143), T=(83, 119, 96), W=(108, 137, 96), Y=(69, 147, 114), V=(106, 170, 50))


def _win(values: list[float], w: int) -> np.ndarray:
    v = np.asarray(values, dtype=float)
    if len(v) < w:
        return np.full(len(v), v.mean() if len(v) else 0.0)
    k = np.convolve(v, np.ones(w) / w, mode="same")
    return k


def _segments(mask, min_len: int = 1, offset: int = 1) -> list[tuple[int, int]]:
    segs, start = [], None
    for i, m in enumerate(list(mask) + [False]):
        if m and start is None:
            start = i
        elif not m and start is not None:
            if i - start >= min_len:
                segs.append((start + offset, i - 1 + offset))
            start = None
    return segs


# ------------------------------------------------------------ primary structure
def composition(seq: str) -> dict:
    """Amino-acid composition (%) and grouped classes."""
    n = max(len(seq), 1)
    c = Counter(seq)
    groups = {"hydrophobic": "AVILMFWC", "aromatic": "FWYH", "charged": "DEKR",
              "acidic": "DE", "basic": "KR", "polar": "STNQ", "tiny": "AGS"}
    out = {f"pct_{aa}": round(100 * c.get(aa, 0) / n, 1) for aa in "ACDEFGHIKLMNPQRSTVWY"}
    out.update({f"pct_{g}": round(100 * sum(c.get(a, 0) for a in s) / n, 1) for g, s in groups.items()})
    out["n_cys"] = c.get("C", 0)
    return out


def solubility_wh(seq: str) -> dict:
    """Revised Wilkinson-Harrison solubility model (E. coli overexpression)."""
    n = max(len(seq), 1)
    c = Counter(seq)
    cv = 15.43 * (c["N"] + c["G"] + c["P"] + c["S"]) / n \
        - 29.56 * abs((c["R"] + c["K"] - c["D"] - c["E"]) / n - 0.03)
    d = cv - 1.71
    p = 0.4934 + 0.276 * abs(d) - 0.0392 * d * d
    return {"wh_soluble": d <= 0, "wh_probability": round(min(max(p, 0), 1), 2)}


def low_complexity(seq: str, w: int = 12, max_entropy: float = 2.2) -> list[tuple[int, int]]:
    """SEG-like low-complexity regions: windows with Shannon entropy < max_entropy bits."""
    if len(seq) < w:
        return []
    mask = [False] * len(seq)
    for i in range(len(seq) - w + 1):
        cnt = Counter(seq[i:i + w])
        h = -sum(v / w * math.log2(v / w) for v in cnt.values())
        if h < max_entropy:
            for j in range(i, i + w):
                mask[j] = True
    return _segments(mask, min_len=w)


# ---------------------------------------------------------- secondary structure
def secondary_structure(seq: str) -> str:
    """Chou-Fasman style H/E/C string (~55-60 % Q3). Quick estimate only:
    use PSIPRED, S4PRED or the structure-derived string (`ss_from_structure`)."""
    pa = _win([CF.get(a, (100, 100, 100))[0] for a in seq], 6)
    pb = _win([CF.get(a, (100, 100, 100))[1] for a in seq], 5)
    pt = _win([CF.get(a, (100, 100, 100))[2] for a in seq], 4)
    ss = []
    for a, b, t in zip(pa, pb, pt):
        if t > 115 and t > a and t > b:
            ss.append("C")
        elif a > 103 and a >= b:
            ss.append("H")
        elif b > 105 and b > a:
            ss.append("E")
        else:
            ss.append("C")
    s = "".join(ss)
    s = re.sub(r"(?<!H)H{1,3}(?!H)", lambda m: "C" * len(m.group()), s)  # helices >= 4
    s = re.sub(r"(?<!E)E{1,2}(?!E)", lambda m: "C" * len(m.group()), s)  # strands >= 3
    return s


def ss_from_structure(pdb: str | Path) -> str:
    """H/E/C from CA geometry of a model (P-SEA-like distance criteria)."""
    ca = []
    for line in Path(pdb).read_text().splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            ca.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            if len(ca) > 1 and line[21] != prev_chain:
                break
            prev_chain = line[21]
    ca = np.array(ca)
    n = len(ca)
    ss = ["C"] * n
    d = lambda i, j: float(np.linalg.norm(ca[i] - ca[j]))  # noqa: E731
    for i in range(n - 4):
        if 4.6 <= d(i, i + 3) <= 5.8 and 5.6 <= d(i, i + 4) <= 6.8 and 4.8 <= d(i, i + 2) <= 5.9:
            for j in range(i, i + 5):
                ss[j] = "H"
    for i in range(n - 2):
        if ss[i] != "H" and d(i, i + 2) >= 6.2:
            ss[i + 1] = "E" if ss[i + 1] != "H" else "H"
    s = "".join(ss)
    s = re.sub(r"(?<!E)E{1,2}(?!E)", lambda m: "C" * len(m.group()), s)
    return s


def ss_summary(ss: str) -> dict:
    n = max(len(ss), 1)
    return {"ss_helix": round(ss.count("H") / n, 2), "ss_strand": round(ss.count("E") / n, 2),
            "ss_coil": round(ss.count("C") / n, 2), "n_helices": len(re.findall("H+", ss)),
            "n_strands": len(re.findall("E+", ss))}


def supersecondary(ss: str) -> list[dict]:
    """Super-secondary motifs from an H/E/C string: helix-turn-helix, alpha-hairpin,
    beta-hairpin and beta-alpha-beta (Rossmann-like) units."""
    pats = {
        "helix-turn-helix": r"H{7,20}C{2,5}H{7,}",
        "alpha-hairpin": r"H{10,}C{2,8}H{10,}",
        "beta-hairpin": r"E{3,}C{2,5}E{3,}",
        "beta-alpha-beta": r"E{3,}C{1,10}H{8,}C{1,10}E{3,}",
    }
    out = []
    for name, p in pats.items():
        for m in re.finditer(f"(?=({p}))", ss):
            s, e = m.start(1) + 1, m.start(1) + len(m.group(1))
            if not any(o["motif"] == name and o["start"] <= s <= o["end"] for o in out):
                out.append({"motif": name, "start": s, "end": e})
    return out


# ------------------------------------------------------- membrane and targeting
def transmembrane(seq: str, w: int = 19, cutoff: float = 1.6) -> list[tuple[int, int]]:
    """Kyte-Doolittle hydropathy TM helices (window 19, mean >= 1.6).
    Use DeepTMHMM / TMbed for reliable topology."""
    h = _win([KD.get(a, 0) for a in seq], w)
    segs = _segments(h >= cutoff, min_len=12)
    return [(max(1, s - 3), min(len(seq), e + 3)) for s, e in segs]


def signal_peptide(seq: str) -> dict:
    """Heuristic Sec signal peptide (n/h/c regions + AxA cleavage), Tat motif and
    lipobox. Use SignalP 6 for real predictions."""
    n_term = seq[:45]
    res = {"signal_peptide": False, "sp_cleavage": None, "tat_motif": False, "lipobox": False}
    if len(n_term) < 20:
        return res
    pos_n = sum(a in "KR" for a in n_term[1:8]) >= 1
    acidic_ok = sum(a in "DE" for a in n_term[:20]) <= 1   # signal peptides are rarely acidic
    h = _win([KD.get(a, 0) for a in n_term], 8)
    core = [i for i, v in enumerate(h) if v >= 1.6 and 5 <= i <= 25]
    if len(core) >= 3 and pos_n and acidic_ok:
        end_core = core[-1] + 4
        rng = range(max(end_core, 15), min(36, len(n_term) - 1, end_core + 12))
        strict = [c for c in rng if n_term[c - 3] in "AV" and n_term[c - 1] == "A" and n_term[c] != "P"]
        loose = [c for c in rng if n_term[c - 3] in "AGSCTV" and n_term[c - 1] in "AGSC" and n_term[c] != "P"]
        if strict or loose:
            res.update(signal_peptide=True, sp_cleavage=(strict or loose)[0])
    res["tat_motif"] = bool(re.search(r"[ST]RR.[FL][LIMVA]K?", seq[:40]))
    lb = re.search(r"[LVI][ASTVI][GAS]C", seq[10:40])
    res["lipobox"] = bool(lb)
    return res


TOP_IDP = dict(W=-0.884, F=-0.697, Y=-0.510, I=-0.486, M=-0.397, L=-0.326, V=-0.121, N=0.007,
               C=0.02, T=0.059, A=0.06, G=0.166, R=0.180, D=0.192, H=0.303, Q=0.318, K=0.586,
               S=0.341, E=0.736, P=0.987)


def disorder(seq: str, w: int = 31, threshold: float = 0.2, min_len: int = 30) -> dict:
    """TOP-IDP disorder-propensity scale averaged over a window. Use IUPred3 /
    AIUPred, or pLDDT < 50 from a predicted model, for better calls."""
    v = _win([TOP_IDP.get(a, 0) for a in seq], w)
    mask = v > threshold
    return {"disorder_frac": round(float(mask.mean()) if len(seq) else 0, 2),
            "disordered_regions": _segments(mask, min_len=min_len)}


def coiled_coils(seq: str, w: int = 28) -> list[tuple[int, int]]:
    """Heptad-repeat scan: >= 80 % L/I/V/M at a/d, >= 65 % polar/charged at
    b/c/e/f/g, no proline. ~5 % of bacterial proteins hit. Use DeepCoil /
    CoCoNat or Paircoil2 for real predictions."""
    if len(seq) < w:
        return []
    polar = set("EKQRSDTNAH")
    mask = np.zeros(len(seq), bool)
    for i in range(len(seq) - w + 1):
        win = seq[i:i + w]
        if "P" in win:
            continue
        for reg in range(7):
            ad = [win[j] for j in range(w) if (j + reg) % 7 in (0, 3)]
            other = [win[j] for j in range(w) if (j + reg) % 7 in (1, 2, 4, 5, 6)]
            if sum(a in "LIVM" for a in ad) / len(ad) >= 0.8 and \
                    sum(a in polar for a in other) / len(other) >= 0.65:
                mask[i:i + w] = True
                break
    return _segments(mask, min_len=w)


# ------------------------------------------------------------------ disulfides
def disulfides_from_structure(pdb: str | Path, max_sg: float = 2.5, max_cb: float = 4.5) -> list[tuple[int, int]]:
    """Cys pairs bonded (SG-SG <= 2.5 A) or geometrically compatible (CB-CB <= 4.5 A)."""
    sg, cb = {}, {}
    for line in Path(pdb).read_text().splitlines():
        if line.startswith("ATOM") and line[17:20] == "CYS":
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            rn = int(line[22:26])
            if line[12:16].strip() == "SG":
                sg[rn] = xyz
            elif line[12:16].strip() == "CB":
                cb[rn] = xyz
    pairs = []
    ids = sorted(set(sg) | set(cb))
    for i, a in enumerate(ids):
        for b in ids[i + 1:]:
            if b - a < 3:
                continue
            if a in sg and b in sg and np.linalg.norm(sg[a] - sg[b]) <= max_sg:
                pairs.append((a, b))
            elif a in cb and b in cb and np.linalg.norm(cb[a] - cb[b]) <= max_cb:
                pairs.append((a, b))
    return pairs


# --------------------------------------------------------------- all-in-one
def fmt_regions(regs) -> str:
    return ";".join(f"{a}-{b}" for a, b in regs)


def profile(seq: str) -> dict:
    """Every sequence-only feature in one flat dict (for tables)."""
    ss = secondary_structure(seq)
    tm = transmembrane(seq)
    sp = signal_peptide(seq)
    if sp["signal_peptide"]:  # a signal peptide's h-region is not a TM helix
        tm = [t for t in tm if t[0] > sp["sp_cleavage"]]
    dis = disorder(seq)
    cc = coiled_coils(seq)
    lc = low_complexity(seq)
    row = {
        "ss_pred": ss, **ss_summary(ss),
        "supersecondary": ",".join(sorted({m["motif"] for m in supersecondary(ss)})),
        "n_tm": len(tm), "tm_regions": fmt_regions(tm),
        **sp,
        "disorder_frac": dis["disorder_frac"], "disordered_regions": fmt_regions(dis["disordered_regions"]),
        "coiled_coils": fmt_regions(cc), "low_complexity": fmt_regions(lc),
        "n_cys": seq.count("C"),
        "localisation_guess": ("membrane" if len(tm) >= 2 or (tm and not sp["signal_peptide"])
                               else "lipoprotein" if sp["lipobox"] and sp["signal_peptide"]
                               else "secreted/periplasmic" if sp["signal_peptide"] or sp["tat_motif"]
                               else "cytoplasmic"),
    }
    row.update(solubility_wh(seq))
    return row
