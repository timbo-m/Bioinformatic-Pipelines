"""Export a famscan output folder as Markdown, Excel and/or PDF.

* md   - ``report.md`` (always written by every run; listed for completeness)
* xlsx - ``report.xlsx``: a Summary sheet (run parameters + report text), the
         full ``results.tsv`` and every other ``*.tsv`` in the folder as its own
         sheet, with frozen headers, autofilter and sensible column widths.
* pdf  - ``report.pdf``: ``report.md`` rendered (headings, lists, tables).

Works on any finished output folder, so you can also convert old runs:
``famscan export famscan_out --format xlsx,pdf``.
"""
from __future__ import annotations

import html
import json
import re
from pathlib import Path

FORMATS = ("md", "xlsx", "pdf")


def parse_formats(spec: str | list | tuple | None) -> list[str]:
    """'md,xlsx' / ['pdf'] / 'all' -> validated list (md always included)."""
    if not spec:
        return ["md"]
    items = spec if isinstance(spec, (list, tuple)) else re.split(r"[,\s]+", str(spec))
    items = [i.strip().lower().lstrip(".") for i in items if i and i.strip()]
    items = ["xlsx" if i in ("excel", "xls") else i for i in items]
    if "all" in items:
        return list(FORMATS)
    bad = [i for i in items if i not in FORMATS]
    if bad:
        raise ValueError(f"unknown format(s) {bad}; choose from {', '.join(FORMATS)} or 'all'")
    return ["md"] + [f for f in FORMATS[1:] if f in items]


def export(out: str | Path, formats=("md",), log=print) -> list[Path]:
    """Write the requested formats into ``out``; returns the files written."""
    out = Path(out)
    if not (out / "report.md").exists():
        raise FileNotFoundError(f"{out}/report.md not found - run a workflow into this folder first")
    written = [out / "report.md"]
    for fmt in parse_formats(formats)[1:]:
        try:
            path = {"xlsx": to_xlsx, "pdf": to_pdf}[fmt](out)
        except ImportError as e:
            log(f"      {fmt} export needs an extra package ({e.name}): pip install openpyxl fpdf2")
            continue
        written.append(path)
        log(f"      wrote {path}")
    return written


# --------------------------------------------------------------------- Excel
def _sheet_name(name: str, used: set) -> str:
    base = re.sub(r"[\[\]:*?/\\]", "_", name)[:31] or "sheet"
    n, i = base, 2
    while n.lower() in used:
        n = f"{base[:28]}_{i}"
        i += 1
    used.add(n.lower())
    return n


def to_xlsx(out: Path) -> Path:
    import pandas as pd
    from openpyxl.styles import Alignment, Font
    from openpyxl.utils import get_column_letter

    dest = out / "report.xlsx"
    tsvs = sorted(out.glob("*.tsv"), key=lambda p: (p.name != "results.tsv", p.name))
    used: set = set()
    with pd.ExcelWriter(dest, engine="openpyxl") as xw:
        # Summary: parameters from run.json, then the report text line by line
        rows = []
        run = out / "run.json"
        if run.exists():
            meta = json.loads(run.read_text())
            if "workflow" in meta:
                rows.append(("workflow", meta["workflow"]))
            for k, v in (meta.get("params") or meta.get("config") or {}).items():
                if v not in (None, "", [], {}):
                    rows.append((k, v if isinstance(v, (int, float)) else str(v)))
            for k, v in (meta.get("timings") or {}).items():
                rows.append((f"time: {k}", v))
        summary = pd.DataFrame(rows, columns=["item", "value"])
        report = pd.DataFrame({"report.md": (out / "report.md").read_text(encoding="utf-8").splitlines()})
        name = _sheet_name("Summary", used)
        summary.to_excel(xw, sheet_name=name, index=False)
        report.to_excel(xw, sheet_name=name, index=False, startrow=len(summary) + 2)
        ws = xw.sheets[name]
        ws.column_dimensions["A"].width = 28
        ws.column_dimensions["B"].width = 80
        for row in ws.iter_rows(min_row=len(summary) + 4):
            for c in row:
                if isinstance(c.value, str) and c.value.startswith("#"):
                    c.font = Font(bold=True, size=13 if c.value.startswith("# ") else 11)

        for tsv in tsvs:
            try:
                df = pd.read_csv(tsv, sep="\t")
            except Exception:  # empty or malformed table - skip it
                continue
            if df.empty and not len(df.columns):
                continue
            name = _sheet_name(tsv.stem, used)
            df.to_excel(xw, sheet_name=name, index=False)
            ws = xw.sheets[name]
            ws.freeze_panes = "B2"
            ws.auto_filter.ref = ws.dimensions
            for c in ws[1]:
                c.font = Font(bold=True)
                c.alignment = Alignment(wrap_text=True, vertical="top")
            for j, col in enumerate(df.columns, 1):
                sample = [str(col)] + [str(v) for v in df[col].head(200).tolist()]
                ws.column_dimensions[get_column_letter(j)].width = min(60, max(8, max(map(len, sample)) + 2))
    return dest


# ----------------------------------------------------------------------- PDF
_FONT_CANDIDATES = [
    ("C:/Windows/Fonts/arial.ttf", "C:/Windows/Fonts/arialbd.ttf",
     "C:/Windows/Fonts/ariali.ttf", "C:/Windows/Fonts/consola.ttf"),
    ("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
     "/usr/share/fonts/truetype/dejavu/DejaVuSans-Oblique.ttf",
     "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf"),
    ("/System/Library/Fonts/Supplemental/Arial.ttf", "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
     "/System/Library/Fonts/Supplemental/Arial Italic.ttf", "/System/Library/Fonts/Supplemental/Courier New.ttf"),
]


def _inline(text: str) -> str:
    """Markdown inline markup -> the small HTML subset fpdf2 understands."""
    t = html.escape(text, quote=False)
    t = re.sub(r"`([^`]+)`", r"<font face='mono'>\1</font>", t)
    t = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", t)
    t = re.sub(r"(?<![\w*])[_*]([^_*]+)[_*](?![\w*])", r"<i>\1</i>", t)
    t = re.sub(r"\[([^\]]+)\]\((https?://[^)]+)\)", r"<a href='\2'>\1</a>", t)
    t = re.sub(r"(?<!['>])(https?://[^\s<)]+)", r"<a href='\1'>\1</a>", t)
    return t


def md_to_html(md: str, max_cols: int = 12) -> str:
    out, para, lst, i = [], [], [], 0
    lines = md.splitlines()

    def flush():
        if para:
            out.append(f"<p>{_inline(' '.join(para))}</p>")
            para.clear()
        if lst:
            out.append("<ul>" + "".join(f"<li>{_inline(x)}</li>" for x in lst) + "</ul>")
            lst.clear()

    while i < len(lines):
        s = lines[i].rstrip()
        if not s.strip():
            flush()
        elif m := re.match(r"(#{1,4})\s+(.*)", s):
            flush()
            n = min(len(m.group(1)), 3)
            out.append(f"<h{n}>{_inline(m.group(2))}</h{n}>")
        elif s.lstrip().startswith("|"):
            flush()
            rows = []
            while i < len(lines) and lines[i].lstrip().startswith("|"):
                cells = [c.strip() for c in lines[i].strip().strip("|").split("|")]
                if not all(re.fullmatch(r":?-{2,}:?", c) for c in cells if c):
                    rows.append(cells)
                i += 1
            if rows:
                head, body = rows[0][:max_cols], [r[:max_cols] for r in rows[1:]]
                note = f" (first {max_cols} of {len(rows[0])} columns; all in report.xlsx/results.tsv)" \
                    if len(rows[0]) > max_cols else ""
                out.append("<table border='1' width='100%'><thead><tr>"
                           + "".join(f"<th>{_inline(h)}</th>" for h in head) + "</tr></thead><tbody>"
                           + "".join("<tr>" + "".join(f"<td>{_inline(c)}</td>" for c in r + [""] * (len(head) - len(r)))
                                     + "</tr>" for r in body) + "</tbody></table>")
                if note:
                    out.append(f"<p><i>{html.escape(note.strip())}</i></p>")
            continue
        elif m := re.match(r"\s*[-*+]\s+(.*)", s):
            if para:
                flush()
            lst.append(m.group(1))
        elif m := re.match(r"\s*>\s?(.*)", s):
            flush()
            out.append(f"<blockquote>{_inline(m.group(1))}</blockquote>")
        else:
            if lst:
                flush()
            para.append(s.strip())
        i += 1
    flush()
    return "\n".join(out)


def to_pdf(out: Path) -> Path:
    import logging

    from fpdf import FPDF

    logging.getLogger("fontTools").setLevel(logging.ERROR)  # noisy font-subsetting chatter
    dest = out / "report.pdf"
    md = (out / "report.md").read_text(encoding="utf-8")
    pdf = FPDF(orientation="landscape", format="A4")
    pdf.set_auto_page_break(True, margin=12)
    pdf.set_margins(12, 12, 12)
    family, mono = "helvetica", "courier"
    for reg, bold, ital, mon in _FONT_CANDIDATES:
        if all(Path(p).exists() for p in (reg, bold, ital, mon)):
            pdf.add_font("body", "", reg)
            pdf.add_font("body", "B", bold)
            pdf.add_font("body", "I", ital)
            pdf.add_font("body", "BI", bold)
            pdf.add_font("mono", "", mon)
            family, mono = "body", "mono"
            break
    else:  # core fonts are latin-1 only
        md = md.encode("latin-1", "replace").decode("latin-1")
    pdf.set_font(family, size=9)
    pdf.add_page()
    body = md_to_html(md)
    if mono != "mono":
        body = body.replace("face='mono'", "face='courier'")
    pdf.write_html(body, font_family=family,
                   tag_styles=None, table_line_separators=True)
    pdf.output(str(dest))
    return dest
