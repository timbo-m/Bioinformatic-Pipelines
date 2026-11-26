"""Export a finished run as Markdown, Excel and/or PDF.

    md    report.md (always written)
    xlsx  report.xlsx: Summary sheet, Results (one row per protein), and one
          sheet per stage table (domains, motifs, variants, conservation, ...)
    pdf   report.pdf: the Markdown report typeset with headings, tables and code
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

FORMATS = ("md", "xlsx", "pdf")


def parse_formats(value) -> list[str]:
    if not value:
        return ["md"]
    items = value if isinstance(value, (list, tuple)) else re.split(r"[,\s]+", str(value))
    items = [i.strip().lower().lstrip(".") for i in items if i.strip()]
    items = ["xlsx" if i in ("excel", "xls") else i for i in items]
    if "all" in items:
        return list(FORMATS)
    bad = [i for i in items if i not in FORMATS]
    if bad:
        raise ValueError(f"unknown format(s) {bad}; choose from {', '.join(FORMATS)} or all")
    return sorted(set(items) | {"md"}, key=FORMATS.index)


def export(out: Path, formats, log=print) -> list[Path]:
    out = Path(out)
    written = []
    for f in parse_formats(formats):
        if f == "xlsx":
            written.append(to_excel(out))
        elif f == "pdf":
            written.append(to_pdf(out))
    for w in written:
        log(f"      wrote {w}")
    return written


# --------------------------------------------------------------------- Excel
def _summary_rows(md: str) -> list[dict]:
    rows, section = [], "Run"
    for line in md.splitlines():
        if line.startswith("## "):
            section = line[3:].strip()
        elif line.startswith("- ") and section in ("Run",):
            k, _, v = line[2:].partition(":")
            rows.append({"item": k.strip(), "value": v.strip().strip("`")})
        elif line.startswith("> "):
            rows.append({"item": "note", "value": line[2:].replace("**", "")})
        elif section not in ("Run", "Files", "Final candidates") and line and not line.startswith(("|", "```", "#")):
            if not rows or rows[-1].get("item") != section:
                rows.append({"item": section, "value": _plain(line)})
    return rows


def _plain(s: str) -> str:
    return re.sub(r"[*`]", "", s)


def to_excel(out: Path) -> Path:
    path = out / "report.xlsx"
    md = (out / "report.md").read_text() if (out / "report.md").exists() else ""
    sheets: dict[str, pd.DataFrame] = {"Summary": pd.DataFrame(_summary_rows(md))}
    if (out / "results.tsv").exists():
        sheets["Results"] = pd.read_csv(out / "results.tsv", sep="\t")
    for tsv in sorted(out.glob("*.tsv")):
        if tsv.name == "results.tsv" or tsv.stat().st_size == 0:
            continue
        try:
            sheets[tsv.stem[:31]] = pd.read_csv(tsv, sep="\t")
        except pd.errors.EmptyDataError:
            continue
    with pd.ExcelWriter(path, engine="openpyxl") as xw:
        for name, df in sheets.items():
            df.to_excel(xw, sheet_name=name, index=False)
            ws = xw.sheets[name]
            ws.freeze_panes = "A2"
            if len(df.columns):
                ws.auto_filter.ref = ws.dimensions
            for i, col in enumerate(df.columns, 1):
                width = max([len(str(col))] + [len(str(v)) for v in df[col].head(200)])
                ws.column_dimensions[ws.cell(1, i).column_letter].width = min(max(width + 2, 8), 60)
            for c in ws[1]:
                c.font = c.font.copy(bold=True)
    return path


# ----------------------------------------------------------------------- PDF
def _inline(s: str) -> str:
    s = s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    s = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", s)
    s = re.sub(r"`(.+?)`", r"<font face='Courier'>\1</font>", s)
    s = re.sub(r"(?<![\w*])_(.+?)_(?![\w*])", r"<i>\1</i>", s)
    return s


def to_pdf(out: Path) -> Path:
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4, landscape
    from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
    from reportlab.lib.units import mm
    from reportlab.platypus import (Paragraph, Preformatted, SimpleDocTemplate, Spacer, Table,
                                    TableStyle)

    path = out / "report.pdf"
    md = (out / "report.md").read_text()
    ss = getSampleStyleSheet()
    body = ParagraphStyle("b", parent=ss["BodyText"], fontSize=9, leading=12)
    cell = ParagraphStyle("c", parent=body, fontSize=6.5, leading=8)
    head = ParagraphStyle("h", parent=cell, textColor=colors.white, fontName="Helvetica-Bold")
    code = ParagraphStyle("code", parent=body, fontName="Courier", fontSize=6, leading=7)
    pagesize = landscape(A4)
    width = pagesize[0] - 30 * mm
    story, table, codeblock, in_code = [], [], [], False

    def flush_table():
        if not table:
            return
        rows = [r for r in table if not re.fullmatch(r"\|?(\s*-+\s*\|)+\s*-*\s*", r.replace(":", "-"))]
        data = [[c.strip() for c in r.strip().strip("|").split("|")] for r in rows]
        ncol = max(len(r) for r in data)
        data = [r + [""] * (ncol - len(r)) for r in data]
        cells = [[Paragraph(_inline(v), head if i == 0 else cell) for v in r] for i, r in enumerate(data)]
        t = Table(cells, colWidths=[width / ncol] * ncol, repeatRows=1)
        t.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#2f5d62")),
            ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#eef3f3")]),
            ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#b0c4c4")),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ]))
        story.extend([t, Spacer(1, 6)])
        table.clear()

    for line in md.splitlines():
        if line.startswith("```"):
            if in_code:
                story.append(Preformatted("\n".join(codeblock), code))
                codeblock.clear()
            in_code = not in_code
            continue
        if in_code:
            codeblock.append(line)
            continue
        if line.startswith("|"):
            table.append(line)
            continue
        flush_table()
        if line.startswith("# "):
            story.append(Paragraph(_inline(line[2:]), ss["Title"]))
        elif line.startswith("## "):
            story.append(Paragraph(_inline(line[3:]), ss["Heading2"]))
        elif line.startswith("- "):
            story.append(Paragraph("&bull; " + _inline(line[2:]), body))
        elif line.startswith("> "):
            story.append(Paragraph(_inline(line[2:]), ParagraphStyle("q", parent=body, leftIndent=10,
                                                                      textColor=colors.HexColor("#8a4b00"))))
        elif line.strip():
            story.append(Paragraph(_inline(line), body))
        else:
            story.append(Spacer(1, 4))
    flush_table()
    SimpleDocTemplate(str(path), pagesize=pagesize, leftMargin=15 * mm, rightMargin=15 * mm,
                      topMargin=12 * mm, bottomMargin=12 * mm, title="famscan report").build(story)
    return path
