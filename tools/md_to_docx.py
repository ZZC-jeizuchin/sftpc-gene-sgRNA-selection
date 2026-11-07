#!/usr/bin/env python3
# 本项目由华南师范大学附属中学知识城校区高二一班同学制作
# Copyright (C) 2026 华南师范大学附属中学知识城校区高二一班
# SPDX-License-Identifier: AGPL-3.0-or-later
"""把项目里的 Markdown 文档转成排版规范的 Word（.docx）。

只处理本项目用到的那几种 Markdown 语法：标题、段落、表格、引用、代码块、
无序/有序列表、粗体、行内代码、水平线。
中文字体通过 w:eastAsia 显式指定（python-docx 只设 font.name 对汉字无效）。
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Pt, RGBColor

BODY_EA = "宋体"
HEAD_EA = "黑体"
LATIN = "Times New Roman"
MONO = "Consolas"

INK = RGBColor(0x1A, 0x1A, 0x1A)
GREY = RGBColor(0x55, 0x5F, 0x6D)
ACCENT = RGBColor(0x1F, 0x4E, 0x79)


# --------------------------------------------------------------------------
# 底层工具
# --------------------------------------------------------------------------

def set_fonts(run, latin=LATIN, ea=BODY_EA):
    """同时设置西文与中文字体。"""
    run.font.name = latin
    rpr = run._element.get_or_add_rPr()
    rfonts = rpr.find(qn("w:rFonts"))
    if rfonts is None:
        rfonts = OxmlElement("w:rFonts")
        rpr.append(rfonts)
    rfonts.set(qn("w:ascii"), latin)
    rfonts.set(qn("w:hAnsi"), latin)
    rfonts.set(qn("w:eastAsia"), ea)


def shade(element, hex_fill):
    shd = OxmlElement("w:shd")
    shd.set(qn("w:val"), "clear")
    shd.set(qn("w:color"), "auto")
    shd.set(qn("w:fill"), hex_fill)
    element.append(shd)


def para_bottom_border(paragraph, color="D0D7DE", size=6):
    ppr = paragraph._p.get_or_add_pPr()
    borders = OxmlElement("w:pBdr")
    bottom = OxmlElement("w:bottom")
    bottom.set(qn("w:val"), "single")
    bottom.set(qn("w:sz"), str(size))
    bottom.set(qn("w:space"), "1")
    bottom.set(qn("w:color"), color)
    borders.append(bottom)
    ppr.append(borders)


def add_page_number(section):
    footer = section.footer
    p = footer.paragraphs[0] if footer.paragraphs else footer.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = p.add_run()
    fld_begin = OxmlElement("w:fldChar")
    fld_begin.set(qn("w:fldCharType"), "begin")
    instr = OxmlElement("w:instrText")
    instr.set(qn("xml:space"), "preserve")
    instr.text = " PAGE "
    fld_end = OxmlElement("w:fldChar")
    fld_end.set(qn("w:fldCharType"), "end")
    run._r.append(fld_begin)
    run._r.append(instr)
    run._r.append(fld_end)
    set_fonts(run)
    run.font.size = Pt(9)
    run.font.color.rgb = GREY


# --------------------------------------------------------------------------
# 行内富文本
# --------------------------------------------------------------------------

INLINE_RE = re.compile(
    r"(\*\*.+?\*\*|`[^`]+`|\[[^\]]+\]\([^)]+\))"
)


def add_rich(paragraph, text, size=10.5, ea=BODY_EA, color=INK, bold=False):
    """把含 **粗体**、`代码`、[链接](url) 的文本写进段落。"""
    text = text.replace("\\|", "|")
    for part in INLINE_RE.split(text):
        if not part:
            continue
        if part.startswith("**") and part.endswith("**") and len(part) > 4:
            run = paragraph.add_run(part[2:-2])
            run.bold = True
            set_fonts(run, ea=ea)
            run.font.size = Pt(size)
            run.font.color.rgb = color
        elif part.startswith("`") and part.endswith("`") and len(part) > 2:
            run = paragraph.add_run(part[1:-1])
            set_fonts(run, latin=MONO, ea=ea)
            run.font.size = Pt(size - 0.5)
            run.font.color.rgb = RGBColor(0xB0, 0x30, 0x60)
        elif part.startswith("[") and "](" in part:
            m = re.match(r"\[([^\]]+)\]\(([^)]+)\)", part)
            label, url = (m.group(1), m.group(2)) if m else (part, "")
            run = paragraph.add_run(label)
            run.font.underline = True
            run.font.color.rgb = RGBColor(0x25, 0x63, 0xA8)
            set_fonts(run, ea=ea)
            run.font.size = Pt(size)
            if url and url.startswith("http"):
                run2 = paragraph.add_run(f"（{url}）")
                set_fonts(run2, ea=ea)
                run2.font.size = Pt(size - 1.5)
                run2.font.color.rgb = GREY
        else:
            run = paragraph.add_run(part)
            run.bold = bold
            set_fonts(run, ea=ea)
            run.font.size = Pt(size)
            run.font.color.rgb = color


# --------------------------------------------------------------------------
# 表格
# --------------------------------------------------------------------------

def split_row(line):
    line = line.strip()
    if line.startswith("|"):
        line = line[1:]
    if line.endswith("|"):
        line = line[:-1]
    # 处理转义竖线
    parts, buf, i = [], "", 0
    while i < len(line):
        if line[i] == "\\" and i + 1 < len(line) and line[i + 1] == "|":
            buf += "|"
            i += 2
            continue
        if line[i] == "|":
            parts.append(buf)
            buf = ""
            i += 1
            continue
        buf += line[i]
        i += 1
    parts.append(buf)
    return [p.strip() for p in parts]


def add_table(doc, rows):
    header, body = rows[0], rows[1:]
    ncol = len(header)
    table = doc.add_table(rows=1, cols=ncol)
    table.style = "Table Grid"
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.autofit = True

    hdr = table.rows[0].cells
    for i, text in enumerate(header):
        hdr[i].text = ""
        p = hdr[i].paragraphs[0]
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        add_rich(p, text, size=9, bold=True, color=RGBColor(0xFF, 0xFF, 0xFF))
        shade(hdr[i]._tc.get_or_add_tcPr(), "1F4E79")

    for row in body:
        cells = table.add_row().cells
        for i in range(ncol):
            cells[i].text = ""
            p = cells[i].paragraphs[0]
            val = row[i] if i < len(row) else ""
            add_rich(p, val, size=9)

    doc.add_paragraph()


# --------------------------------------------------------------------------
# Markdown -> docx
# --------------------------------------------------------------------------

CJK_RANGES = (
    (0x3000, 0x303F), (0x4E00, 0x9FFF), (0x3400, 0x4DBF),
    (0xFF00, 0xFFEF), (0x2000, 0x206F),
)


def _is_cjk(ch: str) -> bool:
    code = ord(ch)
    return any(lo <= code <= hi for lo, hi in CJK_RANGES)


def join_soft_lines(parts):
    """把 Markdown 里的软换行接回一段。

    中文之间直接相连（不插空格），中英边界插一个空格。
    """
    out = parts[0] if parts else ""
    for seg in parts[1:]:
        if not seg:
            continue
        if not out:
            out = seg
            continue
        if _is_cjk(out[-1]) and _is_cjk(seg[0]):
            out += seg
        else:
            out += " " + seg
    return out


BLOCK_START_RES = [
    re.compile(r"^```"), re.compile(r"^\|"), re.compile(r"^#{1,6}\s"),
    re.compile(r"^(-{3,}|\*{3,})$"), re.compile(r"^>"),
    re.compile(r"^[-*+]\s+"), re.compile(r"^\d+[.、)]\s+"),
]


def is_block_start(line: str) -> bool:
    st = line.strip()
    if not st:
        return True
    return any(r.match(st) for r in BLOCK_START_RES)


def convert(md_path: Path, docx_path: Path):
    lines = md_path.read_text(encoding="utf-8").splitlines()
    doc = Document()

    sec = doc.sections[0]
    sec.page_width, sec.page_height = Cm(21.0), Cm(29.7)
    sec.top_margin = Cm(2.4)
    sec.bottom_margin = Cm(2.2)
    sec.left_margin = Cm(2.6)
    sec.right_margin = Cm(2.4)
    add_page_number(sec)

    normal = doc.styles["Normal"]
    normal.font.name = LATIN
    normal.font.size = Pt(10.5)
    normal.element.rPr.rFonts.set(qn("w:eastAsia"), BODY_EA)
    normal.paragraph_format.line_spacing_rule = WD_LINE_SPACING.MULTIPLE
    normal.paragraph_format.line_spacing = 1.35
    normal.paragraph_format.space_after = Pt(5)

    for lvl, size in ((1, 15), (2, 12.5), (3, 11)):
        st = doc.styles[f"Heading {lvl}"]
        st.font.name = LATIN
        st.font.size = Pt(size)
        st.font.bold = True
        st.font.color.rgb = ACCENT
        st.element.rPr.rFonts.set(qn("w:eastAsia"), HEAD_EA)
        st.paragraph_format.space_before = Pt(14 if lvl == 1 else 10)
        st.paragraph_format.space_after = Pt(6)

    i, n = 0, len(lines)
    first_h1 = None
    while i < n:
        raw = lines[i]
        line = raw.rstrip()
        stripped = line.strip()

        if not stripped:
            i += 1
            continue

        # 代码块
        if stripped.startswith("```"):
            i += 1
            buf = []
            while i < n and not lines[i].strip().startswith("```"):
                buf.append(lines[i])
                i += 1
            i += 1
            p = doc.add_paragraph()
            p.paragraph_format.left_indent = Cm(0.6)
            p.paragraph_format.space_before = Pt(4)
            p.paragraph_format.space_after = Pt(8)
            p.paragraph_format.line_spacing = 1.05
            shade(p._p.get_or_add_pPr(), "F4F6F8")
            for k, code_line in enumerate(buf):
                run = p.add_run(code_line if k == 0 else "\n" + code_line)
                set_fonts(run, latin=MONO, ea=BODY_EA)
                run.font.size = Pt(8.5)
            continue

        # 表格
        if stripped.startswith("|") and i + 1 < n and re.match(
                r"^\s*\|[\s:|-]+\|\s*$", lines[i + 1]):
            rows = [split_row(stripped)]
            i += 2
            while i < n and lines[i].strip().startswith("|"):
                rows.append(split_row(lines[i].strip()))
                i += 1
            add_table(doc, rows)
            continue

        # 水平线
        if re.fullmatch(r"-{3,}|\*{3,}", stripped):
            p = doc.add_paragraph()
            p.paragraph_format.space_before = Pt(2)
            p.paragraph_format.space_after = Pt(8)
            para_bottom_border(p)
            i += 1
            continue

        # 标题
        m = re.match(r"^(#{1,6})\s+(.*)$", stripped)
        if m:
            level, text = len(m.group(1)), m.group(2).strip()
            if level == 1 and first_h1 is None:
                first_h1 = text
                p = doc.add_paragraph()
                p.alignment = WD_ALIGN_PARAGRAPH.CENTER
                p.paragraph_format.space_before = Pt(6)
                p.paragraph_format.space_after = Pt(10)
                add_rich(p, text, size=19, ea=HEAD_EA, color=ACCENT, bold=True)
                i += 1
                continue
            # 副标题：紧跟主标题、以 —— 开头
            if level == 2 and text.startswith("——") and i < 4:
                p = doc.add_paragraph()
                p.alignment = WD_ALIGN_PARAGRAPH.CENTER
                p.paragraph_format.space_after = Pt(16)
                add_rich(p, text, size=13, ea=HEAD_EA, color=GREY)
                i += 1
                continue
            lvl = min(level - 1, 3) if first_h1 else min(level, 3)
            p = doc.add_heading(level=lvl)
            add_rich(p, text, size={1: 15, 2: 12.5, 3: 11}[lvl], ea=HEAD_EA,
                     color=ACCENT, bold=True)
            i += 1
            continue

        # 引用
        if stripped.startswith(">"):
            buf = []
            while i < n and lines[i].strip().startswith(">"):
                buf.append(lines[i].strip().lstrip(">").strip())
                i += 1
            p = doc.add_paragraph()
            p.paragraph_format.left_indent = Cm(0.7)
            p.paragraph_format.space_before = Pt(4)
            p.paragraph_format.space_after = Pt(8)
            ppr = p._p.get_or_add_pPr()
            borders = OxmlElement("w:pBdr")
            left = OxmlElement("w:left")
            left.set(qn("w:val"), "single")
            left.set(qn("w:sz"), "18")
            left.set(qn("w:space"), "8")
            left.set(qn("w:color"), "1F4E79")
            borders.append(left)
            ppr.append(borders)
            shade(ppr, "F7F9FC")
            add_rich(p, join_soft_lines(buf), size=10, color=GREY)
            continue

        # 无序列表（手动圆点 + 悬挂缩进，不依赖模板里的编号定义）
        m = re.match(r"^[-*+]\s+(.*)$", stripped)
        if m:
            buf = [m.group(1)]
            i += 1
            while i < n and not is_block_start(lines[i]):
                buf.append(lines[i].strip())
                i += 1
            p = doc.add_paragraph()
            p.paragraph_format.left_indent = Cm(0.9)
            p.paragraph_format.first_line_indent = Cm(-0.45)
            p.paragraph_format.space_after = Pt(3)
            add_rich(p, "• " + join_soft_lines(buf))
            continue

        # 有序列表（保留原文编号，避免跨列表自动续号）
        m = re.match(r"^(\d+)[.、)]\s+(.*)$", stripped)
        if m:
            buf = [m.group(2)]
            i += 1
            while i < n and not is_block_start(lines[i]):
                buf.append(lines[i].strip())
                i += 1
            p = doc.add_paragraph()
            p.paragraph_format.left_indent = Cm(0.95)
            p.paragraph_format.first_line_indent = Cm(-0.6)
            p.paragraph_format.space_after = Pt(3)
            add_rich(p, m.group(1) + ". " + join_soft_lines(buf))
            continue

        # 普通段落：把软换行接回去
        buf = [stripped]
        i += 1
        while i < n and not is_block_start(lines[i]):
            buf.append(lines[i].strip())
            i += 1
        p = doc.add_paragraph()
        p.paragraph_format.first_line_indent = Cm(0)
        add_rich(p, join_soft_lines(buf))

    doc.save(docx_path)
    return docx_path


def main(argv):
    base = Path(argv[1]) if len(argv) > 1 else Path(__file__).resolve().parent.parent
    # (Markdown 源, 目标 docx)。源和产物可以在不同目录。
    pairs = [
        ("docs_src/参赛报告_终稿.md",        "参赛报告_终稿.docx"),
        ("docs_src/课题全解_思路与方法.md",   "课题全解_思路与方法.docx"),
        ("申报材料/01_查新报告_模板.md",      "申报材料/01_查新报告.docx"),
        ("申报材料/02_研究日志_模板.md",      "申报材料/02_研究日志.docx"),
    ]
    for src, dst in pairs:
        sp, dp = base / src, base / dst
        if not sp.is_file():
            print(f"[skip] {src}")
            continue
        convert(sp, dp)
        print(f"[ok] {dst}  ({dp.stat().st_size/1024:.1f} KB)")


if __name__ == "__main__":
    main(sys.argv)
