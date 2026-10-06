"""把 CHOPCHOP 的 results.tsv 解析成结构化记录。

表头（实测）：
    Rank  Target sequence  Genomic location  Strand  GC content (%)
    Self-complementarity  MM0  MM1  MM2  MM3  Efficiency

按表头名取列，官方加列不会解析错。
"""

from __future__ import annotations

import re
from dataclasses import dataclass, fields

_LOC_RE = re.compile(r"^\s*([^:]+)\s*:\s*([\d,]+)\s*$")


def _to_num(v):
    if v is None:
        return None
    v = str(v).strip()
    if v in ("", "-", "NA", "nan", "None"):
        return None
    try:
        return int(v)
    except ValueError:
        pass
    try:
        return float(v)
    except ValueError:
        return None


@dataclass
class GuideRecord:
    """一条候选 guide（CHOPCHOP results.tsv 一行）。"""

    rank: int | None
    target_seq: str          # 23 nt = 20 nt spacer + PAM
    chrom: str               # "chr8"
    location: int | None     # 23mer 正链最左端 1-based 坐标
    strand: str              # "+" / "-"
    gc_content: float | None
    self_comp: int | None
    mm0: int | None
    mm1: int | None
    mm2: int | None
    mm3: int | None
    efficiency: float | None
    extras: dict

    @property
    def offtarget_count(self) -> int | None:
        vals = [v for v in (self.mm0, self.mm1, self.mm2, self.mm3) if v is not None]
        return sum(vals) if vals else None

    @property
    def genomic_location(self) -> str:
        return f"{self.chrom}:{self.location}" if self.location else ""

    def as_row(self) -> dict:
        row = {f.name: getattr(self, f.name) for f in fields(self) if f.name != "extras"}
        row.update(self.extras)
        return row


_COLUMN_MAP = {
    "Rank": "rank",
    "Target sequence": "target_seq",
    "Genomic location": "location",
    "Strand": "strand",
    "GC content (%)": "gc_content",
    "Self-complementarity": "self_comp",
    "MM0": "mm0",
    "MM1": "mm1",
    "MM2": "mm2",
    "MM3": "mm3",
    "Efficiency": "efficiency",
}

_NUMERIC = {"rank", "gc_content", "self_comp", "mm0", "mm1", "mm2", "mm3", "efficiency"}


def parse_results_tsv(text: str) -> list[GuideRecord]:
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines:
        raise ValueError("results.tsv 内容为空")
    header = [h.strip() for h in lines[0].lstrip("#").split("\t")]
    if "Genomic location" not in header:
        raise ValueError(f"这不是 CHOPCHOP 的 results.tsv；表头={header[:5]}")

    out: list[GuideRecord] = []
    for line in lines[1:]:
        cells = line.split("\t")
        if len(cells) < len(header):
            cells += [""] * (len(header) - len(cells))

        known: dict = {}
        extras: dict = {}
        for col, raw in zip(header, cells):
            name = _COLUMN_MAP.get(col)
            if name is None:
                extras[col] = raw.strip()
                continue
            known[name] = _to_num(raw) if name in _NUMERIC else raw.strip()

        chrom, loc = "", None
        m = _LOC_RE.match(known.get("location", ""))
        if m:
            chrom = m.group(1).strip()
            loc = int(m.group(2).replace(",", ""))

        out.append(GuideRecord(
            rank=known.get("rank"),
            target_seq=known.get("target_seq", ""),
            chrom=chrom,
            location=loc,
            strand=known.get("strand", "+"),
            gc_content=known.get("gc_content"),
            self_comp=known.get("self_comp"),
            mm0=known.get("mm0"),
            mm1=known.get("mm1"),
            mm2=known.get("mm2"),
            mm3=known.get("mm3"),
            efficiency=known.get("efficiency"),
            extras=extras,
        ))
    return out
