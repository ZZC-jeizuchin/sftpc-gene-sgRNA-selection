"""把 CRISPOR 的 guides TSV 解析成结构化记录。

TSV 表头（第一行以 # 开头），实测两种格式：

  download=guides&format=tsv
      #guideId targetSeq mitSpecScore cfdSpecScore offtargetCount
      targetGenomeGeneLocus Doench '16-Score Moreno-Mateos-Score
      Doench-RuleSet3-Score Out-of-Frame-Score Lindel-Score GrafEtAlStatus grafType

  showAllScores=1&download=guides&format=tsv
      在上面基础上多出 Chari-Score Xu-Score Wu-Crispr-Score
      Doench '14-Score Wang-Score CCTop-Score

按「表头名字」取列，不按位置取，所以官方加列也不会解析错。
"""

from __future__ import annotations

from dataclasses import dataclass, asdict, fields


def _to_num(value: str):
    """'' / '-' / 'NA' -> None；其余转 int 或 float。"""
    if value is None:
        return None
    v = value.strip()
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
    """一条候选 guide。前 6 个字段来自基础表，后面是各评分模型。"""

    guide_id: str            # 如 "265rev"：位置 + 链
    target_seq: str          # 23 nt = 20 nt spacer + PAM
    mit_spec: int | None     # MIT 特异性分，0-100，越大越好
    cfd_spec: int | None     # CFD 特异性分，0-100，越大越好  <-- 数据2
    offtarget_count: int | None
    locus: str               # 如 "exon:HR"
    doench16: int | None     # Doench '16 效率分，0-100，越大越好 <-- 数据1
    moreno_mateos: int | None
    rs3: int | None          # Doench-RuleSet3，-200 ~ +200
    out_of_frame: int | None
    lindel: int | None
    extras: dict             # showAllScores 时的其它列

    # ---- 两个"数据"的便捷读取 ----

    @property
    def cut_position(self) -> int | None:
        """PAM 在输入序列上的位置（guideId 的数字部分）。"""
        digits = "".join(ch for ch in self.guide_id if ch.isdigit())
        return int(digits) if digits else None

    @property
    def strand(self) -> str:
        return "rev" if self.guide_id.lower().endswith("rev") else "fw"

    def as_row(self) -> dict:
        """导出用：把 extras 展平。"""
        row = {f.name: getattr(self, f.name) for f in fields(self) if f.name != "extras"}
        row.update(self.extras)
        return row


# TSV 表头名 -> GuideRecord 字段名
_COLUMN_MAP = {
    "guideId": "guide_id",
    "targetSeq": "target_seq",
    "mitSpecScore": "mit_spec",
    "cfdSpecScore": "cfd_spec",
    "offtargetCount": "offtarget_count",
    "targetGenomeGeneLocus": "locus",
    "Doench '16-Score": "doench16",
    "Moreno-Mateos-Score": "moreno_mateos",
    "Doench-RuleSet3-Score": "rs3",
    "Out-of-Frame-Score": "out_of_frame",
    "Lindel-Score": "lindel",
}

_NUMERIC_FIELDS = {
    "mit_spec", "cfd_spec", "offtarget_count", "doench16",
    "moreno_mateos", "rs3", "out_of_frame", "lindel",
}


def parse_guides_tsv(text: str) -> list[GuideRecord]:
    """解析 guides TSV 文本，返回 GuideRecord 列表（保持原顺序）。"""
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines:
        raise ValueError("TSV 内容为空")

    header = lines[0].lstrip("#").split("\t")
    header = [h.strip() for h in header]
    if "guideId" not in header:
        raise ValueError(
            "这不是 CRISPOR 的 guides TSV（表头缺少 guideId）；"
            f"实际表头：{header[:5]}"
        )

    records: list[GuideRecord] = []
    for line in lines[1:]:
        cells = line.split("\t")
        # 行比表头短时右侧补空
        if len(cells) < len(header):
            cells += [""] * (len(header) - len(cells))

        known: dict = {}
        extras: dict = {}
        for col, raw in zip(header, cells):
            field_name = _COLUMN_MAP.get(col)
            if field_name is None:
                extras[col] = raw.strip()
                continue
            known[field_name] = _to_num(raw) if field_name in _NUMERIC_FIELDS else raw.strip()

        records.append(
            GuideRecord(
                guide_id=known.get("guide_id", ""),
                target_seq=known.get("target_seq", ""),
                mit_spec=known.get("mit_spec"),
                cfd_spec=known.get("cfd_spec"),
                offtarget_count=known.get("offtarget_count"),
                locus=known.get("locus", ""),
                doench16=known.get("doench16"),
                moreno_mateos=known.get("moreno_mateos"),
                rs3=known.get("rs3"),
                out_of_frame=known.get("out_of_frame"),
                lindel=known.get("lindel"),
                extras=extras,
            )
        )
    return records


def records_to_tsv(records: list[GuideRecord]) -> str:
    """把记录写回 TSV（含 extras 列），用于落盘。"""
    if not records:
        return ""
    rows = [r.as_row() for r in records]
    columns: list[str] = []
    for row in rows:
        for key in row:
            if key not in columns:
                columns.append(key)
    out = ["\t".join(columns)]
    for row in rows:
        out.append("\t".join("" if row.get(c) is None else str(row.get(c)) for c in columns))
    return "\n".join(out) + "\n"
