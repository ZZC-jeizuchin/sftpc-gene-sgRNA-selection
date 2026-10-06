# 本项目由华南师范大学附属中学知识城校区高二一班同学制作
# Made by students of Class 1, Grade 11, The Affiliated High School of SCNU (Knowledge City Campus)
# Copyright (C) 2026 华南师范大学附属中学知识城校区高二一班
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for more
# details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""第三个数据（CRISPOR 版）：切割位点与 chr8:22,163,096 的距离。

和 CHOPCHOP 版**共用同一套坐标逻辑**（genome_utils），保证两边口径一致：
    CRISPOR 的 Position/Strand = **PAM 区间**在输入序列上的左端（anchor="pam"）
    CHOPCHOP 的 Genomic location = **23mer** 在基因组上的最左端（anchor="23mer"）
    fw/正链上两者相差 20 bp，所以必须区分，不能共用同一个偏移。

切点：正链 pam_start-3.5，反链 pam_start+5.5（SpCas9 切在 PAM 前 3 bp）
序列模式换算：基因组坐标 = genome_offset + 序列内位置
"""

from __future__ import annotations

from genome_utils import (
    DEFAULT_TARGET_POS, DISTANCE_BANDS, bands_json, cut_site, distance_band_score,
    distance_to, band_name, to_genomic, pam_start,
)

__all__ = [
    "DEFAULT_TARGET_POS", "DISTANCE_BANDS", "cut_site", "distance_to",
    "distance_band_score", "band_name", "build_distance_payload",
]


def build_distance_payload(records, target_pos: int = DEFAULT_TARGET_POS,
                           offset: int | None = None, chrom: str | None = None) -> dict:
    """把 CRISPOR 的 52/83... 条 guide 打包成第三个数据的 JSON 块。"""
    values = []
    for r in records:
        loc = r.cut_position                        # 输入序列内 1-based 坐标
        strand = r.strand                           # "fw" / "rev"
        genomic = to_genomic(loc, "seq", offset)    # 序列坐标 -> 基因组坐标
        # CRISPOR 的 Position/Strand = PAM 左端 -> anchor="pam"
        cut = cut_site(loc, strand, "seq", offset, anchor="pam")
        dist = distance_to(loc, strand, target_pos, "seq", offset, anchor="pam")
        score = distance_band_score(dist)
        values.append({
            "guide_id": r.guide_id,
            "strand": strand,
            "anchor": "pam",                        # CRISPOR 的号码是 PAM 左端
            "reported_position": loc,               # 网站原样给的数字
            "seq_position": loc,
            "genomic_location": genomic,            # 换算到基因组（= PAM 起点）
            "pam_start": genomic,                   # 归一化锚点：PAM 区间左端（基因组，可能为 None）
            "seq_pam_start": loc,                   # 同一个锚点，但在**输入序列**坐标系里（永远有值）
            "target_seq": r.target_seq,
            "cut_site": cut,
            "distance": dist,
            "score": score,
            "band": band_name(score),
            "mit_spec": r.mit_spec,
            "cfd_spec": r.cfd_spec,
            "doench16": r.doench16,
        })

    numeric = [v["distance"] for v in values if v["distance"] is not None]
    scores = [v["score"] for v in values if v["score"] is not None]
    return {
        "id": "distance_to_22163096",
        "name": "距离（切点 → chr8:22163096）",
        "source": "crispor",
        "source_url": "https://crispor.gi.ucsc.edu/",
        "source_column": "Position/Strand（切点由本程序计算）",
        "description": (
            "切割位点到 SFTPC c.218 / I73T 位点(hg38 chr8:22,163,096)的距离，"
            "按文档分档：<20bp=1.0, 20-50=0.9, 50-100=0.7, 100-200=0.4, >200=0.1。"
            "CRISPOR 给的是输入序列内坐标，需 genome_offset 才能换成基因组坐标。"
        ),
        "value_type": "float",
        "value_min": 0,
        "value_max": None,
        "direction": "lower_is_better",
        "unit": "bp",
        "target_pos": target_pos,
        "genome_chrom": chrom,
        "genome_offset": offset,
        "coordinate_mode": "sequence" if offset is not None else "sequence(no-offset)",
        "needs_genome_offset": offset is None,
        "bands": bands_json(),
        "count": len(values),
        "stats": {
            "min": min(numeric) if numeric else None,
            "max": max(numeric) if numeric else None,
            "mean": round(sum(numeric) / len(numeric), 2) if numeric else None,
            "missing": len(values) - len(numeric),
            "within_20bp": sum(1 for d in numeric if d < 20),
            "score_sum": round(sum(scores), 2) if scores else None,
        },
        "values": values,
    }
