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
"""第三个数据（CHOPCHOP 版）：切割位点与 chr8:22,163,096 的距离。

坐标逻辑全部来自 genome_utils（与 CRISPOR 版共用，保证口径一致）。
"""
"""第三个数据：**切割位点与 chr8:22,163,096 的距离**。

文档里的定义（四个数据的定义.txt 第 5 行）：
    距离 第二个网站 22163096取绝对值，然后打分
    (小于20为1.0, 20-50=0.9, 50-100=0.7, 100-200=0.4, 200以上=0.1)

为什么要自己算切点、而不是直接用 Genomic location：
    CHOPCHOP 的 "Genomic location" 是 **23 bp 靶位点（20nt spacer + PAM）在正链上的
    最左端坐标**，不是切割位点。SpCas9 切在 PAM 前 3 bp，所以：
        + 链：切点 = loc + 16.5
        - 链：切点 = loc + 5.5
    两条链差 11 bp。直接拿 loc 去减，会让 + 链和 - 链的 guide 拿到系统性不同的距离。
    （坐标约定已用 hg38 全序列对 161 条 guide 实测，全部通过。）

目标位点 22163096 是什么：
    hg38 chr8:22,163,096 = SFTPC NM_003018 的 **c.218**，密码子 73 的第 2 位
    （codon 73 = ATT = Ile），即最常见的 chILD 突变 **I73T (c.218T>C)** 的位点。
"""


from genome_utils import (
    DEFAULT_TARGET_POS, DISTANCE_BANDS, bands_json, cut_site, distance_band_score,
    distance_to, band_name, to_genomic, pam_start,
)

# 兼容旧名字
distance_of = distance_to
DEFAULT_TARGET_POS = DEFAULT_TARGET_POS


def build_distance_payload(records, target_pos: int = DEFAULT_TARGET_POS,
                           offset: int | None = None) -> dict:
    """把第三个数据打包成 JSON（与另外两个数据同构）。

    offset = 序列提交时的基因组 0-based 起点；用基因名提交时为 None（本来就是基因组坐标）。
    """
    values = []
    seq_mode = any((r.chrom or "").lower() == "seq" for r in records)
    for r in records:
        genomic = to_genomic(r.location, r.chrom, offset)
        cut = cut_site(r.location, r.strand, r.chrom, offset)
        dist = distance_of(r.location, r.strand, target_pos, r.chrom, offset)
        score = distance_band_score(dist)
        values.append({
            "guide_id": f"{r.location}{r.strand}" if r.location else r.target_seq[:12],
            "rank": r.rank,
            "chrom": r.chrom,
            "location": r.location,
            "anchor": "23mer",                      # CHOPCHOP 的号码是 23mer 最左端
            "reported_position": r.location,
            "genomic_location": genomic,            # 换算到基因组（= 23mer 最左端）
            "pam_start": pam_start(genomic, r.strand, "23mer"),  # 归一化锚点（基因组，可能为 None）
            "seq_pam_start": pam_start(r.location, r.strand, "23mer"),  # 同一锚点的序列坐标版
            "strand": r.strand,
            "target_seq": r.target_seq,
            "cut_site": cut,
            "distance": dist,
            "score": score,
            "band": band_name(score),
            # 对照：直接拿 23mer 左端当切点会算成什么（用来验证偏差）
            "distance_naive_nostrand": (
                None if genomic is None else abs(genomic - target_pos)
            ),
        })

    numeric = [v["distance"] for v in values if v["distance"] is not None]
    scores = [v["score"] for v in values if v["score"] is not None]
    return {
        "id": "distance_to_22163096",
        "name": "距离（切点 → chr8:22163096）",
        "source": "chopchop",
        "source_url": "https://chopchop.cbu.uib.no/",
        "source_column": "Genomic location + Strand（切点由本程序计算）",
        "description": (
            "切割位点到 SFTPC c.218 / I73T 位点(hg38 chr8:22,163,096)的距离，"
            "按文档分档打分：<20bp=1.0, 20-50=0.9, 50-100=0.7, 100-200=0.4, >200=0.1。"
            "切点 = Genomic location + (16.5 若 + 链 / 5.5 若 - 链)，因为 SpCas9 切在 PAM 前 3bp，"
            "而 CHOPCHOP 给的是 23mer 最左端坐标。"
        ),
        "value_type": "float",
        "value_min": 0,
        "value_max": None,
        "direction": "lower_is_better",
        "unit": "bp",
        "target_pos": target_pos,
        "bands": bands_json(),
        "genome_chrom": None,          # 序列模式建议同时给出染色体名（CHOPCHOP 只回 "seq"）
        "coordinate_mode": "sequence" if seq_mode else "genomic",
        "genome_offset": offset,
        "needs_genome_offset": bool(seq_mode and offset is None),
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
