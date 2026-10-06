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
"""把 config / cache / client / parser 串成一条高层流水线。

缓存优先级（越靠前越省时间）：
    1. 语义缓存里已有该输入的 TSV  -> 直接返回，0 次网络请求
    2. 语义缓存里已有 batchId      -> 跳过提交，直接轮询旧批次
    3. 什么都没有                  -> 正常提交
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from .cache import CrisporCache, DEFAULT_CACHE_DIR
from .client import BASE_URL, CrisporClient
from .config import CrisporInput
from .parser import GuideRecord, parse_guides_tsv, records_to_tsv


@dataclass
class CrisporResult:
    """一次查询的完整结果。"""

    batch_id: str
    input: CrisporInput
    records: list[GuideRecord]
    raw_tsv: str = ""
    from_cache: bool = False      # True = TSV 直接来自磁盘缓存
    reused_batch: bool = False    # True = 复用了旧 batchId，但重新下载了结果

    @property
    def url(self) -> str:
        return f"{BASE_URL}?batchId={self.batch_id}"

    def distance_payload(self) -> dict:
        """第三个数据（切点 -> chr8:22163096 的距离）。"""
        from .distance import build_distance_payload
        return build_distance_payload(self.records, self.input.target_pos,
                                      offset=self.input.genome_offset,
                                      chrom=self.input.genome_chrom)

    def two_data(self) -> list[dict]:
        """只保留文档里定义的两个数据：Doench '16 与 CFD Spec. score。"""
        return [
            {
                "guide_id": r.guide_id,
                "strand": r.strand,
                "cut_position": r.cut_position,
                "target_seq": r.target_seq,
                "doench16": r.doench16,
                "cfd_spec": r.cfd_spec,
                "mit_spec": r.mit_spec,
                "offtarget_count": r.offtarget_count,
            }
            for r in self.records
        ]

    def save(self, out_dir: str | Path = ".") -> dict:
        """落盘：原始 TSV + 提取后的两数据 TSV。返回写出的文件路径。"""
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        base = f"crispor_{self.input.name}_{stamp}"

        raw_path = out_dir / f"{base}_all.tsv"
        raw_path.write_text(self.raw_tsv or records_to_tsv(self.records), encoding="utf-8")

        two_path = out_dir / f"{base}_two_data.tsv"
        rows = self.two_data()
        cols = list(rows[0].keys()) if rows else [
            "guide_id", "strand", "cut_position", "target_seq",
            "doench16", "cfd_spec", "mit_spec", "offtarget_count",
        ]
        lines = ["\t".join(cols)]
        for row in rows:
            lines.append("\t".join("" if row.get(c) is None else str(row.get(c)) for c in cols))
        two_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

        return {"all_scores": str(raw_path), "two_data": str(two_path)}


def fetch(
    inp: CrisporInput,
    client: CrisporClient | None = None,
    on_progress=None,
    use_cache: bool | None = None,
    cache_dir: str | Path = DEFAULT_CACHE_DIR,
    from_batch: str | None = None,
) -> CrisporResult:
    """执行一次完整抓取：缓存 -> 提交 -> 轮询 -> 下载 -> 解析。

    use_cache / from_batch 为 None 时取 CrisporInput 上的同名字段。
    """
    def log(msg: str) -> None:
        if on_progress:
            on_progress(msg)

    inp.validate()

    # 0. 只给碱基序列时，自动定位回基因组（这样 22163096 这个基因组目标位点才用得上）
    locate_info = None
    if inp.sequence and inp.genome_offset is None and inp.auto_locate:
        from genome_utils import DEFAULT_LOCUS, DEFAULT_LOCUS_CACHE, locate_sequence
        loc_kw = {}
        if inp.locus_chrom:
            loc_kw = {"chrom": inp.locus_chrom,
                      "start": inp.locus_start if inp.locus_start is not None else DEFAULT_LOCUS["start"],
                      "end": inp.locus_end if inp.locus_end is not None else DEFAULT_LOCUS["end"]}
        locate_info = locate_sequence(inp.sequence, cache_dir=DEFAULT_LOCUS_CACHE, **loc_kw)
        if locate_info.get("found"):
            inp.genome_offset = locate_info["offset"]
            inp.genome_chrom = locate_info.get("chrom")
            log(f"自动定位：序列在 {locate_info['chrom']} 上，genome_offset="
                f"{locate_info['offset']}（{locate_info['note']}）")
        else:
            log(f"⚠ 自动定位失败：{locate_info.get('note')}；距离分将算不出来")

    if use_cache is None:
        use_cache = inp.use_cache
    if from_batch is None:
        from_batch = inp.batch_id
    cache = CrisporCache(cache_dir, enabled=use_cache)
    client = client or CrisporClient(timeout=60, retries=inp.retries, cache=cache)
    client.cache = cache

    key = CrisporCache.key_from_input(inp)

    # ---- 1. 结果已经在缓存里 ----
    if use_cache:
        cached_tsv = cache.get_tsv(key)
        if cached_tsv:
            batch_id = cache.get_batch_id(key) or "cached"
            log(f"缓存命中：直接读取本地结果（batchId={batch_id}），0 次网络请求")
            records = parse_guides_tsv(cached_tsv)
            log(f"解析到 {len(records)} 条候选 guide")
            return CrisporResult(batch_id, inp, records, cached_tsv,
                                 from_cache=True, reused_batch=True)
    else:
        log("已禁用缓存，全部重新请求")

    if from_batch:
        log(f"使用指定的 batchId = {from_batch}（跳过提交）")
        batch_id = from_batch
        cache.put_batch_id(key, batch_id, describe={"input": inp.describe()})
    else:
        log("Step1/2/3 输入：\n" + inp.describe())
        log("提交任务…")
        batch_id, reused = client.submit_cached(inp, on_progress=log)
        if not reused:
            log(f"已排队，batchId = {batch_id}")

    # 关键：不轮询 265 KB 的结果页，直接重试 6 KB 的 TSV 下载端点
    log("等待计算完成（直接探测 TSV 下载端点，绕开结果页刷新）…")
    raw_tsv = client.fetch_guides_when_ready(
        batch_id,
        all_scores=inp.show_all_scores,
        poll_interval=inp.poll_interval,
        timeout=inp.timeout,
        on_progress=log,
    )
    cache.put_tsv(key, batch_id, raw_tsv, describe={"input": inp.describe()})

    records = parse_guides_tsv(raw_tsv)
    log(f"解析到 {len(records)} 条候选 guide，已写入缓存")
    return CrisporResult(batch_id, inp, records, raw_tsv, from_cache=False, reused_batch=False)
