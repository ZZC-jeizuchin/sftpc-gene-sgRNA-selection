#!/usr/bin/env python3
# 本项目由华南师范大学附属中学知识城校区高二一班同学制作
# Copyright (C) 2026 华南师范大学附属中学知识城校区高二一班
# SPDX-License-Identifier: AGPL-3.0-or-later
"""把「只差一个碱基」的突变型序列也提交给 CRISPOR，用于算 ΔScore。

输入序列与第一次完全相同（同长度、同起点），只把 chr8:22163096 的 T 换成 C。
两次运行的 guide 位置一一对应，按 (链向, seq_pam_start) 配对即可求差。
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from genome_utils import DEFAULT_TARGET_POS  # noqa: E402
from sgweb.allele_specific import (  # noqa: E402
    MUTATION, build_mutant, check_ref, mutation_index,
)

DATA = ROOT / "data"
SEQ_FILE = ROOT / "SFTPC_测试序列_hg38_2213bp.txt"


def say(msg):
    print(f"[mutant] {msg}", flush=True)


def main() -> int:
    wt = "".join(SEQ_FILE.read_text(encoding="utf-8").split()).upper()

    # 从野生型那次运行里拿到准确的基因组起点，保证两次坐标完全一致
    cri = json.loads((DATA / "crispor_payload.json").read_text(encoding="utf-8"))
    offset = cri.get("genome_offset")
    chrom = cri.get("genome_chrom")
    if offset is None:
        say("野生型 payload 里没有 genome_offset，无法构造突变型序列")
        return 1
    say(f"genome_offset={offset}  chrom={chrom}")

    idx = mutation_index(wt, offset, DEFAULT_TARGET_POS)
    if not check_ref(wt, idx):
        say(f"参考碱基不符：第 {idx} 位是 {wt[idx]}，突变记录要求 {MUTATION['ref']}")
        return 1
    mut = build_mutant(wt, idx)
    say(f"突变位点：序列第 {idx} 位（基因组 1-based {offset + idx + 1}），"
        f"{wt[idx]} -> {mut[idx]}；序列长度 {len(mut)}")

    # 突变型序列落盘（答辩/复核都用得上）
    (DATA / "SFTPC_I73T_mutant_2213bp.txt").write_text(mut, encoding="utf-8")

    t0 = time.time()
    say("CRISPOR(突变型): 提交中…")
    from crispor_scraper import CrisporInput, fetch as crispor_fetch
    r = crispor_fetch(CrisporInput(
        sequence=mut, genome="hg38", pam="NGG", name="sftpc_I73T",
        target_pos=DEFAULT_TARGET_POS, poll_interval=3.0, timeout=1200.0,
        genome_offset=offset, genome_chrom=chrom or "chr8", auto_locate=False,
    ), on_progress=say)
    payload = r.distance_payload()
    (DATA / "crispor_mutant_payload.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    say(f"CRISPOR(突变型): {len(r.records)} 条 guide, batch={r.batch_id}, "
        f"用时 {time.time() - t0:.1f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
