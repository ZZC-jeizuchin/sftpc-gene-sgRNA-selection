#!/usr/bin/env python3
# 本项目由华南师范大学附属中学知识城校区高二一班同学制作
# Copyright (C) 2026 华南师范大学附属中学知识城校区高二一班
# SPDX-License-Identifier: AGPL-3.0-or-later
"""抓取两个网站的真实原始数据，落盘为 JSON，供离线重算评分。

只负责「采集」，不负责「评分」——评分逻辑改动后可以直接重跑 score.py，
不必再等 CRISPOR 排队。
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from genome_utils import DEFAULT_TARGET_POS  # noqa: E402

OUT = ROOT / "data"
OUT.mkdir(exist_ok=True)

TARGET_POS = DEFAULT_TARGET_POS  # 22_163_096


def say(msg):
    print(f"[acquire] {msg}", flush=True)


def main() -> int:
    seq = (ROOT / "SFTPC_测试序列_hg38_2213bp.txt").read_text(encoding="utf-8")
    seq = "".join(seq.split()).upper()
    say(f"序列 {len(seq)} bp，目标位点 hg38 chr8:{TARGET_POS}")

    t0 = time.time()

    # ---- 网站1：CRISPOR -> Doench'16 + CFD Spec ----
    crispor_payload = None
    try:
        say("CRISPOR: 提交中…")
        from crispor_scraper import CrisporInput, fetch as crispor_fetch
        r1 = crispor_fetch(CrisporInput(
            sequence=seq, genome="hg38", pam="NGG", name="sftpc_full",
            target_pos=TARGET_POS, poll_interval=3.0, timeout=1200.0,
        ), on_progress=say)
        crispor_payload = r1.distance_payload()
        say(f"CRISPOR: {len(r1.records)} 条 guide, batch={r1.batch_id}")
    except Exception as exc:  # noqa: BLE001
        say(f"CRISPOR 失败: {exc!r}")

    # ---- 网站2：CHOPCHOP -> 切点 / 距离 ----
    chopchop_payload = None
    try:
        say("CHOPCHOP: 提交中…")
        from chopchop_scraper import ChopchopInput, fetch as chopchop_fetch
        r2 = chopchop_fetch(ChopchopInput(
            sequence=seq, genome="hg38", sequence_name="sftpc_full",
            poll_interval=3.0, timeout=1200.0,
        ), on_progress=say)
        chopchop_payload = r2.distance_payload(TARGET_POS)
        say(f"CHOPCHOP: {len(r2.records)} 条 guide, job={r2.job_id}")
    except Exception as exc:  # noqa: BLE001
        say(f"CHOPCHOP 失败: {exc!r}")

    (OUT / "crispor_payload.json").write_text(
        json.dumps(crispor_payload, ensure_ascii=False, indent=2), encoding="utf-8")
    (OUT / "chopchop_payload.json").write_text(
        json.dumps(chopchop_payload, ensure_ascii=False, indent=2), encoding="utf-8")

    say(f"完成，用时 {time.time() - t0:.1f}s -> {OUT}")
    if crispor_payload:
        say(f"  crispor values={len(crispor_payload.get('values') or [])}")
    if chopchop_payload:
        say(f"  chopchop values={len(chopchop_payload.get('values') or [])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
