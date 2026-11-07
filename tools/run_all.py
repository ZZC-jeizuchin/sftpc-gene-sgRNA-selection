#!/usr/bin/env python3
# 本项目由华南师范大学附属中学知识城校区高二一班同学制作
# Copyright (C) 2026 华南师范大学附属中学知识城校区高二一班
# SPDX-License-Identifier: AGPL-3.0-or-later
"""一键复现：把整个研究流程从头到尾跑一遍。

    python tools/run_all.py              # 全部跑（有缓存，通常几分钟）
    python tools/run_all.py --skip-net   # 只跑离线分析，不联网取数
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PY = sys.executable

# (相对 root 的脚本路径, 说明)
NET_STEPS = [
    ("tools/acquire.py", "抓 CRISPOR(野生型) + CHOPCHOP"),
    ("tools/acquire_mutant.py", "抓 CRISPOR(突变型)"),
]
STEPS = [
    ("analysis/score_v2.py", "四维评分 + v2 排名表"),
    ("analysis/weight_analysis.py", "CRITIC 客观赋权 + 消融 + 稳定性"),
    ("analysis/validate_external.py", "CFD 外部佐证"),
    ("analysis/simulate_experiment.py", "体外切割模拟 + 凝胶图"),
    ("analysis/case_compare.py", "多基因案例对比"),
]


def run(script: str) -> tuple[bool, float]:
    t0 = time.time()
    r = subprocess.run([PY, script], cwd=ROOT,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    dt = time.time() - t0
    for line in r.stdout.decode("utf-8", "replace").strip().splitlines()[-3:]:
        print("      " + line)
    return r.returncode == 0, dt


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-net", action="store_true", help="跳过联网取数")
    args = ap.parse_args()

    print("=" * 70)
    print("sgRNA 智能筛选工具 —— 全流程复现")
    print("=" * 70)

    steps = STEPS if args.skip_net else NET_STEPS + STEPS

    failed = []
    for script, desc in steps:
        if not (ROOT / script).is_file():
            print(f"[跳过] {script}（文件不存在）")
            continue
        print(f"\n▶ {desc}")
        print(f"  $ python {script}")
        ok, dt = run(script)
        print(f"  {'✔ 完成' if ok else '✘ 失败'}  ({dt:.1f}s)")
        if not ok:
            failed.append(script)

    print("\n" + "=" * 70)
    if failed:
        print("以下步骤失败：")
        for f in failed:
            print("  - " + f)
        return 1
    print("全部完成。产物在 output/ 下：")
    for p in sorted((ROOT / "output").glob("*")):
        if p.is_file():
            print(f"  - output/{p.name}  ({p.stat().st_size/1024:.1f} KB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
