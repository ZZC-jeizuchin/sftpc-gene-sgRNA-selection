#!/usr/bin/env python3
# 本项目由华南师范大学附属中学知识城校区高二一班同学制作
# Copyright (C) 2026 华南师范大学附属中学知识城校区高二一班
# SPDX-License-Identifier: AGPL-3.0-or-later
"""通用版取数：任意基因 / 任意单碱基突变。

用来回答「你的工具是不是只为 SFTPC 调的参」。
输入：染色体区间 + 突变坐标 + 参考/突变碱基；
输出：data/cases/<name>/ 下三份 payload（CRISPOR 野生型、CRISPOR 突变型、CHOPCHOP）。

用法：
  python acquire_case.py --name HBB --chrom chr11 --start 5226300 --end 5227500 \
      --target 5227002 --ref T --alt A
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

UCSC = "https://api.genome.ucsc.edu/getData/sequence"


def say(m):
    print(f"[case] {m}", flush=True)


def fetch_seq(genome: str, chrom: str, start: int, end: int,
              cache: Path) -> str:
    cache.mkdir(parents=True, exist_ok=True)
    p = cache / f"{genome}_{chrom}_{start}_{end}.txt"
    if p.is_file():
        return p.read_text(encoding="utf-8").strip().upper()
    url = f"{UCSC}?" + urllib.parse.urlencode(
        {"genome": genome, "chrom": chrom, "start": start, "end": end})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                dna = json.loads(r.read().decode())["dna"].upper()
            p.write_text(dna, encoding="utf-8")
            return dna
        except Exception as exc:  # noqa: BLE001
            say(f"UCSC 失败({attempt+1}/3): {exc}")
            time.sleep(2 ** attempt)
    raise RuntimeError("拉取序列失败")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--genome", default="hg38")
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--start", type=int, required=True, help="0-based 起点")
    ap.add_argument("--end", type=int, required=True, help="0-based 终点（不含）")
    ap.add_argument("--target", type=int, required=True, help="突变位点，1-based")
    ap.add_argument("--ref", required=True, help="正链参考碱基")
    ap.add_argument("--alt", required=True, help="正链突变碱基")
    ap.add_argument("--pam", default="NGG")
    args = ap.parse_args()

    out_dir = ROOT / "data" / "cases" / args.name
    out_dir.mkdir(parents=True, exist_ok=True)

    seq = fetch_seq(args.genome, args.chrom, args.start, args.end,
                    ROOT / ".genome_cache" / "case")
    idx = (args.target - 1) - args.start
    if not (0 <= idx < len(seq)):
        say(f"突变位点不在序列内（idx={idx}, len={len(seq)}）")
        return 2
    base = seq[idx]
    say(f"序列 {len(seq)} bp  chr{args.chrom}:{args.start}-{args.end}  "
        f"突变位点 {args.chrom}:{args.target} 碱基={base}")
    if base != args.ref.upper():
        say(f"⚠ 参考碱基不符：期望 {args.ref}，实得 {base}。"
            f"请确认链方向（负链基因的 c. 编号与正链碱基互为互补）")
        return 2

    mut = seq[:idx] + args.alt.upper() + seq[idx + 1:]
    (out_dir / f"{args.name}_wt.txt").write_text(seq, encoding="utf-8")
    (out_dir / f"{args.name}_mut.txt").write_text(mut, encoding="utf-8")

    common = dict(genome=args.genome, pam=args.pam, target_pos=args.target)
    loc = dict(genome_offset=args.start, genome_chrom=args.chrom, auto_locate=False)

    t0 = time.time()

    # 1) CRISPOR 野生型
    try:
        from crispor_scraper import CrisporInput, fetch as cfetch
        say("CRISPOR(野生型) 提交中…")
        r1 = cfetch(CrisporInput(sequence=seq, name=f"{args.name}_wt",
                                 poll_interval=3.0, timeout=1200.0, **common, **loc),
                    on_progress=say)
        (out_dir / "crispor_payload.json").write_text(
            json.dumps(r1.distance_payload(), ensure_ascii=False, indent=2), encoding="utf-8")
        say(f"CRISPOR(野生型)：{len(r1.records)} 条")
    except Exception as exc:  # noqa: BLE001
        say(f"CRISPOR(野生型) 失败：{exc!r}")

    # 2) CRISPOR 突变型
    try:
        from crispor_scraper import CrisporInput, fetch as cfetch
        say("CRISPOR(突变型) 提交中…")
        r2 = cfetch(CrisporInput(sequence=mut, name=f"{args.name}_mut",
                                 poll_interval=3.0, timeout=1200.0, **common, **loc),
                    on_progress=say)
        (out_dir / "crispor_mutant_payload.json").write_text(
            json.dumps(r2.distance_payload(), ensure_ascii=False, indent=2), encoding="utf-8")
        say(f"CRISPOR(突变型)：{len(r2.records)} 条")
    except Exception as exc:  # noqa: BLE001
        say(f"CRISPOR(突变型) 失败：{exc!r}")

    # 3) CHOPCHOP
    try:
        from chopchop_scraper import ChopchopInput, fetch as hfetch
        say("CHOPCHOP 提交中…")
        r3 = hfetch(ChopchopInput(sequence=seq, genome=args.genome,
                                  sequence_name=f"{args.name}_wt",
                                  poll_interval=3.0, timeout=1200.0,
                                  genome_offset=args.start,
                                  genome_chrom=args.chrom, auto_locate=False),
                    on_progress=say)
        (out_dir / "chopchop_payload.json").write_text(
            json.dumps(r3.distance_payload(args.target), ensure_ascii=False, indent=2),
            encoding="utf-8")
        say(f"CHOPCHOP：{len(r3.records)} 条")
    except Exception as exc:  # noqa: BLE001
        say(f"CHOPCHOP 失败：{exc!r}")

    meta = {"name": args.name, "genome": args.genome, "chrom": args.chrom,
            "start": args.start, "end": args.end, "target": args.target,
            "ref": args.ref, "alt": args.alt, "pam": args.pam,
            "offset_in_seq": idx}
    (out_dir / "case.json").write_text(json.dumps(meta, ensure_ascii=False, indent=2),
                                       encoding="utf-8")
    say(f"完成，用时 {time.time() - t0:.1f}s -> {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
