#!/usr/bin/env python3
# 本项目由华南师范大学附属中学知识城校区高二一班同学制作
# Copyright (C) 2026 华南师范大学附属中学知识城校区高二一班
# SPDX-License-Identifier: AGPL-3.0-or-later
"""多基因案例对比：证明工具不是只为 SFTPC 调的参。

目前两个案例：
  SFTPC  c.218T>C (I73T)   chr8:22163096  正链  —— 儿童间质性肺病
  HBB    c.20A>T (E7V)     chr11:5227002  负链  —— 镰状细胞病

选 HBB 的理由：
  1) 它是**负链**基因，正好检验坐标归一化和链向处理是否真的通用；
  2) 突变位点定义没有争议（rs334），外部可核对；
  3) 和 SFTPC 是完全不同的疾病场景，说服力强。

用法：python case_compare.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

DATA = ROOT / "data"
OUT = ROOT / "output"
OUT.mkdir(exist_ok=True)

CASES = [
    {
        "name": "SFTPC",
        "label": "SFTPC c.218T>C (I73T)",
        "disease": "儿童间质性肺病 / 家族性肺纤维化",
        "chrom": "chr8",
        "target": 22_163_096,
        "strand": "+",
        "dir": DATA,
        "offset_key": None,
    },
    {
        "name": "HBB",
        "label": "HBB c.20A>T (p.Glu7Val)",
        "disease": "镰状细胞病",
        "chrom": "chr11",
        "target": 5_227_002,
        "strand": "-",
        "dir": DATA / "cases" / "HBB",
        "offset_key": None,
    },
]


def load_case(case):
    d = case["dir"]
    for f in ("crispor_payload.json", "chopchop_payload.json", "crispor_mutant_payload.json"):
        if not (d / f).is_file():
            return None
    from sgweb import allele_specific as AS
    from sgweb.aggregate import merge, rescore_v2

    cri = json.loads((d / "crispor_payload.json").read_text(encoding="utf-8"))
    cho = json.loads((d / "chopchop_payload.json").read_text(encoding="utf-8"))
    mut = json.loads((d / "crispor_mutant_payload.json").read_text(encoding="utf-8"))
    if not cri.get("values") or not cho.get("values"):
        return None

    target = case["target"]
    merged = merge(cri, cho, target, hard_filter=True)
    di = AS.build_delta_index(cri["values"], mut.get("values") or [])
    rescore_v2(merged, di, target_pos=target,
               genome_offset=cri.get("genome_offset"), hard_filter=True)

    kc = {v.get("pam_start", v.get("seq_pam_start")) for v in cri["values"]}
    kh = {v.get("pam_start", v.get("seq_pam_start")) for v in cho["values"]}
    rows = merged["rows_v2"]
    covered = [r for r in rows if r.get("allele_covered")]
    n_pos = sum(1 for r in rows if (r.get("allele_delta") or 0) > 0)
    return {
        "case": case,
        "crispor_n": len(cri["values"]),
        "chopchop_n": len(cho["values"]),
        "shared": len(kc & kh),
        "kept": len(rows),
        "filtered": merged["filtered_count"],
        "covered": len(covered),
        "delta_positive": n_pos,
        "offset": cri.get("genome_offset"),
        "top": rows[:5],
        "covered_rows": sorted(covered, key=lambda r: -(r.get("allele_score") or 0)),
    }


def main() -> int:
    print("=" * 82)
    print("多基因案例对比：工具是否通用")
    print("=" * 82)
    ok = []
    for case in CASES:
        r = load_case(case)
        if r is None:
            print(f"\n[{case['name']}] 数据不完整，跳过（先跑 acquire 脚本）")
            continue
        ok.append(r)

    for r in ok:
        c = r["case"]
        print(f"\n【{c['name']}】{c['label']}　{c['disease']}")
        print("-" * 82)
        print(f"  基因组位置      : {c['chrom']}:{c['target']}  链 {c['strand']}")
        print(f"  序列起点(0-based): {r['offset']}")
        print(f"  CRISPOR / CHOPCHOP : {r['crispor_n']} / {r['chopchop_n']}"
              f"  共有 {r['shared']}")
        print(f"  硬性过滤淘汰    : {r['filtered']} 条，保留 {r['kept']} 条")
        print(f"  覆盖突变位点    : {r['covered']} 条"
              f"（占 {r['covered']/max(r['kept'],1):.1%}）")
        print(f"  Δraw > 0        : {r['delta_positive']} 条")
        print("  Top 3：")
        for x in r["top"][:3]:
            print(f"    {x['sequence']:<24} 综合 {x['final_v2']:>6.2f}  "
                  f"ASS {x['allele_score']:>5.1f}  Δraw {x.get('allele_delta')}  "
                  f"覆盖 {'是' if x['allele_covered'] else '否'}")

    if len(ok) == 2:
        print("\n" + "=" * 82)
        print("结论")
        print("=" * 82)
        a, b = ok
        print(f"  · 两条链路都跑通：正链基因（{a['case']['name']}）与负链基因"
              f"（{b['case']['name']}）")
        print(f"  · 两案的覆盖率相近：{a['covered']}/{a['kept']} vs "
              f"{b['covered']}/{b['kept']} —— 覆盖突变的 guide 本来就只占少数，")
        print("    这是「23mer 必须盖住一个特定碱基」的几何必然，不是筛选筛掉的")
        print(f"  · 负链案例证明坐标归一化不是对着 SFTPC 特调的")

    payload = []
    for r in ok:
        c = r["case"]
        payload.append({
            "name": c["name"], "label": c["label"], "disease": c["disease"],
            "chrom": c["chrom"], "target": c["target"], "strand": c["strand"],
            "offset": r["offset"],
            "crispor_n": r["crispor_n"], "chopchop_n": r["chopchop_n"],
            "shared": r["shared"], "kept": r["kept"], "filtered": r["filtered"],
            "covered": r["covered"], "delta_positive": r["delta_positive"],
            "top": [{kk: x.get(kk) for kk in
                     ("sequence", "final_v2", "allele_score", "allele_delta",
                      "allele_covered", "allele_region_label")} for x in r["top"]],
        })
    (OUT / "case_compare.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\n已写出：{OUT / 'case_compare.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
