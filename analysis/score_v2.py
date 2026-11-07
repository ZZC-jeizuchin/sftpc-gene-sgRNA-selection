#!/usr/bin/env python3
# 本项目由华南师范大学附属中学知识城校区高二一班同学制作
# Copyright (C) 2026 华南师范大学附属中学知识城校区高二一班
# SPDX-License-Identifier: AGPL-3.0-or-later
"""v2 评分：把「等位基因特异性」作为第四维，四个维度统一 0-100 加权排序。

    v1  0.412×效率 + 0.353×脱靶 + 0.235×距离                 (单序列)
    v2  0.35×等位基因特异性 + 0.30×效率 + 0.25×脱靶 + 0.10×距离  (双序列)

数据全部来自真实抓取（data/），不联网也能重跑。
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from genome_utils import DEFAULT_TARGET_POS  # noqa: E402
from sgweb import allele_specific as AS  # noqa: E402
from sgweb.aggregate import merge, rescore_v2  # noqa: E402

DATA = ROOT / "data"
OUT = ROOT / "output"
OUT.mkdir(exist_ok=True)


def load(name):
    p = DATA / name
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None


def fmt(v, nd=2):
    if v is None:
        return "-"
    if isinstance(v, float):
        return f"{v:.{nd}f}"
    return str(v)


def main() -> int:
    cri = load("crispor_payload.json")
    cho = load("chopchop_payload.json")
    cri_mut = load("crispor_mutant_payload.json")
    if not cri or not cho:
        print("缺少 data/ 里的抓取结果，先跑 acquire.py", file=sys.stderr)
        return 1
    if cri_mut is None:
        print("缺少 data/crispor_mutant_payload.json，先跑 acquire_mutant.py", file=sys.stderr)
        return 1

    target = DEFAULT_TARGET_POS
    offset = cri.get("genome_offset")

    # 1) v1 合并（候选交集 / 坐标 / 距离 / 硬性过滤）
    merged = merge(cri, cho, target, hard_filter=True)

    # 2) 双序列 ΔScore
    delta_index = AS.build_delta_index(
        cri["values"], cri_mut["values"], score_key="doench16")
    n_both = sum(1 for v in delta_index.values() if v["delta_source"] == "both")
    n_gained = sum(1 for v in delta_index.values() if v["delta_source"].startswith("mut_only"))
    n_lost = sum(1 for v in delta_index.values() if v["delta_source"].startswith("wt_only"))

    # 3) v2 四维评分
    rescore_v2(merged, delta_index, target_pos=target, genome_offset=offset,
               hard_filter=True)

    rows = merged["rows_v2"]
    covered = [r for r in rows if r.get("allele_covered")]
    nonzero_delta = [r for r in rows
                     if (r.get("allele_delta") or 0) > 0]

    print("=" * 78)
    print("数据规模")
    print("=" * 78)
    kc = {v.get("pam_start", v.get("seq_pam_start")) for v in cri["values"]}
    kh = {v.get("pam_start", v.get("seq_pam_start")) for v in cho["values"]}
    print(f"CRISPOR guide 数            : {len(cri['values'])}")
    print(f"CHOPCHOP guide 数           : {len(cho['values'])}")
    print(f"CRISPOR(突变型) guide 数    : {len(cri_mut['values'])}")
    print(f"两站共有（按 pam_start 配对）: {len(kc & kh)}  "
          f"（仅 CRISPOR {len(kc - kh)}，仅 CHOPCHOP {len(kh - kc)}）")
    print(f"硬性过滤淘汰                : {merged['filtered_count']} 条  "
          f"（{merged['filtered_reasons']}）")
    print(f"进入 v2 评分的 guide        : {len(rows)}  "
          f"（两站都有 {merged['parts']['part1']} + 仅 CRISPOR {merged['parts']['part2']}）")
    print()
    print(f"双序列配对：两站都出现 {n_both} 条；"
          f"只突变型(新PAM) {n_gained} 条；只野生型(PAM丢失) {n_lost} 条")
    print(f"23mer 覆盖突变位点的 guide  : {len(covered)} 条")
    print(f"ΔScore > 0 的 guide         : {len(nonzero_delta)} 条")
    print()

    # ---- 1. v1 第一名 ----
    v1_top = max((r for r in rows if r.get("rank")), key=lambda r: r["rank"] * -1,
                 default=None)
    if v1_top:
        print("=" * 78)
        print("v1（现行三维修）排名第一")
        print("=" * 78)
        print(f"  {v1_top['sequence']}  rank={v1_top['rank']}  "
              f"final={fmt(v1_top['final_score'], 4)}  distance={fmt(v1_top['distance'],0)}bp")
        print(f"  v2 综合分 = {fmt(v1_top['final_v2'])}  等级 = {v1_top['grade']}  "
              f"v2 排名 = {v1_top['rank_v2']}")
        print(f"  ΔScore = {fmt(v1_top['allele_delta'])}  "
              f"覆盖突变 = {v1_top['allele_covered']}  "
              f"位置 = {v1_top['allele_region_label']}")
        print()

    # ---- 2. v2 top 15 ----
    print("=" * 110)
    print("v2 排名 Top 15（四维 0-100 加权）")
    print("=" * 110)
    hdr = (f"{'#':>3}  {'sgRNA(20nt+PAM)':<24} {'链':<2} {'距离':>5} "
           f"{'Δraw':>5} {'ΔScore':>6} {'位置':<9} {'等位':>6} {'效率':>5} {'脱靶':>5} {'距离分':>5} "
           f"{'综合':>6} {'等级':<6}")
    print(hdr)
    print("-" * 110)
    for r in rows[:15]:
        print(f"{r['rank_v2']:>3}  {r['sequence']:<24} {r['strand']:<2} "
              f"{fmt(r['distance'],0):>5} "
              f"{fmt(r['allele_delta'],1):>5} {fmt(r.get('allele_delta_100'),1):>6} "
              f"{r['allele_region_label']:<9} "
              f"{fmt(r['allele_score'],1):>6} {fmt(r['eff_score'],0):>5} "
              f"{fmt(r['cfd_score'],0):>5} {fmt(r['dist_score'],0):>5} "
              f"{fmt(r['final_v2'],2):>6} {r['grade']:<6}")
    print()

    # ---- 3. 只覆盖突变的那些 ----
    print("=" * 110)
    print(f"覆盖突变位点的全部 guide（{len(covered)} 条）")
    print("=" * 110)
    for r in sorted(covered, key=lambda x: -(x["final_v2"] or 0)):
        print(f"  {r['sequence']:<24} 距离{fmt(r['distance'],0):>4}bp  "
              f"Δraw {fmt(r['allele_delta'],1):>5}  ΔScore {fmt(r.get('allele_delta_100'),1):>5}  "
              f"位置 {r['allele_region_label']:<8}(第{r['allele_pos_from_pam']}位)  "
              f"ASS {fmt(r['allele_score'],1):>5}  综合 {fmt(r['final_v2'],2)}")
    print()

    # ---- 4. 模型说明 ----
    print("=" * 78)
    print("v2 模型")
    print("=" * 78)
    print(f"  综合分 = 0.35×等位基因特异性 + 0.30×切割效率 + 0.25×脱靶特异性 + 0.10×突变距离")
    print(f"  等位基因特异性 = 0.6×ΔScore + 0.4×位置分")
    print(f"  ΔScore = 100 × clamp(Doench'16 双序列差, 0, {AS.DELTA_FULL_SCALE}) / {AS.DELTA_FULL_SCALE}")
    print(f"  位置分：PAM 区 100 / Seed 区 80 / 中间区 40 / PAM 远端 10 / 不覆盖 0")
    print(f"  降权：ΔScore < {AS.DOWNGRADE_DELTA} 时 ×{AS.DOWNGRADE_FACTOR}")
    print(f"  距离维度：≤20bp=100 / ≤50=90 / ≤100=70 / ≤200=40 / >200=10，"
          f"覆盖突变再 +10（封顶 100）")
    print(f"  等级分布：{merged['grade_counts']}")
    print()

    # ---- 5. 权重敏感性：等位基因特异性占多大比重才合理 ----
    print("=" * 110)
    print("权重敏感性：等位基因特异性（ASS）占比变化对 Top5 的影响")
    print("=" * 110)
    print(f"{'ASS权重':>7}  {'Top1':<26} {'分':>6}  {'Top2':<26} {'分':>6}  "
          f"{'≥70分条数':>9}  {'覆盖突变进Top10':>14}")
    # 其余三维按 30:25:10 的比例分配剩下的权重
    rest = {"doench16": 30 / 65.0, "cfd_spec": 25 / 65.0, "distance": 10 / 65.0}
    for w in (0.0, 0.10, 0.20, 0.35):
        scored = []
        for r in rows:
            others = sum(rest[k] * float(r[f]) for k, f in
                         (("doench16", "eff_score"), ("cfd_spec", "cfd_score"),
                          ("distance", "dist_score")) if r.get(f) is not None)
            val = w * float(r["allele_score"]) + (1 - w) * others
            scored.append((val, r))
        scored.sort(key=lambda t: -t[0])
        top10 = [r for _, r in scored[:10]]
        n70 = sum(1 for v, _ in scored if v >= 70)
        cov10 = sum(1 for r in top10 if r.get("allele_covered"))
        print(f"{w*100:>6.0f}%  {scored[0][1]['sequence']:<26} {scored[0][0]:>6.2f}  "
              f"{scored[1][1]['sequence']:<26} {scored[1][0]:>6.2f}  "
              f"{n70:>9}  {cov10:>13}/10")
    print()
    print("  说明：ASS 权重越高，能区分两条等位基因（覆盖突变）的 guide 越占优；")
    print("        ASS = 0 就退化成原来的三维修（只是距离权重被压到 10%）。")
    print()

    # ---- 写文件 ----
    write_xlsx(rows, merged, delta_index)
    (OUT / "ranking_v2.json").write_text(
        json.dumps({"model": "v2", "weights": merged["weights_v2"],
                    "allele_model": merged["allele_model"],
                    "count": len(rows),
                    "rows": rows}, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"已写出：{OUT / 'ranking_v2.xlsx'}")
    print(f"已写出：{OUT / 'ranking_v2.json'}")
    return 0


def write_xlsx(rows, merged, delta_index):
    from openpyxl import Workbook
    from openpyxl.styles import Alignment, Font, PatternFill
    from openpyxl.utils import get_column_letter

    wb = Workbook()
    head_fill = PatternFill("solid", fgColor="1F4E79")
    head_font = Font(color="FFFFFF", bold=True)
    hot_fill = PatternFill("solid", fgColor="E8F3E8")

    def style(ws, widths):
        for c in range(1, ws.max_column + 1):
            cell = ws.cell(row=1, column=c)
            cell.fill = head_fill
            cell.font = head_font
            cell.alignment = Alignment(horizontal="center", vertical="center",
                                       wrap_text=True)
        for i, w in enumerate(widths, 1):
            ws.column_dimensions[get_column_letter(i)].width = w
        ws.freeze_panes = "A2"

    # Sheet1 排名
    ws = wb.active
    ws.title = "v2排名"
    cols = ["排名", "sgRNA(20nt+PAM)", "spacer(20nt)", "PAM", "链", "切点(基因组)",
            "距突变(bp)", "覆盖突变", "突变位置", "距PAM第几位",
            "Doench16(野生型)", "Doench16(突变型)", "Δraw(分)", "ΔScore(0-100)",
            "等位基因特异性", "切割效率", "脱靶特异性", "突变距离分",
            "综合分", "等级", "v1综合分", "v1排名"]
    ws.append(cols)
    for r in rows:
        seq = r["sequence"] or ""
        ws.append([
            r["rank_v2"], seq, seq[:20], seq[20:23] if len(seq) >= 23 else "", r["strand"],
            r.get("cut_site"), r.get("distance"),
            "是" if r.get("allele_covered") else "否",
            r.get("allele_region_label"), r.get("allele_pos_from_pam"),
            r.get("eff_score"), r.get("eff_score_mutant"),
            r.get("allele_delta"), r.get("allele_delta_100"),
            r.get("allele_score"), r.get("eff_score"), r.get("cfd_score"),
            r.get("dist_score"), r.get("final_v2"), r.get("grade"),
            r.get("final_score"), r.get("rank"),
        ])
    style(ws, [6, 26, 24, 6, 5, 13, 11, 9, 12, 12, 16, 16, 10, 13,
               14, 10, 12, 12, 9, 10, 10, 8])
    for row in ws.iter_rows(min_row=2):
        for c in row[12:19]:
            c.number_format = "0.00"
        if row[7].value == "是":
            for c in row:
                c.fill = hot_fill

    # Sheet2 前后对比
    ws2 = wb.create_sheet("前后对比")
    ws2.append(["sgRNA(20nt+PAM)", "v1排名", "v2排名", "排名变化",
                "v1综合分", "v2综合分", "v2-四项分",
                "ΔScore", "等位基因特异性", "覆盖突变", "突变位置"])
    big = []
    for r in rows:
        v1 = r.get("rank")
        v2 = r.get("rank_v2")
        change = (v1 - v2) if (v1 and v2) else None
        big.append((abs(change) if change is not None else -1, r, v1, v2, change))
    for _, r, v1, v2, change in sorted(big, key=lambda x: -x[0])[:40]:
        ws2.append([
            r["sequence"], v1, v2, change, r.get("final_score"), r.get("final_v2"),
            f"{fmt(r.get('allele_score'),1)}/{fmt(r.get('eff_score'),0)}/"
            f"{fmt(r.get('cfd_score'),0)}/{fmt(r.get('dist_score'),0)}",
            r.get("allele_delta"), r.get("allele_score"),
            "是" if r.get("allele_covered") else "否", r.get("allele_region_label"),
        ])
    style(ws2, [26, 8, 8, 9, 10, 10, 20, 9, 14, 9, 12])

    # Sheet3 模型说明
    ws3 = wb.create_sheet("模型说明")
    meta = merged["allele_model"]
    tm = meta["mutation"]
    ws3.append(["项目", "内容"])
    ws3.append(["综合分公式", "0.35×等位基因特异性 + 0.30×切割效率 + 0.25×脱靶特异性 + 0.10×突变距离"])
    ws3.append(["等位基因特异性", "0.6×ΔScore + 0.4×位置分；ΔScore<20 时 ×0.3"])
    ws3.append(["Δraw 定义", "Doench'16(突变型序列) − Doench'16(野生型序列)，单位是 Doench 分"])
    ws3.append(["ΔScore 归一化",
                f"100 × clamp(Δraw, 0, {AS.DELTA_FULL_SCALE}) / {AS.DELTA_FULL_SCALE}"
                f"（Doench'16 差 {AS.DELTA_FULL_SCALE} 分算满分。"
                f"实测原始差值最大只有 10 分，不归一化会让 35% 的这一维整体归零）"])
    ws3.append(["突变", f"{tm['gene']} {tm['hgvs_c']} {tm['hgvs_p']} / {tm['hgvs_g']}"])
    ws3.append(["参考碱基", f"{tm['ref']} → {tm['alt']}（codon {tm['codon_wt']} → {tm['codon_mut']}）"])
    ws3.append(["位置分", "PAM 区(21-23)=100 / Seed 区(1-8)=80 / 中间区(9-14)=40 / PAM 远端(15-20)=10 / 不覆盖=0"])
    ws3.append(["距离维度", "≤20bp=100 / ≤50=90 / ≤100=70 / ≤200=40 / >200=10；覆盖突变 +10，封顶 100"])
    ws3.append(["进入评分的 guide 数", len(rows)])
    ws3.append(["覆盖突变位点", merged.get("allele_covered_count")])
    ws3.append(["等级分布", json.dumps(merged["grade_counts"], ensure_ascii=False)])
    ws3.append(["数据来源", "CRISPOR（野生型/突变型各一次）+ CHOPCHOP，真实抓取"])
    style(ws3, [24, 100])

    wb.save(OUT / "ranking_v2.xlsx")


if __name__ == "__main__":
    raise SystemExit(main())
