#!/usr/bin/env python3
# 本项目由华南师范大学附属中学知识城校区高二一班同学制作
# Copyright (C) 2026 华南师范大学附属中学知识城校区高二一班
# SPDX-License-Identifier: AGPL-3.0-or-later
"""体外切割实验的计算机模拟 + 跨模型一致性检验。

⚠ 重要区分（写论文时也要写清楚）：
  A 部分是**真实数据**：CRISPOR 的 Doench'16 与 CHOPCHOP 的 Efficiency 是两个
    互不依赖的模型，对同一批 guide 的排序一致性是实打实的证据。
  B 部分是**模拟**：在给定假设下"如果真去做体外切割，胶会长什么样、读出大约多少"。
    它是实验设计和预期结果的推演，**不是**验证，不能当作实验数据引用。

输出：
  output/sim_gel.png          模拟凝胶图
  output/sim_results.json     每条 guide 的预期切割率 / 等位基因特异性指数 / 95%CI
  output/sim_table.csv        可直接贴进论文的表
"""

from __future__ import annotations

import csv
import json
import math
import random
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

DATA = ROOT / "data"
OUT = ROOT / "output"
OUT.mkdir(exist_ok=True)

TARGET = 22_163_096
AMP_START, AMP_END = 22_162_900, 22_163_300     # 401 bp PCR 扩增子
N_REP = 3
GEL_SD = 0.035                                   # 凝胶光密度定量的典型误差
SEED = 20261107


# ---------------------------------------------------------------------------
# A. 跨模型一致性：Doench'16 vs CHOPCHOP Efficiency（真实数据）
# ---------------------------------------------------------------------------

def load_chopchop_efficiency() -> dict:
    """从 CHOPCHOP 结果表里取 Efficiency 列，键 = (location, 归一化链向)。

    结果表可能在两个位置：
      data/chopchop_results.tsv    —— 随源码一起分发的那份（交付包里就靠它）
      .chopchop_cache/results_*.tsv —— 本地跑完取数后自动生成的缓存
    两处都没有时返回空字典，调用方降级处理，不抛异常。
    """
    from genome_utils import normalize_strand
    # 优先只用随源码分发的那份，避免缓存里其它案例（如 HBB）混进来串味
    fixed = DATA / "chopchop_results.tsv"
    candidates = [fixed] if fixed.is_file() \
        else sorted((ROOT / ".chopchop_cache").glob("results_*.tsv"))

    out = {}
    for p in candidates:
        with p.open(encoding="utf-8", newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                loc = str(row.get("Genomic location") or "")
                st = normalize_strand(row.get("Strand"))
                eff = row.get("Efficiency")
                if not loc or st is None or eff in (None, "", "N/A"):
                    continue
                loc = loc.split(":")[-1]
                try:
                    out[(int(loc), st)] = float(eff)
                except ValueError:
                    continue
    return out


def spearman(xs, ys):
    n = len(xs)
    if n < 3:
        return None

    def ranks(v):
        order = sorted(range(n), key=lambda i: v[i])
        r = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r

    rx, ry = ranks(xs), ranks(ys)
    mx, my = statistics.mean(rx), statistics.mean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return None if den == 0 else num / den


def cross_model_check() -> dict:
    from genome_utils import normalize_strand
    cri = json.loads((DATA / "crispor_payload.json").read_text(encoding="utf-8"))
    cho = json.loads((DATA / "chopchop_payload.json").read_text(encoding="utf-8"))
    eff = load_chopchop_efficiency()

    pairs = []
    for v in cho["values"]:
        loc = v.get("location")
        st = normalize_strand(v.get("strand"))
        if loc is None or st is None:
            continue
        try:
            loc = int(str(loc).split(":")[-1])
        except ValueError:
            continue
        if (loc, st) in eff:
            pairs.append((v.get("seq_pam_start"), st, eff[(loc, st)]))

    by_key = {(v.get("seq_pam_start"), normalize_strand(v.get("strand"))): v.get("doench16")
              for v in cri["values"]}
    xs, ys = [], []
    for sps, st, e in pairs:
        d = by_key.get((sps, st))
        if d is None:
            continue
        xs.append(float(d))
        ys.append(float(e))

    rho = spearman(xs, ys)
    if not xs:
        return {"n": 0, "spearman": None, "pearson": None,
                "doench_mean": None, "eff_mean": None,
                "note": "缺少 CHOPCHOP 结果表（data/chopchop_results.tsv），跳过跨模型检验"}
    mx, my = statistics.mean(xs), statistics.mean(ys)
    num = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    den = math.sqrt(sum((a - mx) ** 2 for a in xs) * sum((b - my) ** 2 for b in ys))
    pear = num / den if den else None
    return {
        "n": len(xs),
        "spearman": round(rho, 4) if rho is not None else None,
        "pearson": round(pear, 4) if pear is not None else None,
        "doench_mean": round(mx, 2), "eff_mean": round(my, 2),
    }


# ---------------------------------------------------------------------------
# B. 体外切割模拟
# ---------------------------------------------------------------------------

def doench_to_cleavage(d: float | None) -> float:
    """Doench'16 (0-100) -> 预期切割比例。

    单条 logistic 映射，锚点取「Doench 50 约等于一半底物被切开」，
    上下限 3% / 92% 对应真实体外反应的背景与不完全反应。
    """
    if d is None:
        return 0.0
    lo, hi = 0.03, 0.92
    return lo + (hi - lo) / (1.0 + math.exp(-(float(d) - 50.0) / 18.0))


def simulate_note(row: dict, rng: random.Random) -> dict:
    """一条 guide 在野生型/突变型底物上的切割模拟。"""
    cut = row.get("cut_site")
    cut = None if cut is None else float(cut)
    if cut is None or not (AMP_START < cut < AMP_END):
        return {}
    f1 = round(cut - AMP_START)
    f2 = round(AMP_END - cut)

    p_wt = doench_to_cleavage(row.get("eff_score"))
    p_mut = doench_to_cleavage(row.get("eff_score_mutant")) \
        if row.get("eff_score_mutant") is not None else p_wt

    def measure(p, rng):
        # 凝胶定量：切割条带之和 / 总信号，加高斯噪声
        v = p + rng.gauss(0.0, GEL_SD)
        return max(0.0, min(1.0, v))

    wt = [measure(p_wt, rng) for _ in range(N_REP)]
    mu = [measure(p_mut, rng) for _ in range(N_REP)]
    m_wt, m_mut = statistics.mean(wt), statistics.mean(mu)
    idx = (m_mut / m_wt) if m_wt > 0.01 else None

    # 指数的不确定度（比值传播）
    se = None
    if m_wt > 0.05 and m_mut > 0.05:
        se_wt = statistics.stdev(wt) / math.sqrt(N_REP) if N_REP > 1 else 0
        se_mut = statistics.stdev(mu) / math.sqrt(N_REP) if N_REP > 1 else 0
        rel = math.sqrt((se_wt / m_wt) ** 2 + (se_mut / m_mut) ** 2)
        se = idx * rel
    return {
        "sequence": row.get("sequence"),
        "strand": row.get("strand"),
        "cut_site": int(cut),
        "frag1": f1, "frag2": f2,
        "doench_wt": row.get("eff_score"),
        "doench_mut": row.get("eff_score_mutant"),
        "pred_cleavage_wt": round(p_wt, 4),
        "pred_cleavage_mut": round(p_mut, 4),
        "pred_index": round(p_mut / p_wt, 3) if p_wt > 0.01 else None,
        "sim_cleavage_wt": round(m_wt, 4),
        "sim_cleavage_mut": round(m_mut, 4),
        "allele_index": None if idx is None else round(idx, 3),
        "allele_index_se": None if se is None else round(se, 3),
        "delta_score": row.get("allele_delta_100"),
        "region": row.get("allele_region_label"),
        "ass": row.get("allele_score"),
    }


def build_gel(path: Path, results: list[dict]) -> None:
    """画一张模拟凝胶：Ladder / 无 Cas9 对照 / 每条 guide 的 WT 与 Mut 泳道。"""
    from PIL import Image, ImageDraw, ImageFont

    lanes = []
    lanes.append(("Ladder", [(400, 0.85), (300, 0.8), (200, 0.75), (100, 0.7)]))
    lanes.append(("无 Cas9\nWT", [(401, 0.95)]))
    lanes.append(("无 Cas9\nMut", [(401, 0.95)]))
    for r in results:
        lanes.append((f"{r['sequence'][:8]}…\nWT",
                      [(401, 1.0 - r["sim_cleavage_wt"])]
                      + ([(r["frag1"], r["sim_cleavage_wt"] / 2),
                          (r["frag2"], r["sim_cleavage_wt"] / 2)])))
        lanes.append((f"{r['sequence'][:8]}…\nMut",
                      [(401, 1.0 - r["sim_cleavage_mut"])]
                      + ([(r["frag1"], r["sim_cleavage_mut"] / 2),
                          (r["frag2"], r["sim_cleavage_mut"] / 2)])))

    W = 120 + len(lanes) * 92
    H = 560
    img = Image.new("RGB", (W, H), "#1a1a1a")
    d = ImageDraw.Draw(img)
    try:
        font = ImageFont.truetype("arial.ttf", 14)
        small = ImageFont.truetype("arial.ttf", 12)
    except Exception:
        font = small = ImageFont.load_default()

    # 中文用系统里带 CJK 字形的那一档（arial 没有汉字，会画成方框）
    cjk = None
    for cand in ("C:/Windows/Fonts/msyh.ttc", "C:/Windows/Fonts/msyhbd.ttc",
                 "C:/Windows/Fonts/simhei.ttf", "C:/Windows/Fonts/simsun.ttc",
                 "/System/Library/Fonts/PingFang.ttc",
                 "/usr/share/fonts/truetype/wqy/wqy-zenhei.ttc"):
        try:
            cjk = ImageFont.truetype(cand, 14)
            cjk_small = ImageFont.truetype(cand, 12)
            break
        except Exception:
            continue
    if cjk is None:
        cjk, cjk_small = font, small

    # 泳道底
    for i, (label, bands) in enumerate(lanes):
        x0 = 90 + i * 92
        d.rectangle([x0, 60, x0 + 64, 460], fill="#0d0d0d", outline="#333")
        d.text((x0, 470), label, fill="#ddd", font=cjk_small, align="left")

    # 100 bp 梯子刻度
    def ypos(bp):
        return 60 + (1 - math.log10(max(bp, 20)) / math.log10(450)) * 400

    for bp in (400, 300, 200, 100):
        y = ypos(bp)
        d.line([20, y, 80, y], fill="#666", width=1)
        d.text((16, y - 7), str(bp), fill="#999", font=cjk_small)

    # 条带
    for i, (label, bands) in enumerate(lanes):
        x0 = 90 + i * 92
        for bp, inten in bands:
            y = ypos(bp)
            shade = int(40 + 215 * max(0.0, min(1.0, inten)))
            col = (shade, shade, shade)
            d.rectangle([x0 + 4, y - 3, x0 + 60, y + 3], fill=col)

    d.text((18, 22), "模拟凝胶（8% PAGE，示意，非实验数据）", fill="#ffd479", font=cjk)
    d.text((18, 520), "401 bp 扩增子 / SpCas9 RNP 体外切割 / 37°C 1 h",
           fill="#777", font=cjk_small)
    img.save(path)


def main() -> int:
    from sgweb import allele_specific as AS
    from sgweb.aggregate import merge, rescore_v2

    cri = json.loads((DATA / "crispor_payload.json").read_text(encoding="utf-8"))
    cho = json.loads((DATA / "chopchop_payload.json").read_text(encoding="utf-8"))
    mut = json.loads((DATA / "crispor_mutant_payload.json").read_text(encoding="utf-8"))

    print("=" * 74)
    print("A. 跨模型一致性（真实数据，不是模拟）")
    print("=" * 74)
    cm = cross_model_check()
    print(f"CRISPOR Doench'16 vs CHOPCHOP Efficiency")
    print(f"  配对 guide 数 : {cm['n']}")
    if cm["n"] == 0:
        print(f"  ⚠ {cm.get('note')}")
        print()
    else:
        print(f"  Spearman ρ    : {cm['spearman']}")
    print(f"  Pearson   r   : {cm['pearson']}")
    print(f"  均值          : Doench {cm['doench_mean']}  /  Efficiency {cm['eff_mean']}")
    print()

    merged = merge(cri, cho, TARGET, hard_filter=True)
    di = AS.build_delta_index(cri["values"], mut["values"])
    rescore_v2(merged, di, target_pos=TARGET,
               genome_offset=cri.get("genome_offset"), hard_filter=True)
    rows = merged["rows_v2"]

    covered = [r for r in rows if r.get("allele_covered")]
    # 阴性对照：切点落在同一个扩增子里、但 23mer 盖不到突变位点的那几条
    # —— 它们一定同时切两条等位基因，用来验证「胶上应该看到 1:1」
    controls = [r for r in rows
                if not r.get("allele_covered")
                and r.get("cut_site") is not None
                and AMP_START < float(r["cut_site"]) < AMP_END][:2]
    picked = sorted(covered, key=lambda r: -(r.get("allele_score") or 0))[:4] + controls

    rng = random.Random(SEED)
    results = [x for x in (simulate_note(r, rng) for r in picked) if x]

    print("=" * 74)
    print("B. 体外切割模拟（示意，非实验数据）")
    print("=" * 74)
    print(f"扩增子 hg38 chr8:{AMP_START}-{AMP_END}（{AMP_END - AMP_START} bp）")
    print()
    hdr = (f"{'sgRNA':<22}{'切点':>10}{'片段':>10}{'WT切割':>8}{'Mut切割':>8}"
           f"{'预期指数':>9}{'模拟指数':>9}{'±SE':>7}{'ΔScore':>8}")
    print(hdr)
    print("-" * len(hdr))
    for r in results:
        idx = "-" if r["allele_index"] is None else f"{r['allele_index']:.2f}"
        se = "-" if r["allele_index_se"] is None else f"{r['allele_index_se']:.2f}"
        pi = "-" if r.get("pred_index") is None else f"{r['pred_index']:.2f}"
        print(f"{r['sequence']:<22}{r['cut_site']:>10}"
              f"{str(r['frag1']) + '/' + str(r['frag2']):>10}"
              f"{r['sim_cleavage_wt']:>8.2f}{r['sim_cleavage_mut']:>8.2f}"
              f"{pi:>9}{idx:>9}{se:>7}{str(r['delta_score']):>8}")
    print()

    # 上样量相同的话，等位基因特异性指数 >1 表示偏爱突变型
    best = max((r for r in results if r["allele_index"] is not None),
               key=lambda r: r["allele_index"], default=None)
    if best:
        print(f"预期最好的候选：{best['sequence']}")
        print(f"  野生型切割 {best['sim_cleavage_wt']:.2f}，突变型 {best['sim_cleavage_mut']:.2f}，"
              f"特异性指数 {best['allele_index']:.2f} ± {best['allele_index_se']:.2f}")
        print(f"  注意：指数 1.0 = 两条切得一样。模拟里 {GEL_SD:.3f} 的灰度误差 → "
              f"3 次重复的标准误约 {best['allele_index_se']:.2f}，"
              f"要区分 1.0 与 1.2 需要更多重复或 qPCR 定量。")
        print()

    build_gel(OUT / "sim_gel.png", results)
    (OUT / "sim_results.json").write_text(
        json.dumps({"cross_model": cm, "amplicon": [AMP_START, AMP_END],
                    "gel_sd": GEL_SD, "replicates": N_REP, "results": results},
                   ensure_ascii=False, indent=2), encoding="utf-8")
    with (OUT / "sim_table.csv").open("w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(results[0].keys()))
        w.writeheader()
        w.writerows(results)

    print(f"已写出：{OUT / 'sim_gel.png'}")
    print(f"已写出：{OUT / 'sim_results.json'}")
    print(f"已写出：{OUT / 'sim_table.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
