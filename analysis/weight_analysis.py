#!/usr/bin/env python3
# 本项目由华南师范大学附属中学知识城校区高二一班同学制作
# Copyright (C) 2026 华南师范大学附属中学知识城校区高二一班
# SPDX-License-Identifier: AGPL-3.0-or-later
"""P0-1：权重的客观赋权 + 消融实验 + 稳定性分析。

现在四个权重（35/30/25/10）是课题组拍的经验值，答辩最容易被问：
「你凭什么是 35 不是 25？」本脚本给三样能直接写进论文的东西：

1) CRITIC 客观赋权（Diakoulaki et al., 1995）
   一种经典的多准则客观赋权法，只用数据本身算权重，不需要人为设定：
       对比强度  σ_j = 第 j 个维度归一化后的标准差
       冲突性    R_j = Σ_k (1 − r_jk)        与其他维度的相关性越低，信息越不重复
       信息量    C_j = σ_j × R_j
       权重      w_j = C_j / Σ C_j
   它和「拍脑袋权重」的差别，就是论文里「方法」那一节的内容。

2) 消融实验（ablation）
   去掉某一维、其余按比例重新归一化，看排序怎么变。
   如果去掉某一维排序几乎不动，说明这一维没在干活 —— 评委一问就能答。

3) 稳定性分析
   在基准权重附近随机扰动（Dirichlet），看 Top-1 命中率和 Top-10 重合度，
   回答「权重调一调结论会不会翻」。

⚠ 说明：CRITIC 给的是「数据驱动的一组参考权重」，**不是**最优权重。
   它的目标不是拟合实验，所以不能宣称「我们的权重更准」；
   它的价值是提供一个独立参照，证明经验权重落在合理区间内。
"""

from __future__ import annotations

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

DIMS = ["allele", "doench16", "cfd_spec", "distance"]
DIM_LABEL = {"allele": "等位基因特异性", "doench16": "切割效率",
             "cfd_spec": "脱靶特异性", "distance": "突变距离"}
BASE = {"allele": 35, "doench16": 30, "cfd_spec": 25, "distance": 10}
N_BOOT = 2000
SEED = 20261107


# ---------------------------------------------------------------------------

def load_rows():
    from genome_utils import DEFAULT_TARGET_POS
    from sgweb import allele_specific as AS
    from sgweb.aggregate import merge, rescore_v2

    cri = json.loads((DATA / "crispor_payload.json").read_text(encoding="utf-8"))
    cho = json.loads((DATA / "chopchop_payload.json").read_text(encoding="utf-8"))
    mut = json.loads((DATA / "crispor_mutant_payload.json").read_text(encoding="utf-8"))
    merged = merge(cri, cho, DEFAULT_TARGET_POS, hard_filter=True)
    di = AS.build_delta_index(cri["values"], mut["values"])
    rescore_v2(merged, di, target_pos=DEFAULT_TARGET_POS,
               genome_offset=cri.get("genome_offset"), hard_filter=True)
    return merged["rows_v2"], cri.get("genome_offset")


def dim_value(row, dim):
    if dim == "allele":
        return row.get("allele_score")
    if dim == "doench16":
        return row.get("eff_score")
    if dim == "cfd_spec":
        return row.get("cfd_score")
    if dim == "distance":
        return row.get("dist_score")
    return None


def matrix(rows, dims=DIMS):
    """返回 {dim: [值...]}，缺失填该维均值（保住样本数一致）。"""
    out = {}
    for d in dims:
        vals = [dim_value(r, d) for r in rows]
        good = [v for v in vals if v is not None]
        fill = statistics.mean(good) if good else 0.0
        out[d] = [float(fill if v is None else v) for v in vals]
    return out


def minmax(v):
    lo, hi = min(v), max(v)
    if hi == lo:
        return [0.5] * len(v)
    return [(x - lo) / (hi - lo) for x in v]


def pearson(a, b):
    n = len(a)
    ma, mb = statistics.mean(a), statistics.mean(b)
    num = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    den = math.sqrt(sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b))
    return 0.0 if den == 0 else num / den


def spearman(xs, ys):
    def ranks(v):
        n = len(v)
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
    return pearson(ranks(xs), ranks(ys))


def critic_weights(rows) -> dict:
    """CRITIC 客观赋权。"""
    m = matrix(rows)
    norm = {d: minmax(m[d]) for d in DIMS}
    sigma = {d: (statistics.stdev(norm[d]) if len(norm[d]) > 1 else 0.0) for d in DIMS}
    conflict = {}
    for j in DIMS:
        s = 0.0
        for k in DIMS:
            if j == k:
                continue
            s += 1.0 - pearson(norm[j], norm[k])
        conflict[j] = s
    info = {d: sigma[d] * conflict[d] for d in DIMS}
    tot = sum(info.values()) or 1.0
    return {d: info[d] / tot for d in DIMS}


def corr_table(rows) -> dict:
    m = matrix(rows)
    norm = {d: minmax(m[d]) for d in DIMS}
    return {j: {k: round(pearson(norm[j], norm[k]), 3) for k in DIMS} for j in DIMS}


def composite(rows, weights) -> list[float]:
    """按给定权重算综合分；weights 里没有的维度自动不参与（消融时用）。"""
    m = matrix(rows)
    tot = sum(weights.values()) or 1.0
    n = len(rows)
    return [sum(w * m[d][i] for d, w in weights.items()) / tot for i in range(n)]


def rank_of(scores):
    order = sorted(range(len(scores)), key=lambda i: -scores[i])
    r = [0] * len(scores)
    for pos, i in enumerate(order, 1):
        r[i] = pos
    return r, order


def topk_overlap(order_a, order_b, k=10) -> float:
    return len(set(order_a[:k]) & set(order_b[:k])) / k


def main() -> int:
    rows, offset = load_rows()
    print("=" * 78)
    print(f"P0-1  权重分析（{len(rows)} 条 guide，{len(DIMS)} 个维度）")
    print("=" * 78)

    base_scores = composite(rows, BASE)
    base_rank, base_order = rank_of(base_scores)

    # ---------- 1. CRITIC ----------
    print("\n1. CRITIC 客观赋权")
    print("-" * 78)
    cw = critic_weights(rows)
    ct = corr_table(rows)
    print(f"{'维度':<14}{'经验权重':>10}{'CRITIC权重':>12}{'标准差(归一)':>14}")
    m = matrix(rows)
    for d in DIMS:
        sd = statistics.stdev(minmax(m[d]))
        print(f"{DIM_LABEL[d]:<14}{BASE[d]/100*100:>9.0f}%{cw[d]*100:>11.1f}%{sd:>14.3f}")
    print("\n   维度间相关系数（归一化后）：")
    print("        " + "".join(f"{DIM_LABEL[k][:6]:>9}" for k in DIMS))
    for j in DIMS:
        print(f"   {DIM_LABEL[j][:6]:<7}" + "".join(f"{ct[j][k]:>9.2f}" for k in DIMS))
    print()

    # ---------- 2. 消融 ----------
    print("2. 消融实验（去掉一维，其余按比例重新归一化）")
    print("-" * 78)
    print(f"{'去掉的维度':<16}{'Top-1 是否变化':>16}{'Top-10 重合':>13}{'与基准排序 ρ':>14}")
    abl = {}
    for drop in DIMS:
        w = {d: BASE[d] for d in DIMS if d != drop}
        sc = composite(rows, w)
        _, order = rank_of(sc)
        same_top1 = rows[order[0]]["sequence"] == rows[base_order[0]]["sequence"]
        ov = topk_overlap(base_order, order, 10)
        rho = spearman(base_scores, sc)
        abl[drop] = {"top1_changed": not same_top1, "top10_overlap": round(ov, 3),
                     "spearman_vs_base": round(rho, 4)}
        print(f"{DIM_LABEL[drop]:<16}{'否' if same_top1 else '是':>16}"
              f"{ov:>12.0%}{rho:>14.3f}")
    print()

    # ---------- 3. 稳定性 ----------
    print("3. 稳定性：在基准权重附近随机扰动")
    print("-" * 78)
    rng = random.Random(SEED)
    scale = 40.0          # Dirichlet 浓度，越小扰动越大
    top1_hits = 0
    ovs = []
    top1_seq = rows[base_order[0]]["sequence"]
    for _ in range(N_BOOT):
        g = [rng.gammavariate(max(BASE[d], 1e-6) / scale * 10.0, 1.0) for d in DIMS]
        s = sum(g) or 1.0
        w = {d: vi / s for d, vi in zip(DIMS, g)}
        sc = composite(rows, w)
        _, order = rank_of(sc)
        if rows[order[0]]["sequence"] == top1_seq:
            top1_hits += 1
        ovs.append(topk_overlap(base_order, order, 10))
    print(f"  扰动次数 {N_BOOT}")
    print(f"  Top-1 仍为「{top1_seq[:20]}…」的比例：{top1_hits / N_BOOT:.1%}")
    print(f"  Top-10 平均重合度：{statistics.mean(ovs):.1%}"
          f"（最低 {min(ovs):.0%}，最高 {max(ovs):.0%}）")
    print()

    # ---------- 4. 结论 ----------
    print("=" * 78)
    print("结论")
    print("=" * 78)
    print(f"  · CRITIC 给出的客观权重："
          + "、".join(f"{DIM_LABEL[d]} {cw[d]*100:.0f}%" for d in DIMS))
    print(f"  · 经验权重：" + "、".join(f"{DIM_LABEL[d]} {BASE[d]}%" for d in DIMS))
    dev = max(abs(cw[d] * 100 - BASE[d]) for d in DIMS)
    print(f"  · 最大偏差 {dev:.0f} 个百分点")
    worst = max(abl, key=lambda d: abl[d]["top10_overlap"])
    print(f"  · 消融影响最小的维度：{DIM_LABEL[worst]}"
          f"（去掉后 Top-10 仍重合 {abl[worst]['top10_overlap']:.0%}）")
    best = min(abl, key=lambda d: abl[d]["top10_overlap"])
    print(f"  · 消融影响最大的维度：{DIM_LABEL[best]}"
          f"（去掉后 Top-10 只重合 {abl[best]['top10_overlap']:.0%}）")
    print()

    (OUT / "weight_analysis.json").write_text(json.dumps({
        "n_guides": len(rows),
        "base_weights": BASE,
        "critic_weights": {d: round(cw[d], 4) for d in DIMS},
        "correlation": ct,
        "ablation": abl,
        "bootstrap": {"n": N_BOOT, "top1_stability": round(top1_hits / N_BOOT, 4),
                      "top10_overlap_mean": round(statistics.mean(ovs), 4)},
    }, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"已写出：{OUT / 'weight_analysis.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
