#!/usr/bin/env python3
# 本项目由华南师范大学附属中学知识城校区高二一班同学制作
# Copyright (C) 2026 华南师范大学附属中学知识城校区高二一班
# SPDX-License-Identifier: AGPL-3.0-or-later
"""P0-2：用**独立模型**（CFD）做外部佐证。

为什么需要这一步
----------------
之前的验证都是「自己跟自己比」：
  · Doench'16 vs CHOPCHOP Efficiency  → ρ=0.9992，其实是同一个模型，不算验证
  · ΔScore 由 Doench 算出来，再用 Doench 解释 → 循环
评委一定会问「你的评分准不准」。没有湿实验时，唯一站得住的打法是：
**换一个独立的、已发表的模型，看它是否支持我们的设计。**

CFD（Cutting Frequency Determination）是 Doench 2016 用实验数据拟合出来的
**错配惩罚矩阵**：20 个位置 × 12 种错配类型，每个格子是「该位置该错配保留多少
切割活性」。它和 Doench'16（预测 on-target 效率）是两个不同的模型，
用途也不同（CFD 用于 off-target）。

本脚本用它回答两个问题
----------------------
Q1（设计层面，n=20）：我们的位置分（PAM 100 / Seed 80 / 中间 40 / 远端 10）
    是否符合 CFD 给出的「PAM 近端错配惩罚重、远端轻」的经验曲线？
Q2（个案层面，n=6）：对真正覆盖突变的那 6 条 guide，
    CFD 判定的「错配杀伤力」与我们的 ΔScore 是否同向？

数据来源
--------
mismatch_score.pkl / pam_scores.pkl —— Doench et al. 2016 补充材料，
经 github.com/MichaelLinn/off_target_prediction (CFDScoring) 公开分发。
pkl 用**受限反序列化**读取，不会执行任意代码。
"""

from __future__ import annotations

import json
import math
import pickle
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

DATA = ROOT / "data"
OUT = ROOT / "output"
OUT.mkdir(exist_ok=True)

# 我们自己的位置分：直接从 allele_specific.POSITION_TABLE 展开，
# 不手抄一份 —— 否则改了表而忘了改这里，验证会静默地验一张过期的表。
from sgweb.allele_specific import POSITION_TABLE  # noqa: E402

OUR_POSITION_SCORE = {
    pos: spec["score"]
    for spec in POSITION_TABLE.values() if spec["lo"] is not None
    for pos in range(spec["lo"], spec["hi"] + 1)
}


class _SafeUnpickler(pickle.Unpickler):
    """只允许基本类型，拒绝任何需要导入类的对象。"""

    def find_class(self, module, name):
        raise pickle.UnpicklingError(f"blocked: {module}.{name}")


def load_cfd():
    with (DATA / "mismatch_score.pkl").open("rb") as fh:
        mm = _SafeUnpickler(fh).load()
    with (DATA / "pam_scores.pkl").open("rb") as fh:
        pam = _SafeUnpickler(fh).load()
    return mm, pam


_COMP = str.maketrans("ACGTU", "TGCAA")


def revcomp_base(b: str) -> str:
    return b.translate(_COMP)


def calc_cfd(wt23: str, query23: str, mm_scores: dict, pam_scores: dict) -> float:
    """完全照 CFDScoring/cfd-score-calculator.py 的算法实现（Python3 化）。"""
    score = 1.0
    wt = wt23.upper().replace("T", "U")
    q = query23.upper().replace("T", "U")
    if len(wt) != len(q):
        raise ValueError("两条序列长度不一致")
    for i, (a, b) in enumerate(zip(wt, q)):
        if a == b:
            continue
        key = f"r{a}:d{revcomp_base(b)},{i + 1}"
        if key not in mm_scores:
            # 含 N 等非常规碱基：按最宽松处理（不惩罚），并记录
            continue
        score *= mm_scores[key]
    pam = query23.upper()[-2:]
    score *= pam_scores.get(pam, 1.0)
    return score


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


def position_curve(mm_scores: dict) -> list[dict]:
    """把 CFD 矩阵按位置汇总：每个位置「错配平均保留多少活性」。

    ⚠ 编号方向（容易搞反，这里写清楚）
        CFD 的 pos=1 是 spacer 的 **5' 端**，也就是 **PAM 远端**；
        CFD 的 pos=20 紧邻 PAM，也就是 **seed 核心**。
        即  距PAM位数 = 21 − pos。

    这个方向有两条独立佐证：
      1) 参考实现里 `sg = off[:-3]`，是 protospacer 5'→3' 顺序；
      2) 算出来的惩罚曲线在 pos 15–20 最重（保留活性最低），
         正好对应「PAM 近端错配最致命」的文献共识 —— 方向反了会得到相反的结论。

    PAM 本身（第 21–23 位）不在错配矩阵里，由单独的 pam_scores 处理，这里不计。
    """
    per: dict[int, list[float]] = {}
    for key, val in mm_scores.items():
        try:
            pos = int(key.rsplit(",", 1)[1])
        except (IndexError, ValueError):
            continue
        per.setdefault(pos, []).append(val)
    rows = []
    for pos in sorted(per):
        vals = per[pos]
        dist_from_pam = 21 - pos          # CFD pos=1 -> 距PAM 20；pos=20 -> 距PAM 1
        rows.append({
            "cfd_position": pos,
            "pos_from_pam": dist_from_pam,
            "mean_retained": round(statistics.mean(vals), 4),
            "mean_penalty": round(1.0 - statistics.mean(vals), 4),
            "n_types": len(vals),
            "ours": OUR_POSITION_SCORE.get(dist_from_pam),
        })
    return rows


def main() -> int:
    from genome_utils import normalize_strand
    from sgweb import allele_specific as AS
    from sgweb.aggregate import merge, rescore_v2

    mm, pam = load_cfd()
    print("=" * 76)
    print("外部佐证：用 CFD（Doench 2016 错配矩阵）检验我们的设计")
    print("=" * 76)
    print(f"CFD 错配矩阵：{len(mm)} 个 (错配类型, 位置)；PAM 矩阵：{len(pam)} 项")
    print()

    # ---------- Q1：位置曲线 ----------
    curve = position_curve(mm)
    spacer = sorted(curve, key=lambda r: r["pos_from_pam"])   # 按「距PAM」排序
    print("Q1  位置分设计是否与 CFD 的经验曲线同向")
    print("     （下表按「距 PAM 的位数」排序，1 = 紧邻 PAM）")
    print("-" * 76)
    print(f"{'距PAM':>5}{'CFD位置':>8}{'错配后平均保留活性':>20}{'平均惩罚':>12}{'我们的位置分':>14}")
    for r in spacer:
        print(f"{r['pos_from_pam']:>5}{r['cfd_position']:>8}{r['mean_retained']:>20.3f}"
              f"{r['mean_penalty']:>12.3f}{r['ours']:>14}")
    print()

    xs = [r["mean_penalty"] for r in spacer]
    ys = [float(r["ours"]) for r in spacer]
    rho_spacer = spearman(xs, ys)
    print(f"spacer 20 个位置：CFD 惩罚 vs 我们的位置分   Spearman ρ = {rho_spacer:.3f}")

    seed_pen = statistics.mean(r["mean_penalty"] for r in spacer if r["pos_from_pam"] <= 8)
    mid_pen = statistics.mean(r["mean_penalty"] for r in spacer if 9 <= r["pos_from_pam"] <= 14)
    far_pen = statistics.mean(r["mean_penalty"] for r in spacer if r["pos_from_pam"] >= 15)
    print(f"  分区平均惩罚：Seed(距PAM 1-8) {seed_pen:.3f} | 中间 {mid_pen:.3f} | "
          f"PAM远端 {far_pen:.3f}")
    print(f"  → Seed 惩罚 / 远端惩罚 = {seed_pen / far_pen:.2f} 倍")
    print()

    # ---------- Q2：6 条覆盖突变的 guide ----------
    cri = json.loads((DATA / "crispor_payload.json").read_text(encoding="utf-8"))
    cho = json.loads((DATA / "chopchop_payload.json").read_text(encoding="utf-8"))
    mutp = json.loads((DATA / "crispor_mutant_payload.json").read_text(encoding="utf-8"))

    wt_seq = "".join((ROOT / "SFTPC_测试序列_hg38_2213bp.txt")
                     .read_text(encoding="utf-8").split()).upper()
    mut_seq = "".join((DATA / "SFTPC_I73T_mutant_2213bp.txt")
                      .read_text(encoding="utf-8").split()).upper()
    offset = cri.get("genome_offset")
    seq_target = AS.MUTATION["pos"] - int(offset)

    merged = merge(cri, cho, AS.MUTATION["pos"], hard_filter=True)
    di = AS.build_delta_index(cri["values"], mutp["values"])
    rescore_v2(merged, di, target_pos=AS.MUTATION["pos"],
               genome_offset=offset, hard_filter=True)

    print("Q2  覆盖突变的 6 条 guide：CFD 判定的错配杀伤力 vs 我们的 ΔScore")
    print("-" * 76)

    def rc23(s):
        return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]

    rows = []
    for r in merged["rows_v2"]:
        if not r.get("allele_covered"):
            continue
        lo, hi = AS.guide_window(r["seq_pam_start"], r["strand"])
        w23 = wt_seq[lo - 1:hi]
        m23 = mut_seq[lo - 1:hi]
        if normalize_strand(r["strand"]) == "-":
            w23, m23 = rc23(w23), rc23(m23)
        cfd_wt = calc_cfd(w23, w23, mm, pam)
        cfd_mut = calc_cfd(w23, m23, mm, pam)
        ratio = (cfd_wt / cfd_mut) if cfd_mut > 0 else float("inf")
        rows.append({
            "sequence": r["sequence"],
            "index_in_23mer": r["allele_pos_from_pam"],
            "region": r["allele_region_label"],
            "cfd_wt": round(cfd_wt, 4),
            "cfd_mut": round(cfd_mut, 4),
            "cfd_discrimination": round(ratio, 3),
            "our_delta": r["allele_delta"],
            "our_delta_100": r["allele_delta_100"],
            "ass": r["allele_score"],
            "final_v2": r["final_v2"],
        })

    print(f"{'sgRNA':<24}{'位置':<10}{'CFD(WT)':>9}{'CFD(Mut)':>9}{'CFD区分度':>10}"
          f"{'Δraw':>7}{'ASS':>7}")
    for r in sorted(rows, key=lambda x: -(x["cfd_discrimination"] or 0)):
        print(f"{r['sequence']:<24}{r['region']:<10}{r['cfd_wt']:>9.3f}"
              f"{r['cfd_mut']:>9.3f}{r['cfd_discrimination']:>10.2f}"
              f"{r['our_delta']:>7.1f}{r['ass']:>7.1f}")
    print()

    xs2 = [r["cfd_discrimination"] for r in rows]
    ys2 = [abs(r["our_delta"] or 0) for r in rows]
    rho2 = spearman(xs2, ys2)
    print(f"CFD 区分度 vs |ΔScore|  Spearman ρ = "
          f"{rho2:.3f}  (n={len(rows)}，样本很小，只能看方向)" if rho2 is not None
          else "样本太少，无法算相关")
    n_pos = sum(1 for r in rows if (r["our_delta"] or 0) > 0)
    print(f"Δraw > 0（偏好突变型）的 guide：{n_pos} / {len(rows)}")
    print()

    # ---------- 结论 ----------
    print("=" * 76)
    print("结论（可直接写进论文）")
    print("=" * 76)
    verdict = []
    if rho_spacer is not None and rho_spacer > 0.5:
        verdict.append(
            f"① 位置分设计得到独立模型支持：CFD 错配惩罚与我们的位置分"
            f"Spearman ρ={rho_spacer:.3f}，Seed 区惩罚是 PAM 远端的 "
            f"{seed_pen / far_pen:.1f} 倍——与「PAM 近端错配最关键」的文献共识一致。")
    else:
        verdict.append(
            f"① 位置分与 CFD 曲线的一致性只有 ρ={rho_spacer}，需要重新审视分档。")
    verdict.append(
        f"② 对 6 条覆盖突变的 guide，CFD 判定的区分度与 |ΔScore| 的方向"
        f"一致性 ρ={rho2:.2f}（n=6，不具统计效力，只作方向性参考）。"
        if rho2 is not None else
        "② 覆盖突变的 guide 仅 6 条，无法做有统计效力的相关分析。")
    verdict.append(
        "③ 局限：目前**没有**公开的、带实测编辑效率的等位基因特异性 sgRNA 数据集，"
        "所以本工具无法做真正的回溯验证；CFD 只能证明「设计原则与文献一致」，"
        "不能证明「我们的排序与实验结果一致」。要补这一环，必须做体外切割或细胞实验。")
    for v in verdict:
        print("  " + v)
    print()

    (OUT / "cfd_validation.json").write_text(json.dumps({
        "source": "Doench et al. 2016 CFD matrix (mismatch_score.pkl / pam_scores.pkl)",
        "n_mismatch_entries": len(mm),
        "position_curve": curve,
        "q1_spearman_spacer": rho_spacer,
        "seed_vs_far_penalty_ratio": round(seed_pen / far_pen, 3),
        "q2_guides": rows,
        "q2_spearman": rho2,
        "verdict": verdict,
    }, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"已写出：{OUT / 'cfd_validation.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
