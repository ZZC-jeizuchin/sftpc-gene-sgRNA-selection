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
"""三个数据合并 + 加权排序。

权重（四个数据的定义.txt）：
    7  Doench '16   均匀映射到 0-1（高考赋分式：按名次均匀铺开）
    6  CFD Spec     除以 100
    4  距离          分档打分（<20bp=1.0, 20-50=0.9, 50-100=0.7, 100-200=0.4, >200=0.1）
                     + 覆盖奖励 0.1，上限 1.0

    final = (7*d16' + 6*cfd/100 + 4*dist) / 17

=====================================================================
负责人定的三部分规则：

    按 20nt spacer / PAM 起点对齐后分三堆：
        第一部分 = 两个网站都有            → 三项齐全
        第二部分 = cris(CRISPOR) 有、跳跳(CHOPCHOP) 没有
        第三部分 = 跳跳有、cris 没有        → 【直接扔掉】

    第二部分的第三项：它没有跳跳的切割位点数据。
        → 能算就算（CRISPOR 其实也给切点，实测和 CHOPCHOP 完全一致）
        → 算不出来才用【第一部分第三项得分的平均值】顶替

    极端情况：第一部分为空（两站完全不重合）
        → 第三项 = 滚木 = null
=====================================================================
"""

from __future__ import annotations

from genome_utils import (
    DEFAULT_TARGET_POS, bands_json, distance_band_score, normalize_strand,
)

WEIGHTS = {"doench16": 7, "cfd_spec": 6, "distance": 4}
TOTAL_WEIGHT = sum(WEIGHTS.values())          # 17

# 覆盖奖励（负责人 2026-10-06 最终定版）：
#   guide 的 20nt 区间压住目标点 (X <= 目标 <= X+20) → 第三项 +0.1
#   没压住                                        → 第三项 −0.1
#   第三项上限 1.0，下限 0.0
#
# 【为什么必须要有那个 −0.1】：
#   切点永远在 [X, X+20] 这个窗口内部（正链偏移 17，反链偏移 6），
#   所以"压住目标"必然意味着切点距离 < 20bp → 分档本来就是 1.0。
#   如果只加不减，1.0+0.1 会被上限 1.0 吃掉，加分永远是废的。
#   加上 −0.1 之后：压住的 = 1.0，没压住的 = 0.9，两者才真正分开。
COVERAGE_BONUS = 0.1      # 压住目标点 → +
COVERAGE_PENALTY = 0.1    # 没压住     → −
THIRD_ITEM_CAP = 1.0
# 下限取 0.1（= 分档表最低那一档），不是 0。
# 原因：扣 0.1 不能把最差档（>200bp = 0.1）压到 0 —— 那样"超出范围"的会变成 0 分，
# 而负责人要求它们保持 0.1。所以扣分只对 >=0.4 的那几档实际起作用。
THIRD_ITEM_FLOOR = 0.1


def uniform_map(values: list[float | None], higher_is_better: bool = True) -> list[float | None]:
    """名次均匀映射到 [0,1]（高考赋分式）。

    **最好的 = 1.0，最差的 = 0.0**，中间按名次等距铺开；并列取平均名次。
    higher_is_better=True 时"最大的最好"（Doench'16 就是这种）。
    """
    idx = [i for i, v in enumerate(values) if v is not None]
    n = len(idx)
    out: list[float | None] = [None] * len(values)
    if n == 0:
        return out
    if n == 1:
        out[idx[0]] = 1.0
        return out

    order = sorted(idx, key=lambda i: values[i], reverse=higher_is_better)
    pos = 0
    while pos < n:
        end = pos
        while end + 1 < n and values[order[end + 1]] == values[order[pos]]:
            end += 1
        avg_rank = (pos + end) / 2.0
        mapped = round(1.0 - avg_rank / (n - 1), 6)
        for k in range(pos, end + 1):
            out[order[k]] = mapped
        pos = end + 1
    return out


def spacer_of(target_seq: str | None) -> str | None:
    """23nt -> 前 20nt（spacer）。两个网站都用 20nt+PAM 的写法，可直接对齐。"""
    if not target_seq:
        return None
    s = str(target_seq).strip().upper()
    return s[:20] if len(s) >= 20 else None


# 覆盖奖励用哪个窗口来判断「压住了目标点」：
#   "guide_start"（默认，负责人 2026-10-06 指定）
#         X = sgRNA 起点（= CHOPCHOP 网页上 chr8 后面那个数字）
#         条件：X <= 目标点 <= X+20
#   "spacer"（原来的做法）
#         用 20nt spacer 区间；正链 [pam-20, pam-1]，反链 [pam+3, pam+22]
COVERAGE_WINDOW = "guide_start"


def covers_target(pam_start, strand, target_pos: int,
                  guide_start=None, mode: str = COVERAGE_WINDOW) -> bool:
    """guide 的 20nt 区间有没有压住目标点。

    mode="guide_start"：X <= target <= X+20，X = sgRNA 起点
    mode="spacer"     ：20nt spacer 区间（正链 [pam-20,pam-1]，反链 [pam+3,pam+22]）
    """
    if mode == "guide_start":
        if guide_start is None:
            return False
        return guide_start <= target_pos <= guide_start + 20

    if pam_start is None:
        return False
    s = normalize_strand(strand)
    if s == "+":
        lo, hi = pam_start - 20, pam_start - 1
    elif s == "-":
        lo, hi = pam_start + 3, pam_start + 22
    else:
        return False
    return lo <= target_pos <= hi


def third_item_value(distance, pam_start, strand, target_pos,
                     guide_start=None, bonus: float = COVERAGE_BONUS,
                     penalty: float = COVERAGE_PENALTY,
                     window_mode: str = COVERAGE_WINDOW):
    """第三项 = 距离分档 + 覆盖奖励。

    压住目标点 → +bonus；没压住 → −penalty。
    最后夹在 [0, 1] 之间（负责人要求「不超过 1」）。
    返回 (最终分, 分档分, 加减分)。
    """
    base = distance_band_score(distance)
    if base is None:
        return None, None, None
    covered = covers_target(pam_start, strand, target_pos,
                            guide_start=guide_start, mode=window_mode)
    delta = bonus if covered else -penalty
    val = max(THIRD_ITEM_FLOOR, min(THIRD_ITEM_CAP, base + delta))
    return round(val, 6), base, round(delta, 6)


def merge(crispor_payload: dict | None, chopchop_payload: dict | None,
          target_pos: int = DEFAULT_TARGET_POS,
          third_mode: str = "auto") -> dict:
    """合并 + 分三部分 + 加权排序。

    third_mode:
        "auto"     能真算就真算，算不出来用第一部分平均值（推荐）
        "average"  完全按负责人原话：第二部分一律用第一部分平均值
        "computed" 第二部分只用真算值（算不出就 null）
    """
    # ---------------------------------------------------------- 1. 建索引
    def index(payload):
        out = {}
        if not payload:
            return out
        for v in payload["values"]:
            key = v.get("pam_start")
            if key is None:
                key = v.get("seq_pam_start")
            if key is None:
                continue
            out[key] = dict(v, _key=key)
        return out

    cri = index(crispor_payload)
    cho = index(chopchop_payload)

    # ---------------------------------------------------------- 2. 分三部分
    part1_keys = sorted(set(cri) & set(cho), key=str)      # 两站都有
    part2_keys = sorted(set(cri) - set(cho), key=str)      # cris 有、跳跳没有
    part3_keys = sorted(set(cho) - set(cri), key=str)      # 跳跳有、cris 没有 → 扔掉

    # ---------------------------------------------------------- 3. 组行
    def make_row(key, kind):
        c = cri.get(key, {})
        h = cho.get(key, {})
        src = c or h
        p = c.get("pam_start") if c.get("pam_start") is not None else h.get("pam_start")
        seqp = c.get("seq_pam_start") if c.get("seq_pam_start") is not None else h.get("seq_pam_start")
        st = normalize_strand(src.get("strand"))
        # sgRNA 起点 = 23mer 在正链上的最左端（= CHOPCHOP 网页上显示的 Genomic location）
        #   + 链：PAM 在右端，spacer 起点 = PAM起点 − 20
        #   − 链：PAM 在左端，起点就是 PAM 起点
        def gstart(pam):
            if pam is None:
                return None
            return pam - 20 if st == "+" else (pam if st == "-" else None)
        return {
            "key": key,
            "part": kind,                       # 1 / 2
            "pam_start": p,
            "seq_pam_start": seqp,
            "guide_start": gstart(p),
            "seq_guide_start": seqp - 20 if (seqp is not None and st == "+") else seqp,
            "strand": src.get("strand"),
            "sequence": src.get("target_seq"),
            "spacer": spacer_of(src.get("target_seq")),
            "cut_site": src.get("cut_site"),
            "distance": src.get("distance"),
            "doench16": c.get("doench16"),
            "cfd_spec": c.get("cfd_spec"),
            "mit_spec": c.get("mit_spec"),
            "guide_id_crispor": c.get("guide_id"),
            "guide_id_chopchop": h.get("guide_id"),
        }

    rows1 = [make_row(k, 1) for k in part1_keys]
    rows2 = [make_row(k, 2) for k in part2_keys]
    all_rows = rows1 + rows2

    # ---------------------------------------------------------- 4. 第三项
    # 第一部分：分档 + 覆盖奖励
    for r in rows1:
        v, base, extra = third_item_value(r["distance"], r["pam_start"], r["strand"], target_pos,
                                              guide_start=r["guide_start"])
        r["third_base"], r["third_bonus"] = base, extra
        r["third_item_raw"] = v
        r["third_item"] = v
        r["third_computed"] = v
        r["third_source"] = "computed" if v is not None else None

    have = [r["third_item_raw"] for r in rows1 if r["third_item_raw"] is not None]
    avg_part1 = round(sum(have) / len(have), 6) if have else None

    # 第二部分：能真算就真算，否则用第一部分平均值（按 third_mode）
    for r in rows2:
        v, base, extra = third_item_value(r["distance"], r["pam_start"], r["strand"], target_pos,
                                              guide_start=r["guide_start"])
        r["third_base"], r["third_bonus"] = base, extra
        r["third_item_raw"] = v
        r["third_computed"] = v
        if third_mode == "average":
            use, src = ((avg_part1, "avg_part1") if avg_part1 is not None else (None, None))
        elif v is not None:
            use, src = v, "computed"
        else:
            use, src = ((avg_part1, "avg_part1") if avg_part1 is not None else (None, None))
        r["third_item"] = use
        r["third_source"] = src

    # ---------------------------------------------------------- 5. 滚木
    #   第一部分为空 = 两个网站完全不重合 → 第三项全部 = null
    kurumi = (len(rows1) == 0)
    if kurumi:
        for r in all_rows:
            r["third_item"] = None
            r["third_source"] = None
        note = ("两个网站没有找到任何共同的 sgRNA（「两站都有」那批为空）→ "
                "突变距离评分 = 滚木(null)。请确认输入的这条基因是否过于变异，"
                "或者根本不是目标基因。")
    else:
        note = None

    # ---------------------------------------------------------- 6. 其他两项
    d16 = uniform_map([r["doench16"] for r in all_rows], higher_is_better=True)
    cfd_norm = [(None if r["cfd_spec"] is None else round(r["cfd_spec"] / 100.0, 6))
                for r in all_rows]

    # ---------------------------------------------------------- 7. 加权
    out = []
    for i, r in enumerate(all_rows):
        parts, missing, used = 0.0, [], 0
        vals = {"doench16": d16[i], "cfd_spec": cfd_norm[i], "distance": r["third_item"]}
        for name in ("doench16", "cfd_spec", "distance"):
            v = vals[name]
            if v is None:
                missing.append(name)
                continue
            parts += WEIGHTS[name] * v
            used += WEIGHTS[name]
        r.update({
            "doench16_mapped": d16[i],
            "cfd_spec_normalized": cfd_norm[i],
            "final_score": round(parts / used, 6) if used else None,
            "missing": missing,
            "weight_used": used,
        })
        out.append(r)

    scored = [r for r in out if r["final_score"] is not None]
    scored.sort(key=lambda r: (-r["final_score"],
                               r["distance"] if r["distance"] is not None else 1e9))
    for n, r in enumerate(scored, 1):
        r["rank"] = n

    return {
        "ok": True,
        "weights": WEIGHTS,
        "total_weight": TOTAL_WEIGHT,
        "target_pos": target_pos,
        "third_mode": third_mode,
        "coverage_bonus": COVERAGE_BONUS,
        "bands": bands_json(),
        "parts": {
            "part1": len(rows1),
            "part2": len(rows2),
            "part3_dropped": len(part3_keys),
        },
        "avg_part1_third": avg_part1,
        "kurumi": kurumi,
        "kurumi_note": note,
        "distance_available": any(r["third_item"] is not None for r in out),
        "count": len(scored),
        "rows": scored,
        "unscored": [r for r in out if r["final_score"] is None],
        "dropped_count": len(part3_keys),
    }
