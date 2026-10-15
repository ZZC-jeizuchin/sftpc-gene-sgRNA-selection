# 本项目由华南师范大学附属中学知识城校区高二一班同学制作
# Copyright (C) 2026 华南师范大学附属中学知识城校区高二一班
# SPDX-License-Identifier: AGPL-3.0-or-later
"""三个数据合并 + 加权排序。

权重 7:6:4（Doench'16 均匀映射到 0-1 / CFD 除以 100 / 距离分档+覆盖修正），
    final = (7*d16 + 6*cfd/100 + 4*dist) / 17

按 PAM 起点对齐后分三堆：两站都有、仅 CRISPOR、仅 CHOPCHOP（丢掉）。
仅 CRISPOR 那批的第三项优先真算，算不出才用两站都有的平均值顶替；
两站都有的那批为空时（滚木），第三项全部置 null。
"""

from __future__ import annotations

from genome_utils import (
    DEFAULT_TARGET_POS, bands_json, distance_band_score, normalize_strand,
)

WEIGHTS = {"doench16": 7, "cfd_spec": 6, "distance": 4}
TOTAL_WEIGHT = sum(WEIGHTS.values())          # 17

# 覆盖奖励：20nt 区间压住目标点则 +0.1，没压住则 -0.1。
# 减分是必须的 —— 压住目标的切点必然 <20bp、分档本来就是 1.0，
# 只加不减的话 1.0+0.1 会被上限吃掉，加分毫无作用。
COVERAGE_BONUS = 0.1      # 压住目标点 -> +
COVERAGE_PENALTY = 0.1    # 没压住     -> -
THIRD_ITEM_CAP = 1.0
# 下限取最低档 0.1，不取 0：扣分不能把 >200bp 那档压成 0 分。
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


# 覆盖判定用哪个窗口："guide_start"（X <= 目标 <= X+20）或 "spacer"
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

    压住目标点 -> +bonus；没压住 -> -penalty。
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


# 硬性过滤规则（6 条），默认关闭，用 merge(hard_filter=True) 开启
HARD_FILTER_RULES = [
    {"id": "length",  "name": "spacer 长度",
     "rule": "必须恰好 20 nt",              "why": "长度 != 20 nt",
     "note": "SpCas9 的 spacer 固定 20 nt，多一个少一个都不行。"},
    {"id": "pam",     "name": "PAM 序列",
     "rule": "必须是 NGG（SpCas9 标准）",    "why": "不包含 NGG",
     "note": "只有 NGG 才是 SpCas9 的标准 PAM。"},
    {"id": "gc",      "name": "GC 含量",
     "rule": "必须介于 40% 到 70%",          "why": "GC < 40% 或 GC > 70%",
     "note": "GC 太低结合不稳，太高容易形成二级结构、脱靶也高。"},
    {"id": "tttt",    "name": "终止信号",
     "rule": "spacer 中不能含连续 4 个 T（TTTT）", "why": "包含 TTTT",
     "note": "TTTT 是转录终止信号，会导致 sgRNA 提前截断。"},
    {"id": "repeat",  "name": "重复序列",
     "rule": "不能含连续 4 个以上相同碱基（如 AAAA）", "why": "连续重复 >= 4",
     "note": "连续同碱基会降低合成质量，也容易滑链错配。"},
    {"id": "dup",     "name": "去重",
     "rule": "相同 spacer 序列只保留一条",   "why": "重复序列",
     "note": "同一个 20 nt spacer 出现在不同位置时只留一条。"},
]


def gc_content(spacer: str | None) -> float | None:
    """spacer 的 GC 百分比（0-100）。"""
    if not spacer:
        return None
    s = spacer.upper()
    return round(100.0 * sum(1 for c in s if c in "GC") / len(s), 4)


def hard_filter_violations(row: dict) -> list[str]:
    """返回这条 guide 违反了哪些硬性规则（空列表 = 通过）。"""
    seq = (row.get("sequence") or "").upper()
    spacer = (row.get("spacer") or seq[:20]).upper()
    bad: list[str] = []

    if len(spacer) != 20:
        bad.append("length")
    if not (len(seq) >= 23 and seq[-2:] == "GG"):
        bad.append("pam")
    g = gc_content(spacer)
    if g is None or not (40.0 <= g <= 70.0):
        bad.append("gc")
    if "TTTT" in spacer:
        bad.append("tttt")
    for base in "ACGT":
        if base * 4 in spacer:
            bad.append("repeat")
            break
    return bad


def _count_reasons(rows: list[dict]) -> dict:
    """统计每条硬性规则各刷掉了多少条。"""
    from collections import Counter
    c = Counter()
    for r in rows:
        for k in r.get("filtered_by", []):
            c[k] += 1
    return dict(c)


def merge(crispor_payload: dict | None, chopchop_payload: dict | None,
          target_pos: int = DEFAULT_TARGET_POS,
          third_mode: str = "auto",
          hard_filter: bool = False) -> dict:
    """合并 + 分三部分 + 加权排序。

    third_mode: "auto"（默认，能真算就真算）/"average"/"computed"
    hard_filter: True 则先用 6 条硬性规则过滤再评分排序
    """
    # 1. 建索引
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

    # 2. 分三部分
    part1_keys = sorted(set(cri) & set(cho), key=str)      # 两站都有
    part2_keys = sorted(set(cri) - set(cho), key=str)      # cris 有、跳跳没有
    part3_keys = sorted(set(cho) - set(cri), key=str)      # 跳跳有、cris 没有 -> 扔掉

    # 3. 组行
    def make_row(key, kind):
        c = cri.get(key, {})
        h = cho.get(key, {})
        src = c or h
        p = c.get("pam_start") if c.get("pam_start") is not None else h.get("pam_start")
        seqp = c.get("seq_pam_start") if c.get("seq_pam_start") is not None else h.get("seq_pam_start")
        st = normalize_strand(src.get("strand"))
# sgRNA 起点 = 23mer 最左端（+ 链为 PAM起点-20，- 链就是 PAM 起点）
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
            "gc_content": gc_content(spacer_of(src.get("target_seq"))),
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

    # 3.5 硬性过滤（可选）
    filtered_rows: list[dict] = []
    if hard_filter:
        kept: list[dict] = []
        seen_spacer: set[str] = set()
        for r in all_rows:
            bad = hard_filter_violations(r)
            sp = (r.get("spacer") or "").upper()
            if not bad and sp in seen_spacer:
                bad = ["dup"]                       # 同 spacer 只留第一条
            if bad:
                r["filtered_by"] = bad
                filtered_rows.append(r)
                continue
            seen_spacer.add(sp)
            kept.append(r)
        all_rows = kept
        rows1 = [r for r in all_rows if r["part"] == 1]
        rows2 = [r for r in all_rows if r["part"] == 2]

    # 4. 第三项
    # 第一部分：分档 + 覆盖奖励
    for r in rows1:
        v, base, extra = third_item_value(r["distance"], r["pam_start"], r["strand"], target_pos,
                                              guide_start=r["guide_start"])
        r["covers_target"] = bool(extra and extra > 0)
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
        r["covers_target"] = bool(extra and extra > 0)
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

    # 5. 滚木
    #   第一部分为空 = 两个网站完全不重合 -> 第三项全部 = null
    kurumi = (len(rows1) == 0)
    if kurumi:
        for r in all_rows:
            r["third_item"] = None
            r["third_source"] = None
        note = ("两个网站没有找到任何共同的 sgRNA（「两站都有」那批为空）-> "
                "突变距离评分 = 滚木(null)。请确认输入的这条基因是否过于变异，"
                "或者根本不是目标基因。")
    else:
        note = None

    # 6. 其他两项
    d16 = uniform_map([r["doench16"] for r in all_rows], higher_is_better=True)
    cfd_norm = [(None if r["cfd_spec"] is None else round(r["cfd_spec"] / 100.0, 6))
                for r in all_rows]

    # 7. 加权
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
        # 硬性过滤
        "hard_filter": hard_filter,
        "hard_filter_rules": HARD_FILTER_RULES,
        "filtered_count": len(filtered_rows),
        "filtered_before": len(all_rows) + len(filtered_rows),
        "filtered_reasons": _count_reasons(filtered_rows),
        "avg_part1_third": avg_part1,
        "kurumi": kurumi,
        "kurumi_note": note,
        "distance_available": any(r["third_item"] is not None for r in out),
        "count": len(scored),
        "rows": scored,
        "unscored": [r for r in out if r["final_score"] is None],
        "dropped_count": len(part3_keys),
    }
