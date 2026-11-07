# 本项目由华南师范大学附属中学知识城校区高二一班同学制作
# Copyright (C) 2026 华南师范大学附属中学知识城校区高二一班
# SPDX-License-Identifier: AGPL-3.0-or-later
"""第四维度：等位基因特异性评分（Allele-Specificity Score, ASS）。

为什么需要这一维
----------------
原来的三项（切割效率 / 脱靶特异性 / 突变距离）都是**单序列**评分：
它们只看参考基因组那条序列。但 SFTPC 相关 chILD 是常染色体显性遗传，
患者只带一个突变等位基因，真正要回答的问题是：

    「这条 sgRNA 能不能只切突变型、不切野生型？」

「切点离突变近」只是这个问题的**代理指标**：切点很近但两条等位基因都切的
sgRNA，做等位基因特异性编辑毫无价值。所以本模块把评分从「单序列」扩展到
「野生型 vs 突变型双序列比较」。

ΔScore（双序列切割效率差）
--------------------------
    原始差 Δ_raw = 切割效率(突变型) − 切割效率(野生型)     范围 −100 … +100

**关键：Δ_raw 必须再归一化到 0–100 才能当 ΔScore 用。**

CRISPOR 的 Doench '16 是 0–100 的整数分，单个碱基错配只会让它变动几个点
（本项目实测：434 条 guide 里只有 7 条 Δ_raw ≠ 0，最大 10 分，中位数 1 分）。
如果直接把 Δ_raw 当成 0–100 的特异性分数，那么「ΔScore < 20 就降权」这条
规则会对**所有** guide 生效，四个维度里 35% 的那一维直接归零，总分被压在
65 分以下，等级表（90/80/70/60）永远走不到。

所以这里用一条固定锚：**Doench '16 差 20 分 = ΔScore 满分**。

    ΔScore = 100 × clamp(Δ_raw, 0, 20) / 20

20 分这个锚的含义是：突变型比野生型的预测切割效率高出 20 分，
已经是很大的量级差，足以认为这条 sgRNA「明显偏爱突变型」。
换成 CRISPRseek 那种「0 = 两条都切、100 = 只切一条」的口径后，
降权阈值 20 才重新有意义（ΔScore 20 ⇔ Δ_raw 4 分）。

ΔScore 越大，说明这条 sgRNA 越偏向切突变型。计算方式：把**同一条输入序列**
和**只改了一个碱基的突变型序列**分别提交给 CRISPOR，拿到每条 sgRNA 的
Doench '16 评分，再按位置配对求差。

这样做的好处是双序列用的是**同一套模型、同一个量纲**，不需要自己拟合。
而且三条边界情况会被自然分开：

    两条都出现      Δ_raw = Doench(突变型) − Doench(野生型)
    只出现在突变型  突变造出了新 PAM —— 最强的特异性，Δ_raw = +100
    只出现在野生型  突变破坏了 PAM —— Δ_raw = −100（夹到 0）

如果这条 sgRNA 的 23mer 根本不覆盖突变位点，它在两条序列上**完全相同**，
Δ_raw 必然是 0 —— 这在生物学上也是对的：不覆盖突变，就不可能区分两条等位基因。

位置分（Position score）
------------------------
突变落在 23mer 的哪个位置，决定了错配被 SpCas9 感知的强度。位置编号统一按
PAM 近端为 1 起算（CRISPOR 的惯例）：

    PAM 区   位置 21–23   100 分   突变直接把 NGG 改掉，区分最彻底
    Seed 区  位置  1– 8    80 分   错配显著破坏配对稳定性
    中间区   位置  9–14    40 分
    PAM 远端 位置 15–20    10 分   错配耐受性高
    不覆盖    —             0 分

合成
----
    ASS = 0.6 × clamp(ΔScore, 0, 100) + 0.4 × 位置分

    降权：ΔScore < 20（两条切得差不多）时整体 ×0.3，
          避免「距离近但两条都切」的候选靠位置分虚高。
"""

from __future__ import annotations

from genome_utils import normalize_strand

# 突变定义（hg38，1-based；与 突变坐标信息173T.xls 一致）
MUTATION = {
    "gene": "SFTPC",
    "hgvs_c": "c.218T>C",
    "hgvs_p": "p.Ile73Thr (I73T)",
    "hgvs_g": "chr8:g.22163096T>C",
    "chrom": "chr8",
    "pos": 22_163_096,
    "ref": "T",
    "alt": "C",
    "codon_wt": "ATT",     # Ile
    "codon_mut": "ACT",    # Thr
}

# 位置分表：键 -> (位置下界, 位置上界, 分数, 中文名)
# 位置按「距 PAM 近端第几位」编号，1 = 紧邻 PAM 的 spacer 碱基
POSITION_TABLE: dict[str, dict] = {
    "pam":  {"lo": 21, "hi": 23, "score": 100.0, "label": "PAM 区"},
    "seed": {"lo": 1,  "hi": 8,  "score": 80.0,  "label": "Seed 区"},
    "mid":  {"lo": 9,  "hi": 14, "score": 40.0,  "label": "中间区"},
    "far":  {"lo": 15, "hi": 20, "score": 10.0,  "label": "PAM 远端"},
    "none": {"lo": None, "hi": None, "score": 0.0, "label": "不覆盖突变"},
}

# 合成权重
DELTA_WEIGHT = 0.6          # ΔScore 占 60%
POSITION_WEIGHT = 0.4       # 位置分 占 40%
DOWNGRADE_DELTA = 20.0      # ΔScore（0-100）低于此值视为「两条都切」
DOWNGRADE_FACTOR = 0.3      # 降权系数
ASS_MAX = 100.0

# Δ_raw -> ΔScore 的锚：Doench '16 差 20 分算满分
DELTA_FULL_SCALE = 20.0


# --------------------------------------------------------------------------
# 双序列构造
# --------------------------------------------------------------------------

def mutation_index(wt_sequence: str, genome_offset: int | None,
                   target_pos: int = MUTATION["pos"]) -> int:
    """突变位点在输入序列里的 0-based 下标。

    genome_offset 是输入序列的基因组 0-based 起点（与 genome_utils 一致），
    所以：下标 = (target_pos − 1) − genome_offset。
    """
    if genome_offset is None:
        raise ValueError("需要 genome_offset 才能把突变位点换算到序列坐标")
    idx = (target_pos - 1) - int(genome_offset)
    if not (0 <= idx < len(wt_sequence)):
        raise ValueError(
            f"突变位点 {target_pos} 不在输入序列范围内"
            f"（序列覆盖 {genome_offset + 1}–{genome_offset + len(wt_sequence)}）")
    return idx


def check_ref(wt_sequence: str, idx: int, ref: str = MUTATION["ref"]) -> bool:
    """参考碱基是否与突变记录一致（防坐标翻车）。"""
    return wt_sequence[idx].upper() == ref.upper()


def build_mutant(wt_sequence: str, idx: int,
                 ref: str = MUTATION["ref"], alt: str = MUTATION["alt"]) -> str:
    """把野生型序列的 idx 位由 ref 换成 alt，返回突变型序列。"""
    wt = wt_sequence.upper()
    base = wt[idx]
    if base != ref.upper():
        raise ValueError(f"第 {idx} 位是 {base}，不是预期的 {ref}——坐标可能不对")
    return wt[:idx] + alt.upper() + wt[idx + 1:]


# --------------------------------------------------------------------------
# 几何：23mer 区间 / 突变在 23mer 里的位置
# --------------------------------------------------------------------------

def guide_window(seq_pam_start: int, strand: str) -> tuple[int, int]:
    """guide 23mer 在**输入序列坐标**（1-based）里的闭区间。

    正链：23mer = [PAM起点−20, PAM起点+2]
    反链：23mer = [PAM起点, PAM起点+22]
    """
    s = normalize_strand(strand)
    p = int(seq_pam_start)
    if s == "+":
        return p - 20, p + 2
    if s == "-":
        return p, p + 22
    raise ValueError(f"未知链向：{strand!r}")


def index_in_23mer(seq_pam_start: int, strand: str, seq_target: int) -> int | None:
    """突变位点是 23mer 里第几个碱基（1-based，从 target_seq 的 5' 端数起）。

    正链：第 i 位对应序列坐标 (PAM起点 − 21 + i)
    反链：第 i 位对应序列坐标 (PAM起点 + 23 − i)   —— 因为报的是反向互补
    """
    s = normalize_strand(strand)
    p = int(seq_pam_start)
    m = int(seq_target)
    lo, hi = guide_window(p, s)
    if not (lo <= m <= hi):
        return None
    return (m - p + 21) if s == "+" else (p + 23 - m)


def position_from_pam(i: int) -> int:
    """23mer 第 i 位 -> 「距 PAM 近端第几位」（1 = 紧邻 PAM）。"""
    return 21 - i


def classify(i: int | None) -> str:
    """23mer 第 i 位落在哪个功能区。"""
    if i is None:
        return "none"
    if 21 <= i <= 23:
        return "pam"
    p = position_from_pam(i)
    if 1 <= p <= 8:
        return "seed"
    if 9 <= p <= 14:
        return "mid"
    if 15 <= p <= 20:
        return "far"
    return "none"


def position_detail(seq_pam_start, strand, seq_target) -> dict:
    """一次算齐：23mer 位置、距 PAM 位数、功能区、位置分。"""
    i = index_in_23mer(seq_pam_start, strand, seq_target)
    region = classify(i)
    return {
        "index_in_23mer": i,
        "pos_from_pam": (None if i is None else position_from_pam(i)),
        "region": region,
        "region_label": POSITION_TABLE[region]["label"],
        "position_score": POSITION_TABLE[region]["score"],
    }


# --------------------------------------------------------------------------
# ΔScore：双序列配对
# --------------------------------------------------------------------------

def _key(v: dict):
    """配对键：(链向, PAM 在输入序列坐标系里的起点)。两条序列等长，可直接对。"""
    return (normalize_strand(v.get("strand")), v.get("seq_pam_start"))


def build_delta_index(wt_values: list[dict], mut_values: list[dict],
                      score_key: str = "doench16") -> dict:
    """把两次运行的 guide 按位置配对，算出每条 guide 的 ΔScore。

    返回 {key: {"delta", "efficiency_wt", "efficiency_mut", "delta_source"}}。
    """
    wt = {_key(v): v for v in wt_values if _key(v)[1] is not None}
    mut = {_key(v): v for v in mut_values if _key(v)[1] is not None}

    out: dict = {}
    for k in set(wt) | set(mut):
        a, b = wt.get(k), mut.get(k)
        sw = None if a is None else a.get(score_key)
        sm = None if b is None else b.get(score_key)
        src = "both"
        if a is not None and b is None:
            # 只出现在野生型：突变破坏了 PAM
            out[k] = {"delta": -100.0, "efficiency_wt": sw, "efficiency_mut": None,
                      "delta_source": "wt_only(pam_lost)"}
            continue
        if b is not None and a is None:
            # 只出现在突变型：突变造出了新 PAM
            out[k] = {"delta": 100.0, "efficiency_wt": None, "efficiency_mut": sm,
                      "delta_source": "mut_only(pam_gained)"}
            continue
        if sw is None or sm is None:
            out[k] = {"delta": None, "efficiency_wt": sw, "efficiency_mut": sm,
                      "delta_source": "missing"}
            continue
        out[k] = {"delta": round(float(sm) - float(sw), 4),
                  "efficiency_wt": sw, "efficiency_mut": sm,
                  "delta_source": src}
    return out


# --------------------------------------------------------------------------
# 合成
# --------------------------------------------------------------------------

def clamp(x: float, lo: float = 0.0, hi: float = 100.0) -> float:
    return max(lo, min(hi, x))


def delta_to_100(delta_raw: float | None, full_scale: float = DELTA_FULL_SCALE) -> float:
    """把 Doench '16 的原始差值归一化成 0–100 的 ΔScore。

    Δ_raw ≤ 0        -> 0      （偏野生型，对等位基因特异性没有正面价值）
    Δ_raw ≥ full_scale -> 100  （对突变型有决定性偏好）
    """
    d = 0.0 if delta_raw is None else float(delta_raw)
    return round(clamp(d, 0.0, full_scale) / full_scale * 100.0, 4)


def combine(delta: float | None, region: str) -> tuple[float, dict]:
    """ASS = 0.6 × ΔScore + 0.4 × 位置分，ΔScore < 20 时整体 ×0.3。

    delta 是**原始** Doench 差值；内部先归一化成 0–100 的 ΔScore 再合成。
    """
    pos = POSITION_TABLE.get(region, POSITION_TABLE["none"])["score"]
    d_raw = 0.0 if delta is None else float(delta)
    d100 = delta_to_100(d_raw)

    base = DELTA_WEIGHT * d100 + POSITION_WEIGHT * pos
    downgraded = d100 < DOWNGRADE_DELTA
    final = base * (DOWNGRADE_FACTOR if downgraded else 1.0)
    return round(clamp(final, 0.0, ASS_MAX), 4), {
        "delta_raw": round(d_raw, 4),
        "delta_100": d100,
        "position_score": pos,
        "ass_raw": round(base, 4),
        "downgraded": downgraded,
    }


def score_guide(seq_pam_start, strand, seq_target, delta: float | None,
                assigned_mutation: bool) -> dict:
    """给一条 guide 打 ASS 全部分项。"""
    det = position_detail(seq_pam_start, strand, seq_target)
    ass, extra = combine(delta, det["region"])
    return {
        "ass": ass,
        "delta": None if delta is None else round(float(delta), 4),
        **det,
        **extra,
        "mutation_window": bool(assigned_mutation),
    }


def weights_json() -> dict:
    return {
        "delta_weight": DELTA_WEIGHT,
        "position_weight": POSITION_WEIGHT,
        "downgrade_delta": DOWNGRADE_DELTA,
        "downgrade_factor": DOWNGRADE_FACTOR,
        "delta_full_scale": DELTA_FULL_SCALE,
    }


def position_bands_json() -> list[dict]:
    return [
        {"region": k, "lo": v["lo"], "hi": v["hi"], "score": v["score"], "label": v["label"]}
        for k, v in POSITION_TABLE.items()
    ]
