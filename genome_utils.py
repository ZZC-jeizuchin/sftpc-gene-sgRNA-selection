"""两个爬虫共用的基因组/坐标工具。

为什么要共用：
    CRISPOR 和 CHOPCHOP 的 "位置" 列用的是**同一个约定** ——
    23 bp 靶位点（20nt spacer + PAM）在正链上的最左端 1-based 坐标。
    所以切点公式、距离分档、序列定位这三件事必须两边完全一致，
    不能各写一份（否则又会像之前那样差 11 bp）。

提供：
    locate_sequence()   裸碱基序列 -> 基因组 offset（本地匹配，不依赖 BLAT）
    cut_site()          位置 + 链 -> 切割位点（正链坐标）
    distance_to()       |切点 - 目标位点|
    distance_band_score() / band_name()   文档规定的区间给分
"""

from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

# ============================================================ 常量

# 目标位点：SFTPC c.218 / I73T，hg38 chr8:22,163,096
DEFAULT_TARGET_POS = 22_163_096

# 文档规定的距离分档：(上界不含, 分数)
DISTANCE_BANDS: list[tuple[float, float]] = [
    (20.0, 1.0),
    (50.0, 0.9),
    (100.0, 0.7),
    (200.0, 0.4),
    (float("inf"), 0.1),
]

# ⚠️ 两个网站的"位置"列锚点**不同**，这是最容易错的地方：
#   CHOPCHOP  Genomic location = 23mer 在正链上的最左端（anchor="23mer"）
#   CRISPOR   Position/Strand  = **PAM 区间**在正链上的左端（anchor="pam"）
#   fw / + 链上两者相差 20 bp。
#
# SpCas9 切在 PAM 5' 侧 3 bp（protospacer 第 17/18 位之间），用 PAM 起点 p 表示：
#   正链: 切点 = p - 3.5      （protospacer 在 p-20..p-1）
#   反链: 切点 = p + 5.5      （protospacer 在 p+3..p+22）
CUT_FROM_PAM = {"+": -3.5, "-": 5.5}
# 23mer 锚点 -> PAM 起点：+ 链 PAM 在右端(+20)，- 链 PAM 在左端(+0)
PAM_FROM_23MER = {"+": 20, "-": 0}

_UCSC_SEQ_API = "https://api.genome.ucsc.edu/getData/sequence"

# 默认搜索窗口：SFTPC 基因座（hg38 chr8 ~22.15–22.17 Mb）
DEFAULT_LOCUS = {"genome": "hg38", "chrom": "chr8", "start": 22_150_000, "end": 22_170_000}

# 基因座序列的本地缓存目录（两个爬虫共用一份，不要在 CWD 里乱丢 locus/）
DEFAULT_LOCUS_CACHE = ".genome_cache"


# ============================================================ 坐标换算

def normalize_strand(strand: str) -> str | None:
    if strand is None:
        return None
    s = str(strand).strip().lower()
    if s in ("+", "fw", "fwd", "forward", "1"):
        return "+"
    if s in ("-", "rev", "reverse", "-1"):
        return "-"
    return None


def pam_start(location, strand, anchor: str = "23mer"):
    """把报告值统一成 **PAM 区间的左端坐标**。

    anchor="23mer" -> CHOPCHOP 约定；anchor="pam" -> CRISPOR 约定
    """
    if location is None:
        return None
    key = normalize_strand(strand)
    if key is None:
        return None
    if anchor == "pam":
        return location
    if anchor == "23mer":
        return location + PAM_FROM_23MER[key]
    raise ValueError(f"未知 anchor：{anchor!r}（应为 '23mer' 或 'pam'）")


def cut_site(location, strand, chrom: str | None = None, offset: int | None = None,
             anchor: str = "23mer"):
    """算切割位点（正链坐标）。

    location : 网站报告的位置
    strand   : "+"/"fw" 或 "-"/"rev"
    anchor   : "23mer"（CHOPCHOP）或 "pam"（CRISPOR）
    chrom    : "seq" 表示序列坐标（CHOPCHOP 贴序列模式），需配合 offset
    换算：基因组坐标 = offset + location
    """
    genomic = to_genomic(location, chrom, offset)
    if genomic is None:
        return None
    p = pam_start(genomic, strand, anchor)
    key = normalize_strand(strand)
    return None if p is None else p + CUT_FROM_PAM[key]


def to_genomic(location, chrom: str | None, offset: int | None):
    if location is None:
        return None
    try:
        loc = int(location)
    except (TypeError, ValueError):
        return None
    if chrom and str(chrom).lower() != "seq":
        return loc
    return None if offset is None else offset + loc


def distance_to(location, strand, target_pos: int = DEFAULT_TARGET_POS,
                chrom: str | None = None, offset: int | None = None,
                anchor: str = "23mer"):
    """|切点 − 目标位点|，单位 bp。"""
    cut = cut_site(location, strand, chrom, offset, anchor)
    return None if cut is None else abs(cut - target_pos)


def distance_band_score(distance):
    """按文档分档给分。"""
    if distance is None:
        return None
    for upper, score in DISTANCE_BANDS:
        if distance < upper:
            return score
    return DISTANCE_BANDS[-1][1]


def band_name(score):
    if score is None:
        return None
    return {1.0: "<=20bp", 0.9: "20-50bp", 0.7: "50-100bp",
            0.4: "100-200bp", 0.1: ">200bp"}.get(score)


def bands_json() -> list[dict]:
    return [{"lt": (None if u == float("inf") else u), "score": s} for u, s in DISTANCE_BANDS]


# ============================================================ 序列定位

def _revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def fetch_locus(genome: str, chrom: str, start: int, end: int,
                cache_dir: str | Path = DEFAULT_LOCUS_CACHE) -> str:
    cache = Path(cache_dir) / "locus"
    cache.mkdir(parents=True, exist_ok=True)
    path = cache / f"{genome}_{chrom}_{start}_{end}.txt"
    if path.is_file():
        return path.read_text(encoding="utf-8")

    url = f"{_UCSC_SEQ_API}?" + urllib.parse.urlencode(
        {"genome": genome, "chrom": chrom, "start": start, "end": end}
    )
    last = None
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as resp:
                dna = json.loads(resp.read().decode())["dna"].upper()
            path.write_text(dna, encoding="utf-8")
            return dna
        except Exception as exc:
            last = exc
            time.sleep(2 ** attempt)
    raise RuntimeError(f"拉取基因座失败：{last}")


def _mismatches(a: str, b: str, limit: int) -> int:
    n = 0
    for x, y in zip(a, b):
        if x != y:
            n += 1
            if n > limit:
                return n
    return n


def locate_sequence(sequence: str,
                    genome: str = DEFAULT_LOCUS["genome"],
                    chrom: str = DEFAULT_LOCUS["chrom"],
                    start: int = DEFAULT_LOCUS["start"],
                    end: int = DEFAULT_LOCUS["end"],
                    max_mismatch: int = 5,
                    cache_dir: str | Path = DEFAULT_LOCUS_CACHE) -> dict:
    """裸碱基序列 -> 基因组位置。

    返回 {"found", "chrom", "offset"(0-based), "strand", "mismatches", "note"}
    换算关系：基因组坐标 = offset + N（N = 序列内 1-based 位置）
    """
    seq = "".join(str(sequence).split()).upper()
    if not seq:
        return {"found": False, "note": "空序列"}

    try:
        locus = fetch_locus(genome, chrom, start, end, cache_dir)
    except Exception as exc:
        return {"found": False, "note": f"拉取基因座失败：{exc}"}

    # 1) 精确匹配
    rc = _revcomp(seq)
    i = locus.find(seq)
    if i >= 0:
        return {"found": True, "chrom": chrom, "offset": start + i,
                "strand": "+", "mismatches": 0, "note": "精确匹配（正链）"}
    i = locus.find(rc)
    if i >= 0:
        return {"found": True, "chrom": chrom, "offset": start + i,
                "strand": "-", "mismatches": 0, "note": "精确匹配（反链）"}

    # 2) 模糊匹配：切成 max_mismatch+1 段，至少一段完全干净
    if len(seq) < (max_mismatch + 1) * 8:
        return {"found": False, "note": "序列太短，无法在允许错配下定位"}

    for query, strand in ((seq, "+"), (rc, "-")):
        n_parts = max_mismatch + 1
        seg = len(query) // n_parts
        for k in range(n_parts):
            s0 = k * seg
            s1 = len(query) if k == n_parts - 1 else (k + 1) * seg
            seed = query[s0:s1]
            if len(seed) < 8:
                continue
            pos = locus.find(seed)
            while pos >= 0:
                off = pos - s0
                if 0 <= off <= len(locus) - len(query):
                    mm = _mismatches(query, locus[off: off + len(query)], max_mismatch)
                    if mm <= max_mismatch:
                        return {"found": True, "chrom": chrom, "offset": start + off,
                                "strand": strand, "mismatches": mm,
                                "note": f"模糊匹配（{mm} 个错配，{strand} 链）"}
                pos = locus.find(seed, pos + 1)

    return {"found": False,
            "note": f"在 {chrom}:{start}-{end} 内找不到这段序列（错配上限 {max_mismatch}）"}
