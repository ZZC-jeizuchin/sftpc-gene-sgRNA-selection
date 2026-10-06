"""CHOPCHOP（网站2）数据抓取包。

用途：抓第三个数据 —— **切割位点到 chr8:22,163,096 的距离**（SFTPC c.218 / I73T）。

CHOPCHOP 有正规 JSON API（比网站1 还干净）：
    POST /                             {"opts":[...], "geneInput":"SFTPC", ...}
    ->  {"jobId": "..."}
    GET  /results/<jobId>/run.info     tab 分隔的运行信息
    GET  /results/<jobId>/results.tsv  结果表（Rank/Target sequence/Genomic location/...）
    GET  /results/<jobId>/cutcoords.json   绘图用坐标
    GET  /results/<jobId>/results.bed

坐标约定（用 hg38 全序列实测 161/161 条通过）：
    "Genomic location" = 23 bp 靶位点（20nt + PAM）在**正链上的最左端 1-based 坐标** loc。
    切割位点（SpCas9，PAM 前 3 bp）：
        + 链： loc + 16 与 loc + 17 之间   →  loc + 16.5
        - 链： loc + 5  与 loc + 6 之间    →  loc + 5.5
"""
# 让本包能 import 到工作区根目录的 genome_utils（两个爬虫共用一套坐标逻辑）
import sys as _sys
import pathlib as _pathlib

_ROOT = _pathlib.Path(__file__).resolve().parents[1]
if str(_ROOT) not in _sys.path:
    _sys.path.insert(0, str(_ROOT))


from .config import ChopchopInput, ChopchopConfigError
from .parser import GuideRecord, parse_results_tsv
from .distance import (
    DEFAULT_TARGET_POS, DISTANCE_BANDS, cut_site, distance_band_score, build_distance_payload,
)
from .cache import ChopchopCache, DEFAULT_CACHE_DIR
from .client import ChopchopClient, ChopchopError
from .pipeline import fetch, ChopchopResult
from . import jsonio

__all__ = [
    "ChopchopInput",
    "ChopchopConfigError",
    "GuideRecord",
    "parse_results_tsv",
    "DEFAULT_TARGET_POS",
    "DISTANCE_BANDS",
    "cut_site",
    "distance_band_score",
    "build_distance_payload",
    "ChopchopCache",
    "DEFAULT_CACHE_DIR",
    "ChopchopClient",
    "ChopchopError",
    "fetch",
    "ChopchopResult",
    "jsonio",
]
