# 本项目由华南师范大学附属中学知识城校区高二一班同学制作
# Copyright (C) 2026 华南师范大学附属中学知识城校区高二一班
# SPDX-License-Identifier: AGPL-3.0-or-later
"""CHOPCHOP（网站2）抓取包 —— 第三个数据：切割位点到 chr8:22163096 的距离。

    POST /                              {"opts": [...]} -> {"jobId": "..."}
    GET  /results/<jobId>/results.tsv   结果表
    GET  /results/<jobId>/run.info      运行信息

「Genomic location」= 23mer 在正链上的最左端坐标，换算见 genome_utils。
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
