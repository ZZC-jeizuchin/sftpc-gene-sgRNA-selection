# 本项目由华南师范大学附属中学知识城校区高二一班同学制作
# Copyright (C) 2026 华南师范大学附属中学知识城校区高二一班
# SPDX-License-Identifier: AGPL-3.0-or-later
"""CRISPOR（网站1）抓取包。

    from crispor_scraper import CrisporInput, fetch
    r = fetch(CrisporInput(sequence="ACGT...", genome="hg38", pam="NGG"))
"""
# 让本包能 import 到工作区根目录的 genome_utils（两个爬虫共用一套坐标逻辑）
import sys as _sys
import pathlib as _pathlib

_ROOT = _pathlib.Path(__file__).resolve().parents[1]
if str(_ROOT) not in _sys.path:
    _sys.path.insert(0, str(_ROOT))


from .config import CrisporInput, CrisporConfigError
from .parser import GuideRecord, parse_guides_tsv
from .cache import CrisporCache, DEFAULT_CACHE_DIR
from .client import CrisporClient, CrisporError
from .pipeline import fetch, CrisporResult
from .datasets import (
    DATASETS, DEFAULT_ORDER, DOENCH16, CFD_SPEC, DatasetSpec, resolve, build_payload,
)
from .distance import build_distance_payload, cut_site, distance_band_score
from . import jsonio

__all__ = [
    "CrisporInput",
    "CrisporConfigError",
    "GuideRecord",
    "parse_guides_tsv",
    "CrisporCache",
    "DEFAULT_CACHE_DIR",
    "CrisporClient",
    "CrisporError",
    "fetch",
    "CrisporResult",
    "DatasetSpec",
    "DATASETS",
    "DEFAULT_ORDER",
    "DOENCH16",
    "CFD_SPEC",
    "resolve",
    "build_payload",
    "jsonio",
]
