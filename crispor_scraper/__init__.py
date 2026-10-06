"""CRISPOR（网站1）数据抓取包。

对外只暴露三层：
    config   —— 定义 Step1/Step2/Step3 三个输入
    client   —— 与 CRISPOR 的 CGI 接口通信（提交 / 轮询 / 下载）
    parser   —— 把 TSV 解析成结构化记录

典型用法：
    from crispor_scraper import CrisporInput, CrisporClient, fetch
    guides = fetch(CrisporInput(sequence="ACGT...", genome="hg38", pam="NGG"))
    for g in guides:
        print(g.guide_id, g.doench16, g.cfd_spec)
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
