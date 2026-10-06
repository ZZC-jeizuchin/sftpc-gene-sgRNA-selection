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
