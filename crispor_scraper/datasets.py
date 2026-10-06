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
"""两个数据的**独立**定义。

每个数据是一个 DatasetSpec：能单独取、单独算、单独输出成 JSON，互不依赖。
以后加第三个数据（距离）只需在这里加一个 spec，不用动 client / parser / cache。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

from .parser import GuideRecord


@dataclass(frozen=True)
class DatasetSpec:
    """一个"数据"的完整定义。"""

    id: str                     # 程序内标识，如 "doench16"
    name: str                   # 网站上的列名，如 "Doench '16"
    source: str                 # 来源网站，如 "crispor"
    source_url: str
    source_column: str          # CRISPOR TSV 里的列名
    description: str
    value_type: str             # "int" / "float"
    value_min: float | None
    value_max: float | None
    direction: str              # "higher_is_better" / "lower_is_better"
    extract: Callable[[GuideRecord], object]

    def value_of(self, record: GuideRecord):
        return self.extract(record)

    def describe(self) -> dict:
        return {
            "id": self.id,
            "name": self.name,
            "source": self.source,
            "source_url": self.source_url,
            "source_column": self.source_column,
            "description": self.description,
            "value_type": self.value_type,
            "range": [self.value_min, self.value_max],
            "direction": self.direction,
        }


# --------------------------------------------------------------------- 数据 1

DOENCH16 = DatasetSpec(
    id="doench16",
    name="Doench '16",
    source="crispor",
    source_url="https://crispor.gi.ucsc.edu/",
    source_column="Doench '16-Score",
    description=(
        "on-target 编辑效率预测。CRISPOR 内部叫 fusi，调的是微软官方 Azimuth 2.0 "
        "（Doench et al. 2016 Rule Set 2，梯度提升回归树）。输入 30mer(-4/+3)，"
        "原始输出 0-1，CRISPOR 已缩放到 0-100 的整数。"
        "适用：U6 启动子、细胞内表达。序列两端 50bp 内的 guide 或含 N 的算不出（值为 null）。"
    ),
    value_type="int",
    value_min=0,
    value_max=100,
    direction="higher_is_better",
    extract=lambda r: _score(r.doench16),
)

# --------------------------------------------------------------------- 数据 2

CFD_SPEC = DatasetSpec(
    id="cfd_spec",
    name="CFD Spec. score",
    source="crispor",
    source_url="https://crispor.gi.ucsc.edu/",
    source_column="cfdSpecScore",
    description=(
        "off-target 特异性预测。逐个脱靶位点算 CFD（Doench 2016 的位置×错配类型权重矩阵 "
        "× PAM 因子），再用 MIT 那套公式聚合：round(10000/(100+ΣCFD))。"
        "注意这是**guide 级聚合特异性分**（越大越安全），不是单个脱靶位点的 CFD 切割分"
        "（那个是 0-1、越大越危险）。全基因组搜到 4 个错配，NGG guide 另计入 NAG/NGA PAM。"
    ),
    value_type="int",
    value_min=0,
    value_max=100,
    direction="higher_is_better",
    extract=lambda r: _score(r.cfd_spec),
)


# --------------------------------------------------------------------- 注册表

def _score(value):
    """CRISPOR 用负数当"算不出"的哨兵值（序列含 N、Cpf1、重复序列等）。
    统一转成 None，让 JSON 里是 null，统计时计入 missing。"""
    if value is None:
        return None
    try:
        v = float(value)
    except (TypeError, ValueError):
        return None
    return None if v < 0 else value


DATASETS: dict[str, DatasetSpec] = {d.id: d for d in (DOENCH16, CFD_SPEC)}

# 默认输出顺序
DEFAULT_ORDER: list[str] = ["doench16", "cfd_spec"]


def resolve(ids: list[str] | None) -> list[DatasetSpec]:
    """把 id 列表解析成 spec 列表；None 表示用默认顺序。"""
    if not ids:
        return [DATASETS[i] for i in DEFAULT_ORDER]
    unknown = [i for i in ids if i not in DATASETS]
    if unknown:
        raise KeyError(
            f"未知的数据 id {unknown}；可用的有 {sorted(DATASETS)}"
        )
    # 去重并保持请求顺序
    seen, out = set(), []
    for i in ids:
        if i not in seen:
            seen.add(i)
            out.append(DATASETS[i])
    return out


def build_payload(spec: DatasetSpec, records: list[GuideRecord]) -> dict:
    """把一个数据单独打包成 JSON 结构：元信息 + 逐 guide 的取值。"""
    values = []
    for r in records:
        values.append({
            "guide_id": r.guide_id,
            "strand": r.strand,
            "cut_position": r.cut_position,     # CRISPOR Position/Strand 的数字部分
            "target_seq": r.target_seq,         # 20nt spacer + PAM
            "value": spec.value_of(r),
        })
    numeric = [v["value"] for v in values if isinstance(v["value"], (int, float))]
    return {
        **spec.describe(),
        "count": len(values),
        "stats": {
            "min": min(numeric) if numeric else None,
            "max": max(numeric) if numeric else None,
            "mean": round(sum(numeric) / len(numeric), 2) if numeric else None,
            "missing": len(values) - len(numeric),
        },
        "values": values,
    }
