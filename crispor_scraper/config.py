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
"""Step1 / Step2 / Step3 三个输入的统一定义与校验。

对应 CRISPOR 首页那张表单：

    Step 1  Target Sequence : 序列 / FASTA 文件 / 染色体区间   -> 表单字段 seq
    Step 2  Select a genome : 基因组，如 hg38 / hg19         -> 表单字段 org
    Step 3  Select a PAM    : PAM 类型，如 NGG / NAG         -> 表单字段 pam

另外 name（可选序列名）在表单里是 Step1 的一部分，这里也一起放进 Step1 组。
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field

# CRISPOR 官网硬限制：输入序列 < 2300 bp
MAX_SEQ_LEN = 2300

# 允许的碱基（其余字符 CRISPOR 会自动剔除，但本地先校验更友好）
_ALLOWED_BASES = re.compile(r"^[ACGTNacgtn]+$")

# CRISPOR 接受的染色体区间写法，例：chr8:22,119,000-22,125,000  或 chr8:22119000-22119500:+
_REGION_RE = re.compile(r"^\s*[A-Za-z0-9_.\-]+\s*:\s*[\d,]+\s*-\s*[\d,]+\s*(?::\s*[+-])?\s*$")

# 官网下拉框里常用的 PAM 取值（Step3 的可选清单，不限于这些）
COMMON_PAMS = {
    "NGG": "20bp-NGG - SpCas9 / SpCas9-HF1 / eSpCas9 1.1",
    "NAG": "SpCas9 非典型 PAM",
    "NGA": "SpCas9 非典型 PAM",
    "NNG": "20bp-NNG - Cas9 S. canis",
    "NGN": "20bp-NGN - SpG",
    "NNGT": "20bp-NNGT - Cas9 S. canis - 高效 PAM",
    "NAA": "20bp-NAA - iSpyMacCas9",
    "TTN": "TTN-23bp - hfCas12Max",
    "TNN": "TNN-23bp - hfCas12Max",
    "NGG-22": "NGG-22bp - eSpOT-ON (ePsCas9)",
}

# 官网下拉框里常用的基因组（Step2 的可选清单）
COMMON_GENOMES = {
    "hg38": "Homo sapiens - GRCh38/hg38（默认）",
    "hg19": "Homo sapiens - GRCh37/hg19",
    "mm39": "Mus musculus - GRCm39/mm39",
    "mm10": "Mus musculus - GRCm38/mm10",
    "danRer11": "Danio rerio - GRCz11/danRer11",
    "dm6": "Drosophila melanogaster - dm6",
    "ce11": "C. elegans - ce11",
    "sacCer3": "S. cerevisiae - sacCer3",
}


class CrisporConfigError(ValueError):
    """三个 Step 的输入不合法时抛出。"""


@dataclass
class CrisporInput:
    """一次 CRISPOR 查询的完整输入（= 首页表单的三个 Step）。

    Step1 三选一（按下面的优先级）：
        fasta_path : FASTA 文件路径（取第一条序列）
        region     : 染色体区间字符串，如 "chr8:22,119,000-22,125,000"
        sequence   : 直接给序列字符串
    """

    # ---- Step 1: Target sequence ----
    sequence: str | None = None
    region: str | None = None
    fasta_path: str | None = None
    name: str = "query"          # 可选：输出里显示的序列名

    # ---- Step 2: Select a genome ----
    genome: str = "hg38"

    # ---- Step 3: Select a PAM ----
    pam: str = "NGG"

    # ---- 运行参数（不影响三个 Step 的语义） ----
    show_all_scores: bool = False   # True 时额外下载 Chari/Xu/Wang/Doench'14 等分
    poll_interval: float = 5.0      # 轮询间隔（秒）
    timeout: float = 900.0          # 单次查询最长等待（秒）
    retries: int = 3                # 网络失败重试次数
    batch_id: str | None = None     # 指定则跳过提交，直接复用这个 batchId
    use_cache: bool = True          # False = 完全不用磁盘缓存

    # ---- 第三个数据（距离）用 ----
    target_pos: int = 22_163_096    # 目标位点（hg38 chr8，SFTPC c.218 / I73T）
    genome_offset: int | None = None  # 输入序列的基因组 0-based 起点；不给则自动定位
    genome_chrom: str | None = None   # 仅显示用，如 chr8
    auto_locate: bool = True          # 只给碱基时，自动去基因组定位 offset
    # ---- 序列定位窗口（换基因时改这三个） ----
    locus_chrom: str | None = None    # 如 "chr8"；None = 用默认 SFTPC 基因座
    locus_start: int | None = None    # 0-based
    locus_end: int | None = None

    extra: dict = field(default_factory=dict)

    # ---------------------------------------------------------------

    def resolve_step1(self) -> str:
        """把 Step1 的三种写法归一成 CRISPOR 要的 seq 字符串。

        返回值要么是纯序列，要么是染色体区间字符串（CRISPOR 服务端自己会去取序列）。
        """
        if self.fasta_path:
            return _read_fasta(self.fasta_path)

        if self.region:
            region = self.region.strip()
            if not _REGION_RE.match(region):
                raise CrisporConfigError(
                    f"Step1 的染色体区间格式不对：{region!r}；应形如 chr8:22,119,000-22,125,000"
                )
            return region.replace(" ", "")

        if self.sequence:
            seq = _clean_sequence(self.sequence)
            if not seq:
                raise CrisporConfigError("Step1 的序列为空")
            if len(seq) > MAX_SEQ_LEN:
                raise CrisporConfigError(
                    f"Step1 序列长度 {len(seq)} bp 超过 CRISPOR 上限 {MAX_SEQ_LEN} bp；"
                    "请只放一个外显子/目标区域，而不是整条基因"
                )
            return seq

        raise CrisporConfigError(
            "Step1 未指定输入：请给出 sequence / region / fasta_path 三者之一"
        )

    def validate(self) -> None:
        """三 step 全量校验，及早报错，不要等提交到服务器才失败。"""
        self.resolve_step1()
        if not self.genome or not self.genome.strip():
            raise CrisporConfigError("Step2 未指定基因组（genome），例如 hg38")
        if not self.pam or not self.pam.strip():
            raise CrisporConfigError("Step3 未指定 PAM（pam），例如 NGG")
        if self.pam.strip() != self.pam:
            raise CrisporConfigError(f"Step3 的 PAM 有多余空格：{self.pam!r}")

    def describe(self) -> str:
        """人类可读的三 step 摘要，提交前打印用。"""
        s1 = self.resolve_step1()
        if _REGION_RE.match(s1):
            s1_desc = f"染色体区间 {s1}"
        else:
            s1_desc = f"{len(s1)} bp 序列（{s1[:20]}{'...' if len(s1) > 20 else ''}）"
        return (
            "  Step1 Target sequence : " + s1_desc + f"   [name={self.name}]\n"
            "  Step2 Genome          : " + self.genome + "\n"
            "  Step3 PAM             : " + self.pam
        )

    def to_form(self) -> dict:
        """转成 POST 给 crispor.py 的表单字段。"""
        self.validate()
        return {
            "seq": self.resolve_step1(),
            "org": self.genome,
            "pam": self.pam,
            "name": self.name,
            "submit": "SUBMIT",
        }


# ------------------------------------------------------------------ helpers

def _read_fasta(path: str) -> str:
    if not os.path.isfile(path):
        raise CrisporConfigError(f"FASTA 文件不存在：{path}")
    chunks: list[str] = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            chunks.append(line)
    seq = _clean_sequence("".join(chunks))
    if not seq:
        raise CrisporConfigError(f"FASTA 文件里没有读到序列：{path}")
    if len(seq) > MAX_SEQ_LEN:
        raise CrisporConfigError(
            f"FASTA 序列长度 {len(seq)} bp 超过 CRISPOR 上限 {MAX_SEQ_LEN} bp"
        )
    return seq


def _clean_sequence(text: str) -> str:
    """去掉空白和换行，保留 ACGTN。"""
    seq = re.sub(r"\s+", "", text or "")
    if not seq:
        return ""
    if not _ALLOWED_BASES.match(seq):
        bad = sorted(set(re.findall(r"[^ACGTNacgtn]", seq)))
        raise CrisporConfigError(f"序列里出现非 ACGTN 字符：{bad}")
    return seq.upper()
