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
"""CHOPCHOP 的输入定义。

对应 CHOPCHOP 首页那四个框：

    Target            -> target      （基因名 / 染色体区间，或 sequence 直接给序列）
    In  (genome)      -> genome      （hg38 / hg19 ...）
    Using (nuclease)  -> nuclease    （CRISPR / CAS13 / TALEN / NICKASE / CPF1）
    For  (application)-> application （knockout / knockin / ...）

外加高级选项（Options 页）：target_region / pam / guide_len / gc 过滤 / scoring。
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field

# CHOPCHOP 的 nuclease 取值（对应 #typeSelect）
NUCLEASES = {
    "CRISPR": "CRISPR/Cas9（默认）",
    "CAS13": "CRISPR/Cas13（knock-down）",
    "CPF1": "CRISPR/Cpf1",
    "NICKASE": "Cas9 Nickase",
    "TALEN": "TALEN",
    "CAS12": "Cas12",
}

# For 下拉（对应 #forSelect）
APPLICATIONS = {
    "knock-out": "基因敲除（官网默认）",
    "knock-in": "敲入 / HDR",
    "knock-down": "knock-down（Cas13）",
    "activation": "激活 CRISPRa",
    "repression": "抑制 CRISPRi",
    "nanoporE": "Nanopore 富集",
}

# 靶向区域（对应 #targetRadio）
TARGET_REGIONS = {
    "WHOLE": "整条基因（含 UTR）",
    "CODING": "只看编码区 CDS",
    "SPLICE": "剪接位点",
    "UTR5": "5' UTR",
    "UTR3": "3' UTR",
    "PROMOTER": "启动子",
}

# 效率评分模型（对应 #scoringMatrix，传给 -scoringMethod）
SCORING_METHODS = {
    "DOENCH_2016": "Doench 2016（默认）",
    "DOENCH_2014": "Doench 2014",
    "CHARI_2015": "Chari 2015",
    "XU_2015": "Xu 2015",
    "MORENO_MATEOS_2015": "Moreno-Mateos 2015",
    "G20": "G20（只看 PAM 前一位是不是 G）",
}

COMMON_GENOMES = {
    "hg38": "Homo sapiens GRCh38/hg38",
    "hg19": "Homo sapiens GRCh37/hg19",
    "mm10": "Mus musculus mm10",
    "mm39": "Mus musculus mm39",
    "danRer11": "Danio rerio GRCz11",
    "dm6": "Drosophila melanogaster dm6",
    "ce11": "C. elegans ce11",
    "SacCer3": "S. cerevisiae SacCer3",
}

_SEQ_RE = re.compile(r"^[ACGTNacgtn\s]+$")
_REGION_RE = re.compile(r"^\s*[A-Za-z0-9_.\-]+\s*:\s*[\d,]+\s*-\s*[\d,]+\s*$")


class ChopchopConfigError(ValueError):
    """输入不合法。"""


@dataclass
class ChopchopInput:
    """一次 CHOPCHOP 查询的输入（= 首页四个框 + Options）。"""

    # ---- Target（二选一）----
    target: str | None = None        # 基因名 "SFTPC"，或区间 "chr8:22157000-22165000"
    sequence: str | None = None      # 或直接给序列（走 fastaInput）
    sequence_name: str = "query"

    # ---- In / Using / For ----
    genome: str = "hg38"
    nuclease: str = "CRISPR"
    application: str = "knock-out"

    # ---- Options ----
    target_region: str = "CODING"   # 官网默认选中 CODING
    pam: str = "NGG"                 # CHOPCHOP 默认 NGG，一般不用改
    guide_len: int = 20
    scoring: str = "DOENCH_2016"
    gc_min: int = 10
    gc_max: int = 90
    five_prime: str = "NN"
    enzyme: str = "BspQI"            # -n，TALEN 的酶切位点用，CRISPR 下无影响
    consensus_union: bool = False
    exclude_snv: bool = False        # -rm1perfOff：排除在 SNP 上有 1 个错配的 guides

    # ---- 序列模式的坐标换算 ----
    # 用 sequence（fastaInput）提交时，CHOPCHOP 给的是 "seq:N"（**序列内坐标**），
    # 需要 genome_offset = 这段序列在基因组上的 0-based 起始坐标，
    # 才能换算：基因组坐标 = genome_offset + N
    genome_offset: int | None = None
    genome_chrom: str | None = None   # 仅用于输出显示，如 "chr8"
    auto_locate: bool = True          # 用序列提交且没给 offset 时，自动去基因组里定位
    # ---- 序列定位窗口（换基因时改这三个） ----
    locus_chrom: str | None = None    # 如 "chr8"；None = 用默认 SFTPC 基因座
    locus_start: int | None = None    # 0-based
    locus_end: int | None = None


    # ---- 运行参数 ----
    poll_interval: float = 3.0
    timeout: float = 600.0
    retries: int = 3
    use_cache: bool = True
    job_id: str | None = None        # 指定则跳过提交，直接复用这个 jobId
    extra: dict = field(default_factory=dict)

    # ---------------------------------------------------------------

    def split_target(self) -> tuple[str, str]:
        """返回 (geneInput, fastaInput)，对应 API 的两个字段。"""
        if self.sequence:
            seq = re.sub(r"\s+", "", self.sequence).upper()
            if not seq:
                raise ChopchopConfigError("sequence 为空")
            bad = sorted(set(re.findall(r"[^ACGTN]", seq)))
            if bad:
                raise ChopchopConfigError(f"序列里出现非 ACGTN 字符：{bad}")
            return "", seq
        if self.target:
            return self.target.strip(), ""
        raise ChopchopConfigError("target 与 sequence 必须给一个")

    def validate(self) -> None:
        gene, fasta = self.split_target()
        if not gene and not fasta:
            raise ChopchopConfigError("Target 为空")
        if self.genome not in COMMON_GENOMES and not re.match(r"^[A-Za-z0-9_.]+$", self.genome):
            raise ChopchopConfigError(f"基因组名不合法：{self.genome!r}")
        if self.genome_offset is not None and self.genome_offset < 0:
            raise ChopchopConfigError("genome_offset 必须是 >=0 的整数（0-based 起点）")
        if self.nuclease not in NUCLEASES:
            raise ChopchopConfigError(
                f"nuclease 必须是 {sorted(NUCLEASES)} 之一，收到 {self.nuclease!r}"
            )
        if self.target_region not in TARGET_REGIONS:
            raise ChopchopConfigError(
                f"target_region 必须是 {sorted(TARGET_REGIONS)} 之一"
            )
        if self.scoring not in SCORING_METHODS:
            raise ChopchopConfigError(f"scoring 必须是 {sorted(SCORING_METHODS)} 之一")
        if not (0 <= self.gc_min < self.gc_max <= 100):
            raise ChopchopConfigError("gc_min/gc_max 不合法")

    def build_opts(self) -> list[str]:
        """拼 CHOPCHOP 的 opts 参数数组（和网页 formController.js 一样）。"""
        opts: list[str] = ["-J", "-BED", "-GenBank", "-G", self.genome,
                           "-filterGCmin", str(self.gc_min),
                           "-filterGCmax", str(self.gc_max)]
        if self.consensus_union:
            opts.append("-consensusUnion")
        if self.exclude_snv:
            opts.append("-rm1perfOff")
        opts += ["-t", self.target_region]
        opts += ["-n", self.enzyme, "-R", "4"]
        if self.nuclease == "CRISPR":
            opts += ["-T", "1", "-g", str(self.guide_len),
                     "-scoringMethod", self.scoring,
                     "-f", self.five_prime, "-v", "3", "-w"]
        return opts

    def to_form(self) -> dict:
        self.validate()
        gene, fasta = self.split_target()
        return {
            "opts": self.build_opts(),
            "fastaInput": fasta,
            "geneInput": gene,
            "isIsoform": self.nuclease == "CAS13",
            "forSelect": self.application,
        }

    def describe(self) -> str:
        gene, fasta = self.split_target()
        target = f"{len(fasta)} bp 序列" if fasta else f"{gene!r}"
        return (
            f"  Target  : {target}   [name={self.sequence_name}]\n"
            f"  Genome  : {self.genome}\n"
            f"  Nuclease: {self.nuclease}\n"
            f"  For     : {self.application}   靶向区域={self.target_region}  评分={self.scoring}"
        )
