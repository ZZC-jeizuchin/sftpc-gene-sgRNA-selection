"""兼容层：实现已经搬到工作区根目录的 genome_utils（两个爬虫共用）。"""

from genome_utils import (  # noqa: F401
    DEFAULT_LOCUS, fetch_locus, locate_sequence,
)

# 旧名字
locate = locate_sequence

__all__ = ["locate", "locate_sequence", "fetch_locus", "DEFAULT_LOCUS"]
