"""CHOPCHOP 的 JSON 输入 / 输出。

输入（一次查询）：
{
  "target": "SFTPC",              // 或 "chr8:22157000-22165000"，或改用 "sequence": "ACGT..."
  "genome": "hg38",
  "nuclease": "CRISPR",
  "application": "knockout",
  "target_region": "WHOLE",
  "sequence_name": "sftpc",
  "target_pos": 22163096,          // 第三个数据的目标位点，缺省 22163096
  "job_id": "a37d4ff6-...",        // 可选，跳过提交
  "options": {"scoring": "DOENCH_2016", "gc_min": 10, "gc_max": 90,
              "cache": true, "timeout": 600, "poll_interval": 3}
}

输出（第三个数据独立成块）：
{
  "ok": true,
  "query":   {...},
  "chopchop": {"job_id": ..., "url": ..., "n_guides": 161, "run_info": {...}},
  "datasets": {"distance_to_22163096": {...}}
}
"""

from __future__ import annotations

import json
from pathlib import Path

from .config import ChopchopConfigError, ChopchopInput
from .distance import DEFAULT_TARGET_POS
from .pipeline import ChopchopResult


def load_queries(source) -> list[tuple[ChopchopInput, int]]:
    """返回 [(输入, target_pos), ...]。"""
    obj = _as_object(source)
    if isinstance(obj, list):
        raws, common = obj, {}
    elif isinstance(obj, dict) and "queries" in obj:
        common = {k: v for k, v in obj.items() if k != "queries"}
        raws = obj["queries"]
    elif isinstance(obj, dict):
        raws, common = [obj], {}
    else:
        raise ChopchopConfigError("输入 JSON 必须是对象、数组，或 {\"queries\": [...]}")
    if not raws:
        raise ChopchopConfigError("输入 JSON 里没有查询")

    out = []
    for i, raw in enumerate(raws):
        if not isinstance(raw, dict):
            raise ChopchopConfigError(f"第 {i + 1} 个查询不是 JSON 对象")
        out.append(parse_query({**common, **raw}))
    return out


def parse_query(obj: dict) -> tuple[ChopchopInput, int]:
    options = obj.get("options") or {}
    if not isinstance(options, dict):
        raise ChopchopConfigError("options 必须是 JSON 对象")
    target_pos = obj.get("target_pos", options.get("target_pos", DEFAULT_TARGET_POS))
    try:
        target_pos = int(target_pos)
    except (TypeError, ValueError):
        raise ChopchopConfigError(f"target_pos 必须是整数，收到 {target_pos!r}")

    inp = ChopchopInput(
        target=obj.get("target"),
        sequence=obj.get("sequence"),
        sequence_name=obj.get("sequence_name", obj.get("name", "query")),
        genome=obj.get("genome", "hg38"),
        nuclease=obj.get("nuclease", "CRISPR"),
        application=obj.get("application", "knock-out"),
        target_region=obj.get("target_region", "CODING"),
        pam=obj.get("pam", "NGG"),
        guide_len=int(options.get("guide_len", 20)),
        scoring=options.get("scoring", "DOENCH_2016"),
        gc_min=int(options.get("gc_min", 10)),
        gc_max=int(options.get("gc_max", 90)),
        five_prime=options.get("five_prime", "NN"),
        consensus_union=bool(options.get("consensus_union", False)),
        exclude_snv=bool(options.get("exclude_snv", False)),
        poll_interval=float(options.get("poll_interval", 3.0)),
        timeout=float(options.get("timeout", 600.0)),
        retries=int(options.get("retries", 3)),
        use_cache=bool(options.get("cache", True)),
        job_id=obj.get("job_id"),
        genome_offset=(int(obj["genome_offset"]) if obj.get("genome_offset") is not None else None),
        genome_chrom=obj.get("genome_chrom"),
        auto_locate=bool(options.get("auto_locate", obj.get("auto_locate", True))),
    )
    inp.validate()
    return inp, target_pos


def _as_object(source):
    if isinstance(source, (dict, list)):
        return source
    text = str(source)
    if Path(text).is_file():
        text = Path(text).read_text(encoding="utf-8")
    try:
        return json.loads(text)
    except json.JSONDecodeError as exc:
        raise ChopchopConfigError(f"输入不是合法 JSON：{exc}") from exc


# ------------------------------------------------------------------ 输出

def query_to_json(inp: ChopchopInput, target_pos: int) -> dict:
    gene, fasta = inp.split_target()
    return {
        "target": gene or None,
        "sequence": (fasta[:40] + "...") if len(fasta) > 40 else (fasta or None),
        "sequence_length": len(fasta) if fasta else None,
        "genome": inp.genome,
        "nuclease": inp.nuclease,
        "application": inp.application,
        "target_region": inp.target_region,
        "sequence_name": inp.sequence_name,
        "target_pos": target_pos,
        "genome_offset": inp.genome_offset,
        "auto_locate": inp.auto_locate,
        "coordinate_mode": "sequence" if fasta else "genomic",
    }


def result_to_json(result: ChopchopResult, target_pos: int = DEFAULT_TARGET_POS) -> dict:
    return {
        "ok": True,
        "query": query_to_json(result.input, target_pos),
        "chopchop": {
            "job_id": result.job_id,
            "url": result.url,
            "n_guides": len(result.records),
            "run_info": result.run_info,
            "from_cache": result.from_cache,
        },
        "datasets": {
            "distance_to_22163096": result.distance_payload(target_pos),
        },
    }


def error_to_json(message: str, query: dict | None = None, kind: str = "error") -> dict:
    return {"ok": False, "error": {"kind": kind, "message": message}, "query": query}


def dumps(obj) -> str:
    return json.dumps(obj, ensure_ascii=False, indent=2)


def write_json(path, obj) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(dumps(obj), encoding="utf-8")
    return path
