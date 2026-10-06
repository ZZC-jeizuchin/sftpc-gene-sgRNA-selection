"""JSON 输入解析 + JSON 输出构造。

输入 JSON（一次查询）：
{
  "step1": {"seq": "ACGT..."}        // 或 {"region": "chr8:..."} 或 {"fasta": "x.fa"}
  "step2": {"genome": "hg38"},
  "step3": {"pam": "NGG"},
  "name": "sftpc",
  "data": ["doench16", "cfd_spec"],  // 要哪几个数据，缺省=全部
  "batch_id": "E7soURXwF4IZYUeAKOwP",// 可选，指定则跳过提交
  "options": {"all_scores": false, "cache": true,
              "timeout": 900, "poll_interval": 3, "retries": 3}
}

也接受：直接一个数组（批量），或 {"queries": [ ... ]}。

输出 JSON（每个数据独立成块，互不干扰）：
{
  "ok": true,
  "query":   { "step1": {...}, "step2": {...}, "step3": {...}, "name": "..." },
  "crispor": { "batch_id": "...", "url": "...", "n_guides": 52, "from_cache": false },
  "datasets": {
      "doench16": { ...元信息..., "count": 52, "stats": {...}, "values": [ {...}, ... ] },
      "cfd_spec": { ... }
  }
}
"""

from __future__ import annotations

import json
from pathlib import Path

from .config import CrisporConfigError, CrisporInput
from .datasets import DatasetSpec, build_payload, resolve
from .pipeline import CrisporResult


# ============================================================ 输入：JSON -> CrisporInput

def load_queries(source: str | dict | list) -> tuple[list[CrisporInput], list[str]]:
    """解析输入 JSON。

    返回 (查询列表, 数据 id 列表)。数据 id 取第一条里写的，缺省=全部。
    """
    obj = _as_object(source)

    if isinstance(obj, list):
        raw_queries = obj
        common: dict = {}
    elif isinstance(obj, dict) and "queries" in obj:
        common = {k: v for k, v in obj.items() if k != "queries"}
        raw_queries = obj["queries"]
    elif isinstance(obj, dict):
        raw_queries = [obj]
        common = {}
    else:
        raise CrisporConfigError("输入 JSON 必须是对象、对象数组，或 {\"queries\": [...]}")

    if not raw_queries:
        raise CrisporConfigError("输入 JSON 里没有任何查询")

    inputs: list[CrisporInput] = []
    data_ids: list[str] = []
    for i, raw in enumerate(raw_queries):
        if not isinstance(raw, dict):
            raise CrisporConfigError(f"第 {i + 1} 个查询不是 JSON 对象")
        merged = {**common, **raw}
        inp, ids = parse_query(merged)
        if ids and not data_ids:
            data_ids = ids
        inputs.append(inp)
    return inputs, data_ids


def parse_query(obj: dict) -> tuple[CrisporInput, list[str]]:
    """单个查询对象 -> CrisporInput。"""
    step1 = obj.get("step1")
    if isinstance(step1, str):          # 容错：直接把 step1 写成字符串
        step1 = {"seq": step1}
    if not isinstance(step1, dict):
        raise CrisporConfigError(
            "缺少 step1；应为 {\"seq\": ...} / {\"region\": ...} / {\"fasta\": ...}"
        )

    step2 = obj.get("step2") or {}
    if isinstance(step2, str):
        step2 = {"genome": step2}
    step3 = obj.get("step3") or {}
    if isinstance(step3, str):
        step3 = {"pam": step3}

    options = obj.get("options") or {}
    if not isinstance(options, dict):
        raise CrisporConfigError("options 必须是 JSON 对象")

    given = sum(1 for k in ("seq", "region", "fasta") if step1.get(k))
    if given == 0:
        raise CrisporConfigError("step1 里必须给 seq / region / fasta 之一")
    if given > 1:
        raise CrisporConfigError("step1 里 seq / region / fasta 只能给一个")

    inp = CrisporInput(
        sequence=step1.get("seq"),
        region=step1.get("region"),
        fasta_path=step1.get("fasta"),
        name=obj.get("name", "query"),
        genome=step2.get("genome", "hg38"),
        pam=step3.get("pam", "NGG"),
        show_all_scores=bool(options.get("all_scores", False)),
        poll_interval=float(options.get("poll_interval", 3.0)),
        timeout=float(options.get("timeout", 900.0)),
        retries=int(options.get("retries", 3)),
        batch_id=obj.get("batch_id"),
        use_cache=bool(options.get("cache", True)),
        target_pos=int(obj.get("target_pos", options.get("target_pos", 22_163_096))),
        genome_offset=(int(obj["genome_offset"]) if obj.get("genome_offset") is not None else None),
        genome_chrom=obj.get("genome_chrom"),
        auto_locate=bool(options.get("auto_locate", obj.get("auto_locate", True))),
    )
    data_ids = obj.get("data") or []
    if not isinstance(data_ids, list):
        raise CrisporConfigError("data 必须是数组，如 [\"doench16\", \"cfd_spec\"]")
    resolve(data_ids)          # 提前校验 id 合法性
    return inp, data_ids


def _as_object(source: str | dict | list):
    if isinstance(source, (dict, list)):
        return source
    text = source
    if isinstance(source, (str, Path)) and Path(str(source)).is_file():
        text = Path(str(source)).read_text(encoding="utf-8")
    try:
        return json.loads(text)
    except json.JSONDecodeError as exc:
        raise CrisporConfigError(f"输入不是合法 JSON：{exc}") from exc


# ============================================================ 输出：CrisporResult -> JSON

def query_to_json(inp: CrisporInput) -> dict:
    """把三个 Step 原样回显（已经过校验/归一化）。"""
    step1 = inp.resolve_step1()
    if inp.fasta_path:
        s1 = {"fasta": inp.fasta_path, "seq": step1, "length": len(step1)}
    elif inp.region:
        s1 = {"region": step1}
    else:
        s1 = {"seq": step1, "length": len(step1)}
    return {
        "step1": s1,
        "step2": {"genome": inp.genome},
        "step3": {"pam": inp.pam},
        "name": inp.name,
        "target_pos": inp.target_pos,
        "genome_offset": inp.genome_offset,
        "genome_chrom": inp.genome_chrom,
    }


def result_to_json(
    result: CrisporResult,
    specs: list[DatasetSpec] | None = None,
) -> dict:
    """把结果按数据拆成 JSON。datasets 里每个 id 一个独立块。"""
    specs = specs or resolve(None)
    datasets = {spec.id: build_payload(spec, result.records) for spec in specs}
    # 只给碱基序列时，附上第三个数据（距离）
    if result.input.sequence or result.input.genome_offset is not None:
        datasets["distance_to_22163096"] = result.distance_payload()
    return {
        "ok": True,
        "query": query_to_json(result.input),
        "crispor": {
            "batch_id": result.batch_id,
            "url": result.url,
            "n_guides": len(result.records),
            "from_cache": result.from_cache,
            "reused_batch": result.reused_batch,
        },
        "datasets": datasets,
    }


def error_to_json(message: str, query: dict | None = None, kind: str = "error") -> dict:
    return {"ok": False, "error": {"kind": kind, "message": message}, "query": query}


def dumps(obj) -> str:
    return json.dumps(obj, ensure_ascii=False, indent=2)


def write_json(path: str | Path, obj) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(dumps(obj), encoding="utf-8")
    return path
