"""CHOPCHOP 高层流水线。"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

from .cache import ChopchopCache, DEFAULT_CACHE_DIR
from .client import BASE_URL, ChopchopClient
from .config import ChopchopInput
from .distance import DEFAULT_TARGET_POS, build_distance_payload
from .parser import GuideRecord, parse_results_tsv


@dataclass
class ChopchopResult:
    job_id: str
    input: ChopchopInput
    records: list[GuideRecord]
    raw_tsv: str = ""
    run_info: dict = field(default_factory=dict)
    from_cache: bool = False

    @property
    def url(self) -> str:
        return f"https://chopchop.cbu.uib.no/results/{self.job_id}/"

    def distance_payload(self, target_pos: int = DEFAULT_TARGET_POS) -> dict:
        """第三个数据（按数据块输出）。

        序列提交时（location 为 "seq:N"）会用 input.genome_offset 换算成基因组坐标；
        没给 offset 时 distance 为 null，并在 needs_genome_offset 里标出来。
        """
        payload = build_distance_payload(self.records, target_pos,
                                         offset=self.input.genome_offset)
        if self.input.genome_chrom:
            payload["genome_chrom"] = self.input.genome_chrom
            for v in payload["values"]:
                if v["genomic_location"] is not None:
                    v["genomic_location_str"] = f"{self.input.genome_chrom}:{v['genomic_location']}"
        return payload

    def save(self, out_dir: str | Path = ".") -> dict:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        base = f"chopchop_{self.input.sequence_name}_{stamp}"
        raw = out_dir / f"{base}_results.tsv"
        raw.write_text(self.raw_tsv, encoding="utf-8")
        return {"results_tsv": str(raw)}


def fetch(
    inp: ChopchopInput,
    client: ChopchopClient | None = None,
    on_progress=None,
    cache_dir: str | Path = DEFAULT_CACHE_DIR,
) -> ChopchopResult:
    """提交 -> 轮询 results.tsv -> 解析。走缓存。"""
    def log(msg: str) -> None:
        if on_progress:
            on_progress(msg)

    inp.validate()
    cache = ChopchopCache(cache_dir, enabled=inp.use_cache)

    # 0. 序列提交但没给 offset：自动定位回基因组（这样"只给碱基"也能算距离）
    locate_info = None
    if inp.sequence and inp.genome_offset is None and inp.auto_locate:
        from genome_utils import DEFAULT_LOCUS, DEFAULT_LOCUS_CACHE, locate_sequence
        loc_kw = {}
        if inp.locus_chrom:
            loc_kw = {"chrom": inp.locus_chrom,
                      "start": inp.locus_start if inp.locus_start is not None else DEFAULT_LOCUS["start"],
                      "end": inp.locus_end if inp.locus_end is not None else DEFAULT_LOCUS["end"]}
        locate_info = locate_sequence(inp.sequence, cache_dir=DEFAULT_LOCUS_CACHE, **loc_kw)
        if locate_info.get("found"):
            inp.genome_offset = locate_info["offset"]
            log(f"自动定位：序列在 {locate_info['chrom']} 上，genome_offset="
                f"{locate_info['offset']}（{locate_info['note']}）")
        else:
            log(f"⚠ 自动定位失败：{locate_info.get('note')}；基因组距离将算不出来")
    client = client or ChopchopClient(timeout=60, retries=inp.retries, cache=cache)
    client.cache = cache
    key = ChopchopCache.key_from_input(inp)

    # 1. 结果已在缓存
    if inp.use_cache:
        cached = cache.get_tsv(key)
        if cached:
            job_id = cache.get_job_id(key) or "cached"
            log(f"缓存命中：直接读本地结果（jobId={job_id}），0 次网络请求")
            records = parse_results_tsv(cached)
            log(f"解析到 {len(records)} 条候选 guide")
            return ChopchopResult(job_id, inp, records, cached,
                                  cache.get_meta(key), from_cache=True)

    # 2. 拿 jobId（提交或复用）
    job_id = inp.job_id or (cache.get_job_id(key) if inp.use_cache else None)
    if job_id:
        log(f"复用 jobId = {job_id}（跳过提交）")
        cache.put_job_id(key, job_id, describe={"input": inp.describe()})
    else:
        log("输入：\n" + inp.describe())
        log("提交任务…")
        job_id, _ = client.submit_cached(inp, on_progress=log)
        log(f"已排队，jobId = {job_id}")

    # 3. 直接重试 results.tsv
    log("等待计算完成（直接探测 results.tsv，不刷结果页）…")
    tsv = client.fetch_results_when_ready(
        job_id, poll_interval=inp.poll_interval, timeout=inp.timeout, on_progress=log
    )

    # 4. 顺手拿 run.info（算坐标偏移用得上：targets / genome / mode / ...）
    meta = {}
    try:
        info = client.fetch_run_info(job_id).strip().split("\t")
        meta = {"run_info_raw": info}
        if len(info) >= 5:
            meta.update({
                "targets": info[0], "genome": info[1], "mode": info[2],
                "unique_method": info[3], "guide_size": info[4],
            })
    except Exception:
        pass

    if locate_info:
        meta["locate"] = locate_info
    cache.put_tsv(key, job_id, tsv, describe={"input": inp.describe()}, meta=meta)
    records = parse_results_tsv(tsv)
    log(f"解析到 {len(records)} 条候选 guide，已写入缓存")
    return ChopchopResult(job_id, inp, records, tsv, meta, from_cache=False)
