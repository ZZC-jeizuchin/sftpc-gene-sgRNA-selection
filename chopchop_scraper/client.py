"""CHOPCHOP HTTP 客户端。

实测到的接口：
    POST /                                  JSON -> {"jobId": "..."}
    GET  /results/<jobId>/run.info          tab 分隔
    GET  /results/<jobId>/query.json        回显请求
    GET  /results/<jobId>/results.tsv       结果表（核心）
    GET  /results/<jobId>/cutcoords.json    绘图坐标
    GET  /results/<jobId>/results.bed

轮询策略与网站1 一致：**不去刷结果 HTML 页**，而是直接重试那个小的 results.tsv。
"""

from __future__ import annotations

import json
import time
import urllib.error
import urllib.parse
import urllib.request

from .cache import ChopchopCache
from .config import ChopchopInput

BASE_URL = "https://chopchop.cbu.uib.no/"
RESULTS_URL = "https://chopchop.cbu.uib.no/results/{job_id}/{name}"

# 结果表第一行的列名（用它判断"算好了没"）
_TSV_HEADER = "Rank\tTarget sequence\tGenomic location"


class ChopchopError(RuntimeError):
    """通信或服务端报错。"""


class ChopchopClient:
    def __init__(
        self,
        base_url: str = BASE_URL,
        timeout: float = 60.0,
        retries: int = 3,
        user_agent: str = "Mozilla/5.0 (chopchop-scraper/1.0; educational project)",
        cache: ChopchopCache | None = None,
    ) -> None:
        self.base_url = base_url
        self.timeout = timeout
        self.retries = max(1, int(retries))
        self.user_agent = user_agent
        self.cache = cache if cache is not None else ChopchopCache()

    # ------------------------------------------------------------- HTTP

    def _request(self, url, data=None, use_cache=True, retries=None, timeout=None,
                 method=None):
        if use_cache:
            hit = self.cache.get_raw(method or ("POST" if data is not None else "GET"), url, data)
            if hit is not None:
                return hit

        tries = self.retries if retries is None else max(1, int(retries))
        tmo = self.timeout if timeout is None else timeout
        last_err = None
        for attempt in range(1, tries + 1):
            req = urllib.request.Request(
                url, data=data, method=method,
                headers={
                    "User-Agent": self.user_agent,
                    "Accept": "application/json,text/plain,text/html,*/*",
                    "Accept-Language": "en-US,en;q=0.9",
                },
            )
            if data is not None:
                req.add_header("Content-Type", "application/json")
            try:
                with urllib.request.urlopen(req, timeout=tmo) as resp:
                    text = resp.read().decode("utf-8", errors="replace")
                if use_cache:
                    self.cache.put_raw(
                        method or ("POST" if data is not None else "GET"), url, data, text
                    )
                return text
            except urllib.error.HTTPError as exc:
                # 404 = 结果还没生成，属于正常轮询状态
                if exc.code == 404:
                    raise ChopchopError(f"404 Not Found（通常表示结果还没生成）：{url}") from exc
                body = ""
                try:
                    body = exc.read().decode("utf-8", errors="replace")[:300]
                except Exception:
                    pass
                last_err = f"HTTP {exc.code}: {body}"
            except (urllib.error.URLError, TimeoutError, OSError) as exc:
                last_err = exc
            if attempt < tries:
                time.sleep(min(2 ** attempt, 10))
        raise ChopchopError(f"请求失败（重试 {tries} 次）：{url} :: {last_err}")

    # ------------------------------------------------------------- 提交

    def submit(self, inp: ChopchopInput) -> str:
        form = inp.to_form()
        payload = json.dumps(form, ensure_ascii=False).encode("utf-8")
        text = self._request(self.base_url, data=payload, use_cache=False, method="POST")
        try:
            obj = json.loads(text)
        except json.JSONDecodeError as exc:
            raise ChopchopError(f"提交返回的不是 JSON：{text[:300]}") from exc
        job_id = obj.get("jobId") or obj.get("jobID")
        if not job_id:
            raise ChopchopError(f"提交后没拿到 jobId：{text[:300]}")
        return job_id

    def submit_cached(self, inp: ChopchopInput, on_progress=None) -> tuple[str, bool]:
        key = ChopchopCache.key_from_input(inp)
        cached = self.cache.get_job_id(key)
        if cached:
            if on_progress:
                on_progress(f"缓存命中：复用旧 jobId {cached}（跳过提交）")
            return cached, True
        job_id = self.submit(inp)
        self.cache.put_job_id(key, job_id, describe={"input": inp.describe()})
        return job_id, False

    # ------------------------------------------------------------- 结果

    def result_url(self, job_id: str, name: str) -> str:
        return RESULTS_URL.format(job_id=job_id, name=name)

    def fetch_run_info(self, job_id: str) -> str:
        return self._request(self.result_url(job_id, "run.info"), retries=1, timeout=45)

    def fetch_query(self, job_id: str) -> dict:
        text = self._request(self.result_url(job_id, "query.json"), retries=1, timeout=45)
        try:
            return json.loads(text)
        except json.JSONDecodeError:
            return {}

    def fetch_cutcoords(self, job_id: str) -> list:
        text = self._request(self.result_url(job_id, "cutcoords.json"), retries=1, timeout=45)
        try:
            return json.loads(text)
        except json.JSONDecodeError:
            return []

    def download_results_tsv(self, job_id: str, fast=False, use_cache=True) -> str:
        text = self._request(
            self.result_url(job_id, "results.tsv"),
            use_cache=use_cache and not fast,
            retries=1 if fast else None,
            timeout=45 if fast else None,
        )
        if _TSV_HEADER not in text:
            raise ChopchopError(
                "拿到的不是 results.tsv（可能还没算完）：" + text[:200].replace("\n", " ")
            )
        return text

    def fetch_results_when_ready(
        self,
        job_id: str,
        poll_interval: float = 3.0,
        timeout: float = 600.0,
        on_progress=None,
    ) -> str:
        """直接重试 results.tsv（小文件），不刷结果页。"""
        started = time.time()
        attempt = 0
        last_error = None
        while True:
            attempt += 1
            try:
                tsv = self.download_results_tsv(job_id, fast=True)
                if on_progress:
                    on_progress(f"第 {attempt} 次探测：算好了（{time.time() - started:.0f}s）")
                return tsv
            except ChopchopError as exc:
                last_error = exc

            elapsed = time.time() - started
            if elapsed > timeout:
                raise ChopchopError(
                    f"等待超时（{timeout:.0f}s），jobId={job_id} 仍未算完；最后探测：{last_error}"
                )
            if on_progress:
                on_progress(f"第 {attempt} 次探测：还在算… 已等待 {elapsed:.0f}s")
            time.sleep(poll_interval)
