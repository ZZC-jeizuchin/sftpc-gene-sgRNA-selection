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
"""与 CRISPOR CGI 接口通信的客户端。

CRISPOR 没有官方 REST API，但它网页前端用的就是一套 CGI 参数接口，
抓这套接口比抓渲染后的 HTML 稳定得多（列顺序、样式变了都不影响）。

实测到的三个端点：

1) 提交任务
   POST https://crispor.gi.ucsc.edu/crispor.py
        seq=<序列或染色体区间>&org=hg38&pam=NGG&name=<名字>&submit=SUBMIT
   -> 返回的 HTML 里含  history.replaceState('crispor.py', document.title, '?batchId=XXXX')

2) 轮询进度
   GET  https://crispor.gi.ucsc.edu/crispor.py?batchId=XXXX
   -> 未完成：页面含 "Status" / "refresh"，无结果表
   -> 已完成：页面含 "MIT Specificity" 与 download 链接

3) 下载结果（关键，直接是 TSV）
   GET  https://crispor.gi.ucsc.edu/crispor.py?batchId=XXXX&download=guides&format=tsv
   GET  ...&showAllScores=1&download=guides&format=tsv
"""

from __future__ import annotations

import re
import time
import urllib.error
import urllib.parse
import urllib.request

from .cache import CrisporCache
from .config import CrisporInput

BASE_URL = "https://crispor.gi.ucsc.edu/crispor.py"

# 提交响应里的 batchId
_BATCH_ID_RE = re.compile(r"batchId=([A-Za-z0-9]+)")

# 结果页完成标志（下载链接 / 结果表头）
_DONE_MARKERS = ("download=guides", "MIT Specificity")

# 结果页报错标志
_ERROR_MARKERS = (
    "We are very sorry",
    "has to be shorter",
    "no PAM sites",
    "abort",
)

_STATUS_MARKERS = ("Status", "refresh")


class CrisporError(RuntimeError):
    """与 CRISPOR 通信或服务端报错。"""


class CrisporClient:
    """一次 CRISPOR 会话的 HTTP 客户端（无状态，可复用）。"""

    def __init__(
        self,
        base_url: str = BASE_URL,
        timeout: float = 60.0,
        retries: int = 3,
        user_agent: str = "crispor-scraper/1.0 (+educational project)",
        cache: CrisporCache | None = None,
    ) -> None:
        self.base_url = base_url
        self.timeout = timeout
        self.retries = max(1, int(retries))
        self.user_agent = user_agent
        self.cache = cache if cache is not None else CrisporCache()

    # ------------------------------------------------------------- 底层 HTTP

    def _request(
        self,
        url: str,
        data: bytes | None = None,
        use_cache: bool = True,
        retries: int | None = None,
        timeout: float | None = None,
    ) -> str:
        """发一次请求，带重试。返回解码后的文本。

        use_cache=True 时先查磁盘缓存；命中则 0 次网络请求。
        轮询必须 use_cache=False，否则永远读到旧的"排队中"页面。
        retries/timeout 可单独覆盖：轮询时用 retries=1 快速失败，不要卡住。
        """
        if use_cache:
            hit = self.cache.get_raw("POST" if data is not None else "GET", url, data)
            if hit is not None:
                return hit

        tries = self.retries if retries is None else max(1, int(retries))
        tmo = self.timeout if timeout is None else timeout
        last_err: Exception | None = None
        for attempt in range(1, tries + 1):
            req = urllib.request.Request(
                url,
                data=data,
                method="POST" if data is not None else "GET",
                headers={
                    "User-Agent": self.user_agent,
                    "Accept": "text/html,text/plain,*/*",
                    "Accept-Language": "en-US,en;q=0.9",
                },
            )
            if data is not None:
                req.add_header("Content-Type", "application/x-www-form-urlencoded")
            try:
                with urllib.request.urlopen(req, timeout=tmo) as resp:
                    raw = resp.read()
                text = raw.decode("utf-8", errors="replace")
                if use_cache:
                    self.cache.put_raw(
                        "POST" if data is not None else "GET", url, data, text
                    )
                return text
            except (urllib.error.URLError, TimeoutError, OSError) as exc:
                last_err = exc
                if attempt < tries:
                    time.sleep(min(2 ** attempt, 10))
        raise CrisporError(f"请求失败（已重试 {tries} 次）：{url} :: {last_err}")

    # ------------------------------------------------------------------ 1. 提交

    def submit(self, inp: CrisporInput) -> str:
        """提交三个 Step 的输入，返回 batchId。"""
        form = inp.to_form()
        payload = urllib.parse.urlencode(form).encode("utf-8")
        html = self._request(self.base_url, data=payload, use_cache=False)

        _raise_if_error(html, stage="提交")

        match = _BATCH_ID_RE.search(html)
        if not match:
            raise CrisporError(
                "提交成功但没在返回页面里找到 batchId；"
                f"页面片段：{_snippet(html)}"
            )
        return match.group(1)

    # ------------------------------------------------- 缓存感知的提交 / 下载

    def submit_cached(self, inp: CrisporInput, on_progress=None) -> tuple[str, bool]:
        """先查缓存，命中就复用旧 batchId，不再提交。

        返回 (batchId, 是否来自缓存)。
        """
        key = CrisporCache.key_from_input(inp)
        cached = self.cache.get_batch_id(key)
        if cached:
            if on_progress:
                on_progress(f"缓存命中：复用旧 batchId {cached}（跳过提交）")
            return cached, True
        batch_id = self.submit(inp)
        self.cache.put_batch_id(key, batch_id, describe={"input": inp.describe()})
        return batch_id, False

    def get_cached_tsv(self, inp: CrisporInput) -> str | None:
        """如果这条输入之前已经抓过结果，直接返回 TSV。"""
        return self.cache.get_tsv(CrisporCache.key_from_input(inp))

    def store_tsv(self, inp: CrisporInput, batch_id: str, tsv: str) -> None:
        self.cache.put_tsv(
            CrisporCache.key_from_input(inp), batch_id, tsv,
            describe={"input": inp.describe()},
        )

    # ------------------------------------------------------------------ 2. 轮询

    def fetch_status_page(self, batch_id: str) -> str:
        """整张结果页（约 265 KB，慢）。只在需要看报错信息时才拉。"""
        url = f"{self.base_url}?{urllib.parse.urlencode({'batchId': batch_id})}"
        return self._request(url, use_cache=False, retries=1, timeout=45)

    def fetch_guides_when_ready(
        self,
        batch_id: str,
        all_scores: bool = False,
        poll_interval: float = 3.0,
        timeout: float = 900.0,
        on_progress=None,
    ) -> str:
        """等任务算完并返回 guides TSV。

        关键优化：**不用结果页轮询**（265 KB / 反复刷新很慢），
        而是直接重试那个 6 KB 的 TSV 下载端点，算好了就一次拿到。
        """
        started = time.time()
        attempt = 0
        while True:
            attempt += 1
            try:
                tsv = self.download_guides(batch_id, all_scores=all_scores, fast=True)
                if on_progress:
                    on_progress(f"第 {attempt} 次探测：算好了（{time.time() - started:.0f}s）")
                return tsv
            except CrisporError as exc:
                last_error = exc

            elapsed = time.time() - started
            if elapsed > timeout:
                # 超时了才去拉一次重的结果页，把服务端真正的报错带出来
                try:
                    html = self.fetch_status_page(batch_id)
                    _raise_if_error(html, stage="计算")
                except CrisporError:
                    pass
                raise CrisporError(
                    f"等待超时（{timeout:.0f}s），batchId={batch_id} 仍未算完；"
                    f"最后一次探测：{last_error}"
                )
            if on_progress:
                on_progress(f"第 {attempt} 次探测：还在算… 已等待 {elapsed:.0f}s")
            time.sleep(poll_interval)

    # ------------------------------------------------------------------ 3. 下载

    def download_guides(
        self,
        batch_id: str,
        all_scores: bool = False,
        fast: bool = False,
        use_cache: bool = True,
    ) -> str:
        """下载 guide 汇总表（TSV 文本）。

        这个是 CRISPOR 的「Download as Excel tables -> Guides」同一份数据，
        列里直接含 cfdSpecScore 与 Doench '16-Score。

        fast=True 时只试一次、超时收紧，用于轮询探测。
        """
        params = {"batchId": batch_id, "download": "guides", "format": "tsv"}
        if all_scores:
            params["showAllScores"] = "1"
        url = f"{self.base_url}?{urllib.parse.urlencode(params)}"
        # 轮询探测时不要写缓存，否则会把"还没算好"的空结果缓存住
        text = self._request(
            url,
            use_cache=use_cache and not fast,
            retries=1 if fast else None,
            timeout=30 if fast else None,
        )
        if "guideId" not in text.split("\n", 1)[0]:
            raise CrisporError(
                "下载到的内容不像 guides TSV，可能任务还没算完或接口变了；"
                f"开头片段：{_snippet(text)}"
            )
        return text

    def download_offtargets(self, batch_id: str) -> str:
        """下载脱靶明细表（TSV）。需要每个位点的 CFD 切割分时用这个。"""
        params = {"batchId": batch_id, "download": "offtargets", "format": "tsv"}
        url = f"{self.base_url}?{urllib.parse.urlencode(params)}"
        return self._request(url)


# ------------------------------------------------------------------ helpers

def _raise_if_error(html: str, stage: str) -> None:
    lowered = html.lower()
    for marker in _ERROR_MARKERS:
        if marker.lower() in lowered:
            raise CrisporError(
                f"{stage}阶段服务端报错（命中标记 {marker!r}）：{_snippet(html)}"
            )


def _snippet(text: str, n: int = 400) -> str:
    text = re.sub(r"<[^>]+>", " ", text)
    text = re.sub(r"\s+", " ", text).strip()
    return text[:n]
