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
"""磁盘缓存：CHOPCHOP 也很慢，同样做两层。

1) 语义缓存 index.json：请求指纹 -> jobId + results.tsv 文件
2) 原始响应缓存 raw/：URL 的 sha1 -> 原始文本（results.tsv 等重复下载瞬间返回）

轮询请求不走原始响应缓存，否则永远读到旧结果。
"""

from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

DEFAULT_CACHE_DIR = ".chopchop_cache"


class ChopchopCache:
    def __init__(self, cache_dir: str | Path = DEFAULT_CACHE_DIR, enabled: bool = True):
        self.enabled = enabled
        self.root = Path(cache_dir)
        self.raw_dir = self.root / "raw"
        self.index_path = self.root / "index.json"
        self._index: dict | None = None

    # ------------------------------------------------------------ 索引

    def _load_index(self) -> dict:
        if self._index is None:
            if self.index_path.is_file():
                try:
                    self._index = json.loads(self.index_path.read_text(encoding="utf-8"))
                except (json.JSONDecodeError, OSError):
                    self._index = {}
            else:
                self._index = {}
        return self._index

    def _save_index(self) -> None:
        if not self.enabled:
            return
        self.root.mkdir(parents=True, exist_ok=True)
        self.index_path.write_text(
            json.dumps(self._load_index(), ensure_ascii=False, indent=2), encoding="utf-8"
        )

    # ------------------------------------------------------------ 指纹

    @staticmethod
    def request_key(form: dict) -> str:
        blob = json.dumps(form, sort_keys=True, ensure_ascii=False)
        return hashlib.sha1(blob.encode("utf-8")).hexdigest()[:16]

    @staticmethod
    def key_from_input(inp) -> str:
        # 名字不参与指纹，改名字不会重新请求
        form = inp.to_form()
        form.pop("sequence_name", None)
        return ChopchopCache.request_key(form)

    def _entry(self, key: str) -> dict:
        return self._load_index().get(key, {})

    # ------------------------------------------------------ 语义 get/put

    def get_job_id(self, key: str) -> str | None:
        return self._entry(key).get("job_id") if self.enabled else None

    def put_job_id(self, key: str, job_id: str, describe: dict | None = None) -> None:
        if not self.enabled:
            return
        idx = self._load_index()
        entry = idx.setdefault(key, {})
        entry["job_id"] = job_id
        entry["updated"] = time.strftime("%Y-%m-%d %H:%M:%S")
        if describe:
            entry["input"] = describe
        self._save_index()

    def get_tsv(self, key: str) -> str | None:
        if not self.enabled:
            return None
        fname = self._entry(key).get("tsv_file")
        if not fname:
            return None
        path = self.root / fname
        return path.read_text(encoding="utf-8", errors="replace") if path.is_file() else None

    def get_meta(self, key: str) -> dict:
        return self._entry(key).get("meta", {})

    def put_tsv(self, key: str, job_id: str, tsv: str,
                describe: dict | None = None, meta: dict | None = None) -> Path:
        self.root.mkdir(parents=True, exist_ok=True)
        fname = f"results_{key}.tsv"
        path = self.root / fname
        path.write_text(tsv, encoding="utf-8")
        if self.enabled:
            idx = self._load_index()
            entry = idx.setdefault(key, {})
            entry.update({
                "job_id": job_id,
                "tsv_file": fname,
                "saved": time.strftime("%Y-%m-%d %H:%M:%S"),
                "n_bytes": len(tsv),
            })
            if describe:
                entry["input"] = describe
            if meta:
                entry["meta"] = meta
            self._save_index()
        return path

    def show_saved(self) -> list[dict]:
        return [{"key": k, **v} for k, v in self._load_index().items()]

    # ------------------------------------------------------- 原始响应缓存

    def raw_path(self, method: str, url: str, body: bytes | None) -> Path:
        digest = hashlib.sha1(
            method.encode() + b"|" + url.encode() + b"|" + (body or b"")
        ).hexdigest()[:20]
        return self.raw_dir / f"{digest}.txt"

    def get_raw(self, method: str, url: str, body: bytes | None) -> str | None:
        if not self.enabled:
            return None
        path = self.raw_path(method, url, body)
        return path.read_text(encoding="utf-8", errors="replace") if path.is_file() else None

    def put_raw(self, method: str, url: str, body: bytes | None, text: str) -> None:
        if not self.enabled:
            return
        path = self.raw_path(method, url, body)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
