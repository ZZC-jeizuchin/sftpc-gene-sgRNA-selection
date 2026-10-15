# 本项目由华南师范大学附属中学知识城校区高二一班同学制作
# Copyright (C) 2026 华南师范大学附属中学知识城校区高二一班
# SPDX-License-Identifier: AGPL-3.0-or-later
"""磁盘缓存，两层：

    index.json  请求指纹 -> batchId + TSV 文件名（已有 TSV 就直接返回）
    raw/        URL 的 sha1 -> 原始响应；只缓存下载类请求，轮询不走这层
"""

from __future__ import annotations

import hashlib
import json
import time
from dataclasses import asdict
from pathlib import Path

DEFAULT_CACHE_DIR = ".crispor_cache"


class CrisporCache:
    def __init__(self, cache_dir: str | Path = DEFAULT_CACHE_DIR, enabled: bool = True):
        self.enabled = enabled
        self.root = Path(cache_dir)
        self.raw_dir = self.root / "raw"
        self.index_path = self.root / "index.json"
        self._index: dict | None = None

    # 索引读写

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
            json.dumps(self._load_index(), ensure_ascii=False, indent=2),
            encoding="utf-8",
        )

    # 请求指纹

    @staticmethod
    def request_key(seq: str, org: str, pam: str, all_scores: bool) -> str:
        blob = f"{seq}|{org}|{pam}|{int(bool(all_scores))}"
        return hashlib.sha1(blob.encode("utf-8")).hexdigest()[:16]

    @staticmethod
    def key_from_input(inp) -> str:
        return CrisporCache.request_key(
            inp.resolve_step1(), inp.genome, inp.pam, inp.show_all_scores
        )

    def _entry(self, key: str) -> dict:
        return self._load_index().get(key, {})

    # 语义缓存 get/put

    def get_batch_id(self, key: str) -> str | None:
        return self._entry(key).get("batch_id") if self.enabled else None

    def put_batch_id(self, key: str, batch_id: str, describe: dict | None = None) -> None:
        if not self.enabled:
            return
        idx = self._load_index()
        entry = idx.setdefault(key, {})
        entry["batch_id"] = batch_id
        entry["updated"] = time.strftime("%Y-%m-%d %H:%M:%S")
        if describe:
            entry["input"] = describe
        self._save_index()

    def get_tsv(self, key: str) -> str | None:
        """命中则返回 TSV 文本。"""
        if not self.enabled:
            return None
        fname = self._entry(key).get("tsv_file")
        if not fname:
            return None
        path = self.root / fname
        if not path.is_file():
            return None
        return path.read_text(encoding="utf-8", errors="replace")

    def put_tsv(self, key: str, batch_id: str, tsv: str, describe: dict | None = None) -> Path:
        self.root.mkdir(parents=True, exist_ok=True)
        fname = f"guides_{key}.tsv"
        path = self.root / fname
        path.write_text(tsv, encoding="utf-8")
        if self.enabled:
            idx = self._load_index()
            entry = idx.setdefault(key, {})
            entry.update({
                "batch_id": batch_id,
                "tsv_file": fname,
                "saved": time.strftime("%Y-%m-%d %H:%M:%S"),
                "n_bytes": len(tsv),
            })
            if describe:
                entry["input"] = describe
            self._save_index()
        return path

    def get_aux(self, name: str) -> str | None:
        path = self.root / name
        return path.read_text(encoding="utf-8", errors="replace") if path.is_file() else None

    def put_aux(self, name: str, text: str) -> Path:
        self.root.mkdir(parents=True, exist_ok=True)
        path = self.root / name
        path.write_text(text, encoding="utf-8")
        return path

    def show_saved(self) -> list[dict]:
        idx = self._load_index()
        rows = []
        for key, entry in idx.items():
            rows.append({"key": key, **entry})
        return rows

    # 原始响应缓存

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
