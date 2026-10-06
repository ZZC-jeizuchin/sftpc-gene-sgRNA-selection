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
"""网页后端：标准库 http.server，零依赖。

    GET  /                前端页面
    GET  /api/health      健康检查
    POST /api/analyze     {"sequence": "...", "genome": "hg38", "pam": "NGG",
                           "target_pos": 22163096, "refresh": false}
                          -> 跑两个爬虫 -> 三数据 -> 加权排序 -> JSON
"""

from __future__ import annotations

import json
import sys
import threading
import time
import traceback
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from genome_utils import DEFAULT_TARGET_POS  # noqa: E402
from sgweb.aggregate import merge  # noqa: E402

INDEX_HTML = _HERE / "index.html"
LAST_RESULT: dict | None = None      # 上次分析结果，供"加载上次结果"用
MAX_SEQ = 10_000


def run_pipeline(sequence: str, genome: str, pam: str, target_pos: int,
                 log=None, genome_offset: int | None = None,
                 genome_chrom: str | None = None,
                 third_mode: str = "auto",
                 locate_window: tuple[str, int, int] | None = None) -> dict:
    """跑两个爬虫 + 合并加权。任一网站失败都能降级出结果。"""
    def say(msg):
        if log:
            log(msg)

    sequence = "".join(sequence.split()).upper()
    result: dict = {"ok": True, "warnings": [], "steps": []}
    t0 = time.time()

    # 定位窗口（默认 SFTPC 基因座；换基因时传进来）
    loc_kw = {}
    if locate_window:
        chrom, lstart, lend = locate_window
        loc_kw = {"locus_chrom": chrom, "locus_start": lstart, "locus_end": lend}
        say(f"定位窗口：{chrom}:{lstart}-{lend}")

    crispor_payload = None
    chopchop_payload = None

    # ---- 网站1：Doench '16 + CFD Spec. (+ 距离回退) ----
    try:
        say("网站1 CRISPOR：提交中…")
        from crispor_scraper import CrisporInput, fetch as crispor_fetch
        inp1 = CrisporInput(
            sequence=sequence, genome=genome, pam=pam, name="web",
            target_pos=target_pos, poll_interval=3.0, timeout=600.0,
            genome_offset=genome_offset, genome_chrom=genome_chrom,
            auto_locate=(genome_offset is None), **loc_kw,
        )
        r1 = crispor_fetch(inp1)
        crispor_payload = r1.distance_payload()
        result["steps"].append({
            "site": "crispor", "ok": True, "batch_id": r1.batch_id,
            "url": r1.url, "n_guides": len(r1.records), "from_cache": r1.from_cache,
        })
        say(f"网站1 完成：{len(r1.records)} 条 guide")
    except Exception as exc:
        result["warnings"].append(f"CRISPOR 失败：{exc}")
        result["steps"].append({"site": "crispor", "ok": False, "error": str(exc)})
        say(f"网站1 失败：{exc}")

    # ---- 网站2：距离 ----
    try:
        say("网站2 CHOPCHOP：提交中…")
        from chopchop_scraper import ChopchopInput, fetch as chopchop_fetch
        inp2 = ChopchopInput(
            sequence=sequence, genome=genome, sequence_name="web",
            poll_interval=3.0, timeout=600.0,
            genome_offset=genome_offset, genome_chrom=genome_chrom,
            auto_locate=(genome_offset is None), **loc_kw,
        )
        r2 = chopchop_fetch(inp2)
        chopchop_payload = r2.distance_payload(target_pos)
        result["steps"].append({
            "site": "chopchop", "ok": True, "job_id": r2.job_id,
            "url": r2.url, "n_guides": len(r2.records), "from_cache": r2.from_cache,
        })
        say(f"网站2 完成：{len(r2.records)} 条 guide")
    except Exception as exc:
        result["warnings"].append(f"CHOPCHOP 失败：{exc}")
        result["steps"].append({"site": "chopchop", "ok": False, "error": str(exc)})
        say(f"网站2 失败：{exc}")

    if crispor_payload is None and chopchop_payload is None:
        return {"ok": False, "error": "两个网站都失败了，无法计算。",
                "warnings": result["warnings"], "steps": result["steps"]}

    merged = merge(crispor_payload, chopchop_payload, target_pos, third_mode=third_mode)

    # ---- 定位情况 ----
    off = genome_offset
    for pl in (crispor_payload, chopchop_payload):
        if pl and off is None:
            off = pl.get("genome_offset")
    merged["resolved_offset"] = off
    merged["resolved_chrom"] = genome_chrom or next(
        (pl.get("genome_chrom") for pl in (crispor_payload, chopchop_payload) if pl and pl.get("genome_chrom")),
        None)

    parts = merged.get("parts", {})

    # ---- 滚木：第一部分为空 ----
    if merged.get("kurumi"):
        merged.setdefault("warnings", []).append(
            "**滚木(null)** —— " + (merged.get("kurumi_note") or
            "两个网站没有共同的 sgRNA，第三项无法计算。")
        )
    # ---- 定位失败 ----
    if off is None and not merged.get("kurumi"):
        merged.setdefault("warnings", []).append(
            "这段序列**不在定位窗口内**，拿不到基因组坐标 → 距离项无法计算。"
            "如果这是别的基因，请在左侧「序列基因组起点」填入起始坐标，"
            "或改大「定位窗口」；如果本来就想测 SFTPC，点下面的按钮换成示例序列。"
        )
        merged["need_offset"] = True

    # ---- 信息提示（不阻止结果） ----
    if parts.get("part3_dropped"):
        merged.setdefault("notes", []).append(
            f"有 {parts['part3_dropped']} 条只出现在 CHOPCHOP 的 sgRNA 已按规则丢弃"
            f"（它们没有 Doench'16 / CFD）。"
        )
    if parts.get("part2") and merged.get("avg_part1_third") is not None:
        n_avg = sum(1 for r in merged["rows"] if r.get("third_source") == "avg_part1")
        if n_avg:
            if third_mode == "average":
                why = "按负责人规则（第二部分一律用第一部分平均值）"
            else:
                why = "这些 guide 两个网站上对不上、算不出切点"
            merged.setdefault("notes", []).append(
                f"第二部分有 {n_avg} 条{why}，第三项用第一部分的平均值 "
                f"({merged['avg_part1_third']}) 顶替。"
            )

    merged["elapsed_sec"] = round(time.time() - t0, 1)
    merged["ok"] = True
    merged["warnings"] = list(result["warnings"]) + list(merged.get("warnings", []))
    merged["steps"] = result["steps"]
    merged["sequence_length"] = len(sequence)
    merged["sequence_preview"] = sequence[:60] + ("..." if len(sequence) > 60 else "")
    return merged


class Handler(BaseHTTPRequestHandler):
    server_version = "sgRNA-Designer/1.0"

    def log_message(self, fmt, *args):
        sys.stderr.write("[web] " + (fmt % args) + "\n")

    # ------------------------------------------------------------ helpers

    def _send(self, code: int, body: bytes, ctype: str = "application/json; charset=utf-8"):
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def _json(self, code: int, obj):
        self._send(code, json.dumps(obj, ensure_ascii=False).encode("utf-8"))

    # ------------------------------------------------------------ routes

    def do_GET(self):
        path = self.path.split("?", 1)[0]
        if path in ("/", "/index.html"):
            try:
                html = INDEX_HTML.read_bytes()
            except OSError:
                self._send(500, b"index.html missing", "text/plain; charset=utf-8")
                return
            self._send(200, html, "text/html; charset=utf-8")
        elif path == "/api/health":
            self._json(200, {"ok": True, "service": "sgRNA Designer"})
        elif path == "/api/last":
            if LAST_RESULT is None:
                self._json(404, {"ok": False, "error": "本次会话还没有结果"})
            else:
                self._json(200, LAST_RESULT)
        else:
            self._json(404, {"ok": False, "error": "not found"})

    def do_POST(self):
        path = self.path.split("?", 1)[0]
        if path != "/api/analyze":
            self._json(404, {"ok": False, "error": "not found"})
            return
        try:
            length = int(self.headers.get("Content-Length") or 0)
            raw = self.rfile.read(length) if length else b"{}"
            params = json.loads(raw.decode("utf-8"))
        except Exception as exc:
            self._json(400, {"ok": False, "error": f"请求体不是合法 JSON：{exc}"})
            return

        sequence = (params.get("sequence") or "").strip()
        if not sequence:
            self._json(400, {"ok": False, "error": "请填写碱基序列。"})
            return
        clean = "".join(sequence.split()).upper()
        bad = sorted(set(c for c in clean if c not in "ACGTN"))
        if bad:
            self._json(400, {"ok": False, "error": f"序列里有非 ACGTN 字符：{bad}"})
            return
        if len(clean) > MAX_SEQ:
            self._json(400, {"ok": False, "error": f"序列过长（{len(clean)} bp），上限 {MAX_SEQ} bp。"})
            return

        genome = (params.get("genome") or "hg38").strip()
        pam = (params.get("pam") or "NGG").strip()
        try:
            target_pos = int(params.get("target_pos") or DEFAULT_TARGET_POS)
        except (TypeError, ValueError):
            self._json(400, {"ok": False, "error": "target_pos 必须是整数。"})
            return

        logs: list[str] = []
        goff = params.get("genome_offset")
        try:
            goff = int(goff) if goff not in (None, "", "auto") else None
        except (TypeError, ValueError):
            self._json(400, {"ok": False, "error": "序列基因组起点必须是整数。"})
            return

        # 定位窗口，形如 "chr8:22150000-22170000"（换基因时用）
        loc_win = None
        raw_win = (params.get("locate_window") or "").strip()
        if raw_win:
            import re as _re
            m = _re.match(r"^\s*([A-Za-z0-9_.\-]+)\s*:\s*([\d,]+)\s*-\s*([\d,]+)\s*$", raw_win)
            if not m:
                self._json(400, {"ok": False,
                                 "error": f"定位窗口格式不对：{raw_win!r}，应形如 chr8:22150000-22170000"})
                return
            loc_win = (m.group(1),
                       int(m.group(2).replace(",", "")),
                       int(m.group(3).replace(",", "")))

        third_mode = (params.get("third_mode") or "auto").strip()
        if third_mode not in ("auto", "average", "computed"):
            self._json(400, {"ok": False, "error": "third_mode 必须是 auto / average / computed"})
            return

        try:
            out = run_pipeline(clean, genome, pam, target_pos, log=logs.append,
                               genome_offset=goff,
                               genome_chrom=(params.get("genome_chrom") or None),
                               third_mode=third_mode,
                               locate_window=loc_win)
        except Exception:
            self._json(500, {"ok": False, "error": "内部错误：\n" + traceback.format_exc()})
            return
        out["logs"] = logs
        global LAST_RESULT
        LAST_RESULT = out
        self._json(200, out)


def serve(port: int = 1145, host: str = "127.0.0.1"):
    httpd = ThreadingHTTPServer((host, port), Handler)
    httpd.daemon_threads = True
    print(f"sgRNA Designer 已启动: http://{host}:{port}", flush=True)
    return httpd


if __name__ == "__main__":
    p = int(sys.argv[1]) if len(sys.argv) > 1 else 1145
    srv = serve(p)
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        print("\nbye")
