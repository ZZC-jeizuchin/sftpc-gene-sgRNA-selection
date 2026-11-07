#!/usr/bin/env python3
# 本项目由华南师范大学附属中学知识城校区高二一班同学制作
# Copyright (C) 2026 华南师范大学附属中学知识城校区高二一班
# SPDX-License-Identifier: AGPL-3.0-or-later
"""端到端验证：启动网页服务 -> POST /api/analyze（v2 四维） -> 检查返回结构。"""

from __future__ import annotations

import json
import sys
import threading
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from sgweb import server  # noqa: E402

SEQ = "".join((ROOT / "SFTPC_测试序列_hg38_2213bp.txt").read_text(encoding="utf-8").split()).upper()

srv = server.ThreadingHTTPServer(("127.0.0.1", 0), server.Handler)
port = srv.server_address[1]
threading.Thread(target=srv.serve_forever, daemon=True).start()
time.sleep(0.4)
print(f"server on {port}", flush=True)

for allele in ("1", "0"):
    body = json.dumps({
        "sequence": SEQ, "genome": "hg38", "pam": "NGG",
        "target_pos": 22163096, "hard_filter": "1", "allele": allele,
    }).encode()
    req = urllib.request.Request(
        f"http://127.0.0.1:{port}/api/analyze", data=body,
        headers={"Content-Type": "application/json"})
    t0 = time.time()
    with urllib.request.urlopen(req, timeout=900) as resp:
        d = json.loads(resp.read())
    print("=" * 70)
    print(f"allele={allele}  ok={d.get('ok')}  elapsed={time.time()-t0:.1f}s  "
          f"model={d.get('model')}  count={d.get('count')}")
    print("warnings:", d.get("warnings"))
    print("steps:", [(s.get("site"), s.get("ok"), s.get("n_guides")) for s in d.get("steps", [])])
    rows = d.get("rows") or []
    if not rows:
        continue
    r = rows[0]
    print("top row keys present:",
          [k for k in ("allele_score", "allele_delta", "allele_region_label",
                       "allele_delta_100", "final_v2", "grade", "rank_v2",
                       "eff_score", "cfd_score", "dist_score") if k in r])
    print("top:", r.get("sequence"), "| v1", r.get("final_score"),
          "| v2", r.get("final_v2"), "| ASS", r.get("allele_score"),
          "| grade", r.get("grade"), "| rank_v2", r.get("rank_v2"))
    print("covered:", d.get("allele_covered_count"), " grade_counts:", d.get("grade_counts"))
    if allele == "1":
        print("weights_v2:", d.get("weights_v2"))
        print("allele_meta:", d.get("allele_meta"))
        (ROOT / "output" / "api_v2_sample.json").write_text(
            json.dumps(d, ensure_ascii=False), encoding="utf-8")
        print("saved output/api_v2_sample.json")

srv.shutdown()
print("DONE")
