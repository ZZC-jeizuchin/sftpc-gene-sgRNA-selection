#!/usr/bin/env python3
# 本项目由华南师范大学附属中学知识城校区高二一班同学制作
# Copyright (C) 2026 华南师范大学附属中学知识城校区高二一班
# SPDX-License-Identifier: AGPL-3.0-or-later
"""启动网页服务。用法： python run_web.py [端口]   默认 1145"""
import sys
from sgweb.server import serve

if __name__ == "__main__":
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 1145
    srv = serve(port)
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        print("\nbye")
