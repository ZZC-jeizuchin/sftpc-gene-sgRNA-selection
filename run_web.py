#!/usr/bin/env python3
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
