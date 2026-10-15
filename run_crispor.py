#!/usr/bin/env python3
# 本项目由华南师范大学附属中学知识城校区高二一班同学制作
# Copyright (C) 2026 华南师范大学附属中学知识城校区高二一班
# SPDX-License-Identifier: AGPL-3.0-or-later
"""入口脚本：python run_crispor.py --region chr8:... --genome hg38 --pam NGG"""

import sys

from crispor_scraper.cli import main

if __name__ == "__main__":
    sys.exit(main())
