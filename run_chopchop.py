#!/usr/bin/env python3
"""入口脚本：python run_chopchop.py --target SFTPC --genome hg38 --out-json d3.json"""

import sys

from chopchop_scraper.cli import main

if __name__ == "__main__":
    sys.exit(main())
