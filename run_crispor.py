#!/usr/bin/env python3
"""入口脚本：python run_crispor.py --region chr8:... --genome hg38 --pam NGG"""

import sys

from crispor_scraper.cli import main

if __name__ == "__main__":
    sys.exit(main())
