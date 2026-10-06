#!/usr/bin/env python3
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
