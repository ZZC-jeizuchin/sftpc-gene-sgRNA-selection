"""CHOPCHOP 命令行：抓第三个数据（切点 → chr8:22163096 的距离）。"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import jsonio
from .client import ChopchopError
from .config import (
    APPLICATIONS, ChopchopConfigError, ChopchopInput, COMMON_GENOMES,
    NUCLEASES, SCORING_METHODS, TARGET_REGIONS,
)
from .distance import DEFAULT_TARGET_POS, DISTANCE_BANDS
from .pipeline import fetch


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="run_chopchop.py",
        description="从 CHOPCHOP（网站2）抓第三个数据：切割位点到 chr8:22163096 的距离。",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "四个输入（对应 CHOPCHOP 首页四个框）：\n"
            "  Target   --target SFTPC   或 --target chr8:22157000-22165000  或 --seq ACGT...\n"
            "  In       --genome hg38\n"
            "  Using    --nuclease CRISPR\n"
            "  For      --application knockout\n\n"
            "示例：\n"
            "  python run_chopchop.py --target SFTPC --genome hg38 --out-json d3.json\n"
            "  python run_chopchop.py --input-json query.json --quiet\n"
            "  python run_chopchop.py --target SFTPC --target-pos 22163096 --quiet\n"
        ),
    )

    gj = p.add_argument_group("JSON 输入 / 输出")
    gj.add_argument("--input-json", metavar="FILE", help="输入 JSON；- 表示 stdin")
    gj.add_argument("--out-json", metavar="FILE", help="把结果写成 JSON")
    gj.add_argument("--out-dir", metavar="DIR", help="把数据块分开写成 JSON（按 id 命名）")
    gj.add_argument("--quiet", action="store_true", help="只输出 JSON")

    g1 = p.add_argument_group("Target（三选一）")
    g1.add_argument("--target", help="基因名（SFTPC）或染色体区间（chr8:22157000-22165000）")
    g1.add_argument("--seq", help="直接给序列")
    g1.add_argument("--name", default="query", help="序列名（仅输出用）")

    g2 = p.add_argument_group("In / Using / For")
    g2.add_argument("--genome", default="hg38", help="基因组，默认 hg38")
    g2.add_argument("--nuclease", default="CRISPR", choices=sorted(NUCLEASES))
    g2.add_argument("--application", default="knock-out", choices=sorted(APPLICATIONS))
    g2.add_argument("--target-region", default="CODING", choices=sorted(TARGET_REGIONS),
                    help="靶向区域，默认 CODING（官网默认）")

    g3 = p.add_argument_group("第三个数据 / Options")
    g3.add_argument("--genome-offset", type=int, default=None,
                    help="用 --seq 提交时，这段序列在基因组上的 0-based 起始坐标；"
                         "不给就算不出基因组距离")
    g3.add_argument("--no-auto-locate", action="store_true",
                    help="用序列提交时不自动去基因组定位 offset")
    g3.add_argument("--genome-chrom", default=None, help="序列模式下的染色体名，如 chr8（仅显示用）")
    g3.add_argument("--target-pos", type=int, default=DEFAULT_TARGET_POS,
                    help=f"目标位点（hg38 正链坐标），默认 {DEFAULT_TARGET_POS} = SFTPC c.218/I73T")
    g3.add_argument("--scoring", default="DOENCH_2016", choices=sorted(SCORING_METHODS))
    g3.add_argument("--gc-min", type=int, default=10)
    g3.add_argument("--gc-max", type=int, default=90)
    g3.add_argument("--guide-len", type=int, default=20)

    g4 = p.add_argument_group("运行参数")
    g4.add_argument("--cache-dir", default=".chopchop_cache")
    g4.add_argument("--job-id", help="复用已有 jobId，跳过提交")
    g4.add_argument("--timeout", type=float, default=600.0)
    g4.add_argument("--poll-interval", type=float, default=3.0)
    g4.add_argument("--retries", type=int, default=3)
    g4.add_argument("--list-cache", action="store_true")
    g4.add_argument("--list-options", action="store_true")
    return p


def _print_distance_table(payload: dict, limit: int | None = None) -> None:
    vals = payload["values"]
    if limit:
        vals = sorted(vals, key=lambda v: (v["distance"] is None, v["distance"]))[:limit]
    header = ("#", "rank", "location", "链", "切点", "距离(bp)", "分", "档", "序列(20+PAM)")
    widths = (4, 5, 11, 3, 12, 9, 5, 10, 24)
    def fmt(cells):
        return "  ".join(str("-" if c is None else c).ljust(w) for c, w in zip(cells, widths))
    print(fmt(header))
    print("-" * (sum(widths) + 2 * (len(widths) - 1)))
    for i, v in enumerate(vals, 1):
        loc = f"{v['chrom']}:{v['location']}" if v["location"] else "-"
        print(fmt([i, v["rank"], loc, v["strand"], v["cut_site"],
                   v["distance"], v["score"], v["band"], v["target_seq"]]))


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)

    if args.list_options:
        print("nuclease:", ", ".join(sorted(NUCLEASES)))
        print("application:", ", ".join(sorted(APPLICATIONS)))
        print("target_region:", ", ".join(sorted(TARGET_REGIONS)))
        print("scoring:", ", ".join(sorted(SCORING_METHODS)))
        print("genome:", ", ".join(sorted(COMMON_GENOMES)))
        print("\n距离分档：", ", ".join(
            f"<{int(u)}bp={s}" if u != float("inf") else f">200bp={s}" for u, s in DISTANCE_BANDS))
        return 0

    if args.list_cache:
        from .cache import ChopchopCache
        cache = ChopchopCache(args.cache_dir)
        rows = cache.show_saved()
        print(f"缓存目录：{cache.root}")
        if not rows:
            print("  （空）")
        for r in rows:
            print(f"  key={r['key']}  jobId={r.get('job_id')}  tsv={r.get('tsv_file','-')}  "
                  f"saved={r.get('saved','-')}")
        return 0

    # ---------------- 收集查询 ----------------
    if args.input_json:
        try:
            text = sys.stdin.read() if args.input_json == "-" else \
                Path(args.input_json).read_text(encoding="utf-8")
            queries = jsonio.load_queries(text)
        except (ChopchopConfigError, OSError) as exc:
            print(jsonio.dumps(jsonio.error_to_json(str(exc), kind="input")), file=sys.stderr)
            return 2
    else:
        if not args.target and not args.seq:
            print("错误：必须给 --input-json，或 --target / --seq。\n", file=sys.stderr)
            build_parser().print_help(sys.stderr)
            return 2
        inp = ChopchopInput(
            target=args.target, sequence=args.seq, sequence_name=args.name,
            genome=args.genome, nuclease=args.nuclease, application=args.application,
            target_region=args.target_region, scoring=args.scoring,
            gc_min=args.gc_min, gc_max=args.gc_max, guide_len=args.guide_len,
            poll_interval=args.poll_interval, timeout=args.timeout,
            retries=args.retries, job_id=args.job_id,
            genome_offset=args.genome_offset, genome_chrom=args.genome_chrom,
            auto_locate=not args.no_auto_locate,
        )
        queries = [(inp, args.target_pos)]

    # ---------------- 抓取 ----------------
    payloads = []
    for inp, target_pos in queries:
        try:
            result = fetch(
                inp,
                on_progress=None if args.quiet else (lambda m: print(f"[chopchop] {m}", flush=True)),
                cache_dir=args.cache_dir,
            )
        except ChopchopConfigError as exc:
            payloads.append(jsonio.error_to_json(str(exc), None, "input"))
            continue
        except ChopchopError as exc:
            payloads.append(jsonio.error_to_json(str(exc), None, "fetch"))
            continue
        payloads.append(jsonio.result_to_json(result, target_pos))
        if not args.quiet:
            _print_distance_table(payloads[-1]["datasets"]["distance_to_22163096"])

    output = payloads[0] if len(payloads) == 1 else payloads

    if args.out_json:
        jsonio.write_json(args.out_json, output)
        if not args.quiet:
            print(f"\n[JSON] {args.out_json}", file=sys.stderr)
    if args.out_dir:
        out = Path(args.out_dir)
        out.mkdir(parents=True, exist_ok=True)
        written = []
        for i, payload in enumerate(payloads, 1):
            if not payload.get("ok"):
                continue
            suffix = "" if len(payloads) == 1 else f"_{i}"
            for did, block in payload["datasets"].items():
                path = out / f"{did}{suffix}.json"
                jsonio.write_json(path, {
                    "ok": True, "query": payload["query"],
                    "chopchop": payload["chopchop"], **block,
                })
                written.append(str(path))
        if written and not args.quiet:
            print("[JSON] 分开保存：" + "  ".join(written), file=sys.stderr)
    if not args.out_json and not args.out_dir:
        print(jsonio.dumps(output))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
