"""命令行入口：把三个 Step 作为参数传进来，抓两个数据并打印。"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import jsonio
from .client import CrisporError
from .config import COMMON_PAMS, CrisporConfigError, CrisporInput
from .datasets import DATASETS, DEFAULT_ORDER, resolve
from .pipeline import fetch


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="run_crispor.py",
        description="从 CRISPOR（网站1）抓取两个数据：Doench '16 与 CFD Spec. score。",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "三个 Step 的输入（对应 CRISPOR 首页表单）：\n"
            "  Step1  --seq / --fasta / --region   目标序列（<2300 bp，通常一个外显子）\n"
            "  Step2  --genome                     基因组，如 hg38 / hg19\n"
            "  Step3  --pam                        PAM，如 NGG / NAG\n\n"
            "输入输出都可以是 JSON：\n"
            "  python run_crispor.py --input-json query.json --out-json out.json\n"
            "  python run_crispor.py --input-json query.json --out-dir json/   # 每个数据一个文件\n"
            "  echo '{\"step1\":{\"region\":\"chr8:22119000-22119400\"},\"step2\":{\"genome\":\"hg38\"},\"step3\":{\"pam\":\"NGG\"}}' \\\n"
            "       | python run_crispor.py --input-json - --quiet\n\n"
            "  python run_crispor.py --region chr8:22119000-22119400 --genome hg38 --pam NGG\n"
            "  python run_crispor.py --fasta sftpc_exon.fa --genome hg38 --pam NGG --data cfd_spec\n"
        ),
    )

    # ---- JSON 输入输出 ----
    gj = p.add_argument_group("JSON 输入 / 输出")
    gj.add_argument("--input-json", metavar="FILE",
                    help="输入 JSON 文件；给 - 表示从 stdin 读")
    gj.add_argument("--out-json", metavar="FILE", help="把完整结果写成单个 JSON 文件")
    gj.add_argument("--out-dir", metavar="DIR",
                    help="把**每个数据分开**各写一个 JSON 文件（doench16.json / cfd_spec.json）")
    gj.add_argument("--data", metavar="IDS",
                    help=f"只要哪些数据，逗号分隔。可选 {sorted(DATASETS)}，默认 {'+'.join(DEFAULT_ORDER)}")
    gj.add_argument("--quiet", action="store_true",
                    help="不打进度和表格，只输出 JSON（给程序调用时用）")

    # ---- Step 1 ----
    g1 = p.add_argument_group("Step 1  Target sequence（三选一）")
    g1.add_argument("--seq", help="直接给 DNA 序列（ACGTN，<2300 bp）")
    g1.add_argument("--fasta", help="FASTA 文件路径（取第一条序列）")
    g1.add_argument("--region", help="染色体区间，如 chr8:22,119,000-22,125,000")
    g1.add_argument("--name", default="query", help="序列名（仅用于输出文件名，默认 query）")

    # ---- Step 2 ----
    g2 = p.add_argument_group("Step 2  Genome")
    g2.add_argument("--genome", default="hg38", help="基因组，如 hg38 / hg19（默认 hg38）")

    # ---- Step 3 ----
    g3 = p.add_argument_group("Step 3  PAM")
    g3.add_argument("--pam", default="NGG", help="PAM 类型，如 NGG / NAG（默认 NGG）")

    # ---- 运行参数 ----
    g4 = p.add_argument_group("运行参数（缓存优先）")
    g4.add_argument("--all-scores", action="store_true",
                    help="同时下载 Chari/Xu/Wang/Doench'14 等附加评分列")
    g4.add_argument("--cache-dir", default=".crispor_cache",
                    help="缓存目录（默认 .crispor_cache）")
    g4.add_argument("--no-cache", action="store_true", help="禁用缓存，强制重新请求")
    g4.add_argument("--from-batch", help="直接复用某个已有的 batchId，跳过提交")
    g4.add_argument("--list-cache", action="store_true", help="列出已缓存的结果后退出")
    g4.add_argument("--list-data", action="store_true", help="列出可用的数据 id 后退出")
    g4.add_argument("--timeout", type=float, default=900.0, help="最长等待秒数（默认 900）")
    g4.add_argument("--poll-interval", type=float, default=3.0, help="轮询间隔秒数（默认 3）")
    g4.add_argument("--retries", type=int, default=3, help="网络重试次数（默认 3）")
    g4.add_argument("--target-pos", type=int, default=22_163_096,
                    help="第三个数据的目标位点（hg38），默认 22163096 = SFTPC c.218/I73T")
    g4.add_argument("--genome-offset", type=int, default=None,
                    help="输入序列的基因组 0-based 起点；不给则自动定位")
    g4.add_argument("--no-auto-locate", action="store_true",
                    help="只给碱基时不自动定位 offset")
    g4.add_argument("--list-options", action="store_true",
                    help="打印常用 PAM / 基因组取值后退出")

    return p


def _print_table(rows: list[dict], limit: int | None = None) -> None:
    if not rows:
        print("（没有解析到任何 guide）")
        return
    shown = rows if limit is None else rows[:limit]
    header = ("rank", "guideId", "strand", "cutPos", "targetSeq(20+PAM)",
              "Doench'16", "CFD_Spec", "MIT_Spec", "offTargets")
    widths = [4, 12, 6, 6, 24, 9, 8, 8, 10]

    def fmt(cells):
        out = []
        for cell, w in zip(cells, widths):
            text = "-" if cell is None else str(cell)
            out.append(text.ljust(w))
        return "  ".join(out)

    print(fmt(header))
    print("-" * (sum(widths) + 2 * (len(widths) - 1)))
    for i, row in enumerate(shown, 1):
        print(fmt([i, row["guide_id"], row["strand"], row["cut_position"],
                   row["target_seq"], row["doench16"], row["cfd_spec"],
                   row["mit_spec"], row["offtarget_count"]]))


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if args.list_options:
        print("常用 Step3 PAM：")
        for k, v in COMMON_PAMS.items():
            print(f"  {k:<8} {v}")
        print("\nStep2 genome 直接填 UCSC 基因组名，如 hg38 / hg19 / mm39 / mm10。")
        return 0

    if args.list_data:
        for i in DEFAULT_ORDER:
            spec = DATASETS[i]
            print(f"  {spec.id:<10} {spec.name:<18} 来源={spec.source:<8} "
                  f"列={spec.source_column:<18} 范围={spec.value_min}-{spec.value_max}")
        print("\n用 --data doench16,cfd_spec 选择；缺省=全部。")
        return 0

    if args.list_cache:
        from .cache import CrisporCache
        cache = CrisporCache(args.cache_dir)
        rows = cache.show_saved()
        if not rows:
            print(f"缓存为空：{cache.root}")
            return 0
        print(f"缓存目录：{cache.root}")
        for row in rows:
            print(f"  key={row['key']}  batchId={row.get('batch_id')}  "
                  f"tsv={row.get('tsv_file', '-')}  saved={row.get('saved', '-')}")
        return 0

    if not args.input_json and not any([args.seq, args.fasta, args.region]):
        print("错误：必须给 --input-json，或 Step1 的 --seq / --fasta / --region 之一。\n",
              file=sys.stderr)
        build_parser().print_help(sys.stderr)
        return 2

    # ---------------- 1. 收集查询 ----------------
    if args.input_json:
        try:
            text = sys.stdin.read() if args.input_json == "-" else \
                Path(args.input_json).read_text(encoding="utf-8")
            inputs, data_ids = jsonio.load_queries(text)
        except (CrisporConfigError, OSError) as exc:
            print(jsonio.dumps(jsonio.error_to_json(str(exc), kind="input")), file=sys.stderr)
            return 2
    else:
        inputs = [CrisporInput(
            sequence=args.seq,
            fasta_path=args.fasta,
            region=args.region,
            name=args.name,
            genome=args.genome,
            pam=args.pam,
            show_all_scores=args.all_scores,
            poll_interval=args.poll_interval,
            timeout=args.timeout,
            retries=args.retries,
            batch_id=args.from_batch,
            use_cache=not args.no_cache,
            target_pos=args.target_pos,
            genome_offset=args.genome_offset,
            auto_locate=not args.no_auto_locate,
        )]
        data_ids = []

    if args.data:
        data_ids = [d.strip() for d in args.data.split(",") if d.strip()]
    try:
        specs = resolve(data_ids)
    except KeyError as exc:
        print(jsonio.dumps(jsonio.error_to_json(str(exc), kind="data")), file=sys.stderr)
        return 2

    # ---------------- 2. 逐个抓取 ----------------
    payloads: list[dict] = []
    for inp in inputs:
        try:
            result = fetch(
                inp,
                on_progress=None if args.quiet else (lambda m: print(f"[crispor] {m}", flush=True)),
                cache_dir=args.cache_dir,
            )
        except CrisporConfigError as exc:
            payloads.append(jsonio.error_to_json(str(exc), jsonio.query_to_json(inp), "input"))
            continue
        except CrisporError as exc:
            payloads.append(jsonio.error_to_json(str(exc), jsonio.query_to_json(inp), "fetch"))
            continue
        payloads.append(jsonio.result_to_json(result, specs))
        if not args.quiet:
            _print_table(result.two_data())

    output = payloads[0] if len(payloads) == 1 else payloads

    # ---------------- 3. 输出 JSON ----------------
    if args.out_json:
        jsonio.write_json(args.out_json, output)
        if not args.quiet:
            print(f"\n[JSON] {args.out_json}", file=sys.stderr)
    if args.out_dir:
        _write_split_json(args.out_dir, payloads, specs, quiet=args.quiet)
    if not args.out_json and not args.out_dir:
        print(jsonio.dumps(output))
    return 0


def _write_split_json(out_dir: str, payloads: list[dict], specs, quiet: bool = False) -> None:
    """把每个数据写成一个独立文件（这就是"分开"）。"""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    written: list[str] = []
    for i, payload in enumerate(payloads, 1):
        if not payload.get("ok"):
            continue
        suffix = "" if len(payloads) == 1 else f"_{i}"
        for spec in specs:
            block = payload["datasets"].get(spec.id)
            if block is None:
                continue
            path = out / f"{spec.id}{suffix}.json"
            jsonio.write_json(path, {
                "ok": True,
                "query": payload["query"],
                "crispor": payload["crispor"],
                **block,
            })
            written.append(str(path))
    if written and not quiet:
        print("[JSON] 分开保存：" + "  ".join(written), file=sys.stderr)


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
