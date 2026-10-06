# CRISPOR 抓取模块（网站1）

从 <https://crispor.gi.ucsc.edu/> 抓取**两个数据**：

| 数据 | CRISPOR 列名 | TSV 字段 | 含义 | 方向 |
|---|---|---|---|---|
| 数据1 | `Doench '16` | `Doench '16-Score` | on-target 编辑效率（Azimuth 2.0 / Rule Set 2），0–100 | 越大越好 |
| 数据2 | `CFD Spec. score` | `cfdSpecScore` | off-target 特异性（CFD 聚合，公式 `round(10000/(100+ΣCFD))`），0–100 | 越大越好 |

## 1. 有没有查询 API？

**没有官方 REST API**，但它前端用的是一套 CGI 参数接口，直接调这套接口比爬渲染后的 HTML 稳得多。
本项目用的是这套（已实测）：

| 用途 | 请求 | 说明 |
|---|---|---|
| 提交 | `POST /crispor.py`，字段 `seq` `org` `pam` `name` `submit=SUBMIT` | 返回 HTML，里面含 `?batchId=XXXX` |
| 取结果 | `GET /crispor.py?batchId=X&download=guides&format=tsv` | **6 KB TSV**，直接含两个数据 |
| 附加分 | 再加 `&showAllScores=1` | 多出 Chari / Xu / Wang / Doench'14 / CCTop |
| 脱靶明细 | `GET /crispor.py?batchId=X&download=offtargets&format=tsv` | 每个脱靶位点的 CFD 切割分 |
| 轻量状态 | `GET /crispor.py?batchId=X&ajaxStatus=1` | ⚠️ 服务器上已坏（HTTP 500），**本项目不用它** |

## 2. 三个 Step 的输入

对应 CRISPOR 首页那张表单，程序里就是 `CrisporInput` 的三个字段：

| Step | 首页字段 | 程序参数 | 例子 |
|---|---|---|---|
| Step 1 | Target sequence | `--seq` / `--fasta` / `--region` | `--region chr8:22,119,000-22,125,000` |
| Step 2 | Select a genome | `--genome` | `hg38`、`hg19`、`mm39` |
| Step 3 | Select a PAM | `--pam` | `NGG`、`NAG`、`NNGT` |

Step1 三选一：直接给序列（≤2300 bp）、给 FASTA 文件、给染色体区间（服务端自己去取序列）。

## 3. 缓存（重要：这个网站在国内很慢）

三层，越靠前越省：

1. **语义缓存命中**：`.crispor_cache/index.json` 里已有该输入指纹的 TSV →
   直接读本地，**0 次网络请求**（实测 0.1 秒，对比联网 4 秒、首次查询 30~60 秒）。
2. **batchId 复用**：有指纹但没 TSV → 跳过"提交"（最慢的一步），直接去取结果。
3. **原始响应缓存**：`.crispor_cache/raw/` 按 URL+参数 的 sha1 存；同一 URL 第二次访问直接读盘。

请求指纹 = `sha1(step1解析结果 | genome | pam | all_scores)`，**不含 `--name`**，所以换个名字不会重新请求。

另外：**不轮询 265 KB 的结果页**，而是每 3 秒重试那个 6 KB 的 TSV 下载端点——
算好了就一次拿到，没算好就快速失败（`retries=1, timeout=30`）。

## 4. 用法

```bash
# 正常查询（第一次慢，之后走缓存）
python run_crispor.py --region chr8:22119000-22119400 --genome hg38 --pam NGG --name sftpc

# 直接给序列
python run_crispor.py --seq ACGT...ACGT --genome hg38 --pam NGG

# 给 FASTA 文件
python run_crispor.py --fasta sftpc_exon.fa --genome hg38 --pam NGG

# 复用某个已知 batchId（跳过提交，最快）
python run_crispor.py --region chr8:22119000-22119400 --genome hg38 --pam NGG --from-batch Ukkl7GWmn9R8eqUN6GYe

# 看缓存里有什么
python run_crispor.py --list-cache

# 强制重新请求
python run_crispor.py --region chr8:22119000-22119400 --genome hg38 --pam NGG --no-cache

# 只要 shell 命令也可跑
python -m crispor_scraper.cli --region chr8:22119000-22119400 --genome hg38 --pam NGG
```

## 4b. JSON 输入 / 输出（两个数据分开）

### 输入 JSON

```json
{
  "step1": {"region": "chr8:22119000-22119400"},
  "step2": {"genome": "hg38"},
  "step3": {"pam": "NGG"},
  "name": "sftpc_exon3",
  "data": ["doench16", "cfd_spec"],
  "batch_id": "Ukkl7GWmn9R8eqUN6GYe",
  "options": {"all_scores": false, "cache": true,
              "timeout": 900, "poll_interval": 3, "retries": 3}
}
```

- `step1` 三选一：`{"seq": "..."}` / `{"region": "chr8:..."}` / `{"fasta": "x.fa"}`
- `step2` / `step3` 也可直接写成字符串：`"step2": "hg38"`、`"step3": "NGG"`
- `data` 缺省 = 全部；`batch_id` 缺省 = 正常提交
- 也接受**数组**（批量）或 `{"queries": [...]}`

### 输出 JSON（每个数据独立成块）

```json
{
  "ok": true,
  "query":   {"step1": {...}, "step2": {"genome": "hg38"}, "step3": {"pam": "NGG"}, "name": "..."},
  "crispor": {"batch_id": "...", "url": "...", "n_guides": 83,
              "from_cache": true, "reused_batch": true},
  "datasets": {
    "doench16": {
      "id": "doench16", "name": "Doench '16",
      "source": "crispor", "source_column": "Doench '16-Score",
      "range": [0, 100], "direction": "higher_is_better",
      "count": 83, "stats": {"min": 11, "max": 72, "mean": 48.08, "missing": 0},
      "values": [
        {"guide_id": "265rev", "strand": "rev", "cut_position": 265,
         "target_seq": "GCAGGTGTGAGCCCGCACCGGGG", "value": 72}
      ]
    },
    "cfd_spec": { "...": "同上，value 换成 CFD Spec. score" }
  }
}
```

### 命令

```bash
# JSON in -> JSON out（打到 stdout）
python run_crispor.py --input-json query.json --quiet

# 写一个完整 JSON
python run_crispor.py --input-json query.json --out-json out.json

# 每个数据分开各写一个文件：doench16.json / cfd_spec.json
python run_crispor.py --input-json query.json --out-dir json/

# 只要其中一个数据
python run_crispor.py --input-json query.json --data cfd_spec --quiet

# stdin
echo '{"step1":{"region":"chr8:22119000-22119400"},"step2":{"genome":"hg38"},"step3":{"pam":"NGG"}}' \
  | python run_crispor.py --input-json - --quiet

# 看有哪些数据可用
python run_crispor.py --list-data
```

`--quiet` 会关掉进度和表格，只吐 JSON，方便被别的程序调用。

## 5. 作为库调用

```python
from crispor_scraper import CrisporInput, fetch, jsonio, resolve

inp = CrisporInput(region="chr8:22119000-22119400", genome="hg38", pam="NGG", name="sftpc")
result = fetch(inp, on_progress=print)      # 自动走缓存

# 原始记录
for g in result.records:
    print(g.guide_id, g.doench16, g.cfd_spec)

# 按数据分开的 JSON
payload = jsonio.result_to_json(result, resolve(["doench16", "cfd_spec"]))
print(payload["datasets"]["doench16"]["values"][0])
print(payload["datasets"]["cfd_spec"]["stats"])

result.save("crispor_out")   # 可选：另存原始 TSV
print(result.url)
```

## 6. 模块结构

```
run_crispor.py              入口脚本
crispor_scraper/
    config.py   Step1/2/3 的定义、校验、转成表单字段
    client.py   HTTP：提交 / 探测 TSV / 下载 / 下载脱靶表
    cache.py    语义缓存(index.json) + 原始响应缓存(raw/)
    parser.py   guides TSV -> GuideRecord（按表头名取列，官方加列不会解析错）
    pipeline.py fetch() 高层流水线 + CrisporResult
    datasets.py **两个数据各自的独立定义**（id / 来源列 / 方向 / 提取函数）
    jsonio.py   JSON 输入解析 + JSON 输出构造（按数据拆块）
    cli.py      命令行
```

### 加第三个数据（距离）时怎么改

只在 `datasets.py` 里加一个 `DatasetSpec`（填 id / 名称 / 来源 / 说明 / 提取函数），
然后 `DATASETS` 注册进去即可。`client`、`parser`、`cache`、`jsonio` 一行都不用动 ——
JSON 输出会自动多出一个块，`--data` 也会自动认这个新 id。

## 7. 已知边界

- 输入序列必须是基因组序列，**不能是 cDNA**；跨外显子边界的 guide 没有真实靶点。
- 序列 ≤ 2300 bp（只放一个外显子/目标区域），序列两端 50 bp 内的 guide 算不出 Doench '16（显示 `--`）。
- 落在重复序列里的 guide，CRISPOR 会把 MIT / CFD 都置 0（不是抓取错误）。
- `Doench '16` 在 CRISPOR 里是 `pam_audit=False`，非 NGG 的 guide 也给分，但模型没在那些 PAM 上验证过。
- 全基因组脱靶只搜到 4 个错配；NGG guide 还会把 NAG/NGA PAM 的位点计入 ΣCFD。
