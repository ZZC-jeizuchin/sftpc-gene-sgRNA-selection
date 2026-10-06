# CHOPCHOP 抓取模块（网站2）—— 第三个数据

从 <https://chopchop.cbu.uib.no/> 抓 **第三个数据：切割位点到 chr8:22,163,096 的距离**。

## 0. 这个位点是什么

```
hg38 chr8:22,163,096
  = SFTPC NM_003018  的 c.218   （exon 3，cDNA 第 249 位）
  = 密码子 73 的第 2 位（codon 73 = ATT = Ile）
  = 最常见的 chILD 突变 I73T (c.218T>C) 的位点
```

⚠️ **这个坐标只对应 hg38**。hg19 的 SFTPC 在 chr8:22,019,310–22,021,991，两者差 14 万 bp，用错基因组所有距离都会落到 `>200bp` 档。

## 1. 接口：CHOPCHOP 有正规 JSON API

比网站1 还干净，全部实测：

| 用途 | 请求 |
|---|---|
| **提交** | `POST /`，`Content-Type: application/json` |
| **运行信息** | `GET /results/<jobId>/run.info` |
| **结果表（核心）** | `GET /results/<jobId>/results.tsv` |
| 请求回显 | `GET /results/<jobId>/query.json` |
| 绘图坐标 | `GET /results/<jobId>/cutcoords.json` |
| BED | `GET /results/<jobId>/results.bed` |

提交体（和网页 `formController.js` 拼出来的一模一样）：

```json
{
  "opts": ["-J","-BED","-GenBank","-G","hg38","-filterGCmin","10","-filterGCmax","90",
           "-t","WHOLE","-n","BspQI","-R","4",
           "-T","1","-g","20","-scoringMethod","DOENCH_2016","-f","NN","-v","3","-w"],
  "fastaInput": "",
  "geneInput": "SFTPC",
  "isIsoform": false,
  "forSelect": "knockout"
}
```

返回 `{"jobId": "a37d4ff6-..."}`，结果在 `/results/<jobId>/`。

`results.tsv` 的列：

```
Rank  Target sequence  Genomic location  Strand  GC content (%)
Self-complementarity  MM0  MM1  MM2  MM3  Efficiency
```

## 2. ⚠️ 坐标约定（关键，已实测 161/161 通过）

> **`Genomic location` = 23 bp 靶位点（20 nt spacer + PAM）在正链上的最左端 1-based 坐标。**

它**不是切割位点**。SpCas9 切在 PAM 前 3 bp，所以：

| 链 | 切割位点（正链坐标） |
|---|---|
| `+` | `loc + 16` 与 `loc + 17` 之间 → **`loc + 16.5`** |
| `−` | `loc + 5` 与 `loc + 6` 之间 → **`loc + 5.5`** |

**两条链差 11 bp。** 文档里写的"22163096 取绝对值"如果直接拿 `loc` 去减，`+` 链和 `−` 链的 guide 会拿到系统性不同的距离。

实测对比（SFTPC / hg38 / WHOLE，161 条）：

| | 距离最近的 5 条 |
|---|---|
| **正确（算切点）** | 22163079+ (0.5) · 22163078+ (1.5) · 22163077+ (2.5) · 22163085+ (5.5) · 22163104− (13.5) |
| 错误（直接用 loc） | 22163103+ (7) · 22163104− (8) · 22163085+ (11) · 22163111− (15) · 22163112− (16) |

**排名完全不同。** JSON 输出里同时给了 `distance`（正确）和 `distance_naive_nostrand`（对照，别用）。

验证方法：拉 hg38 全序列，对每条 guide 取 `loc` 起的 23 bp，`+` 链比对正链、`−` 链比对反向互补 —— 161 条全中。

## 3. 距离分档（照文档）

| 距离 | 分数 |
|---|---|
| < 20 bp | 1.0 |
| 20–50 bp | 0.9 |
| 50–100 bp | 0.7 |
| 100–200 bp | 0.4 |
| > 200 bp | 0.1 |

在 `distance.py` 的 `DISTANCE_BANDS` 里，改一处即可。

## 4. 输入（= CHOPCHOP 首页四个框 + Options）

| 首页 | 参数 | 说明 |
|---|---|---|
| Target | `--target` | 基因名 `SFTPC` 或区间 `chr8:22157000-22165000`；或 `--seq` 直接给序列 |
| In | `--genome` | `hg38`（默认） |
| Using | `--nuclease` | `CRISPR`（默认）/ `CAS13` / `CPF1` / `NICKASE` / `TALEN` |
| For | `--application` | `knockout`（默认） |
| Options | `--target-region` | `WHOLE`（默认）/ `CODING` / `SPLICE` / `UTR5` / `UTR3` / `PROMOTER` |
| | `--scoring` | `DOENCH_2016`（默认）等 |
| | `--target-pos` | 目标位点，默认 `22163096` |

## 5. 用法

```bash
# 抓 SFTPC 全基因的 guides，算到 c.218 的距离
python run_chopchop.py --target SFTPC --genome hg38 --out-json d3.json

# 只看编码区（guide 更少，7 秒出结果）
python run_chopchop.py --target SFTPC --target-region CODING --out-json d3.json

# 复用已有 jobId（跳过提交）
python run_chopchop.py --target SFTPC --job-id a37d4ff6-096e-4653-9202-6ae34e3f1dba --quiet

# 缓存
python run_chopchop.py --list-cache
```

输出（JSON，第三个数据独立成块）：

```json
{
  "ok": true,
  "query": {...},
  "chopchop": {"job_id": "...", "url": "...", "n_guides": 161, "from_cache": true},
  "datasets": {
    "distance_to_22163096": {
      "id": "distance_to_22163096",
      "name": "距离（切点 → chr8:22163096）",
      "source": "chopchop",
      "direction": "lower_is_better",
      "unit": "bp",
      "target_pos": 22163096,
      "bands": [...],
      "count": 161,
      "stats": {"min": 0.5, "max": 1380.5, "mean": 817.21,
                "within_20bp": 9, "score_sum": 39.6},
      "values": [
        {"guide_id": "22163178-", "rank": 1, "chrom": "chr8", "location": 22163178,
         "strand": "-", "target_seq": "GCTGGTAGTCATACACCACGAGG",
         "cut_site": 22163183.5, "distance": 87.5, "score": 0.7, "band": "50-100bp",
         "distance_naive_nostrand": 82}
      ]
    }
  }
}
```

## 6. 实测速度

| 场景 | 耗时 |
|---|---|
| 全新查询（提交 → 轮询 → 下载 → 解析） | **7 秒**（CODING）/ ~10 秒（WHOLE） |
| 复用 jobId 重下 | 2.2 秒 |
| 缓存命中 | **0.06 秒**，0 网络请求 |

轮询策略和网站1 一样：**不刷结果页**，直接每 3 秒重试那个小的 `results.tsv`。

## 7. 模块结构

```
run_chopchop.py               入口脚本
chopchop_scraper/
    config.py    Target/Genome/Nuclease/For + Options 的定义、校验、拼 opts
    client.py    HTTP：POST 提交 / 探测 results.tsv / run.info / query.json
    cache.py     index.json（指纹→jobId）+ raw/（原始响应缓存）
    parser.py    results.tsv -> GuideRecord
    distance.py  **第三个数据**：切点计算 + 区间打分 + JSON 打包
    pipeline.py  fetch() + ChopchopResult
    jsonio.py    JSON 输入解析 / JSON 输出
    cli.py       命令行
```

## 8. 已知边界

- `Genomic location` 是 23 mer 左端，**任何用到"位置"的地方都要先加链偏移**（这是本模块最容易错的一点）。
- CHOPCHOP 只搜到 **3 个错配**（`MM0`–`MM3`），比网站1 浅；它给的是脱靶**计数**，没有 CFD 分。
- `Efficiency` 默认是 Doench 2016，和网站1 的 `Doench '16` **同源但实现不同**，数值不完全一样，别混用。
- 以序列（`fastaInput`）提交时结果里仍会给基因组坐标（若序列能比中），但比不中就没有坐标。
