# SFTPC 基因 sgRNA 筛选工具

面向 **SFTPC 基因 CRISPR/Cas9 编辑**的 sgRNA 筛选工具：输入一段碱基序列，自动从
**CRISPOR** 和 **CHOPCHOP** 两个网站抓取数据，按三数据加权算出综合评分并排序。

课题背景：SFTPC 突变导致儿童间质性肺病（chILD），最常见的是 **I73T（c.218T>C）**。
本工具用于筛选"切得动、切得准、切点离突变位点近"的 sgRNA。

---

## 快速开始

```bash
# 网页版（推荐）
python run_web.py
# 打开 http://127.0.0.1:1145/

# 命令行
python run_crispor.py  --seq "ACGT..." --genome hg38 --pam NGG
python run_chopchop.py --seq "ACGT..." --genome hg38
```

零依赖，只用 Python 标准库（不需要 Flask）。

---

## 三个数据与权重

权重来自 `四个数据的定义.txt`：

| 权重 | 数据 | 来源 | 处理 |
|---|---|---|---|
| **7** | Doench '16 | CRISPOR | **均匀映射到 0–1**（高考赋分式：最好的 = 1.0，最差的 = 0.0，并列取平均名次） |
| **6** | CFD Spec. score | CRISPOR | **÷ 100** |
| **4** | 距离（第三项） | CHOPCHOP 的切割位点 | **分档** + 覆盖奖励，上限 1.0 |

```
总分 = (7 × d16' + 6 × cfd/100 + 4 × 第三项) / 17
```

**第三项**：

```
第三项 = 分档(|切点 − 22163096|) + (spacer 压住目标点 ? 0.1 : 0)      上限 1.0
分档：<20bp=1.0  20-50=0.9  50-100=0.7  100-200=0.4  >200=0.1
```

目标点 `chr8:22,163,096`（hg38）= SFTPC **c.218** = I73T 突变位点。

---

## 三个网站的 guide 怎么合并（负责人定的规则）

按 PAM 起点 / 20nt spacer 对齐后分成三堆：

| 部分 | 定义 | 处理 |
|---|---|---|
| **第一部分** | 两个网站**都有** | 三项正常算 |
| **第二部分** | CRISPOR 有、CHOPCHOP 没有 | 第三项：能真算就真算，否则用**第一部分第三项的平均值** |
| **第三部分** | CHOPCHOP 有、CRISPOR 没有 | **直接扔掉** |

**滚木(null)**：第一部分为空（两站完全没有共同的 sgRNA）时，第三项全部置 `null`，
权重按 7:6 归一化，页面给出提示。

详见 [`负责人指示.md`](负责人指示.md)。

---

## 目录结构

```
run_web.py                网页入口（默认 1145 端口）
run_crispor.py            CRISPOR 命令行
run_chopchop.py           CHOPCHOP 命令行
genome_utils.py           两个爬虫共用的坐标换算 / 序列定位

sgweb/                    网页
    server.py             HTTP + JSON 接口
    aggregate.py          三部分划分 + 加权排序
    index.html            前端（单文件，无外部依赖）

crispor_scraper/          网站1 爬虫
    client.py             CGI 接口：提交 / 探测 TSV / 下载
    parser.py             guides TSV -> 结构化记录
    datasets.py           两个数据的独立定义
    distance.py           第三项（CRISPOR 版）
    cache.py / pipeline.py / jsonio.py / cli.py

chopchop_scraper/         网站2 爬虫
    client.py             JSON API：POST / -> jobId
    parser.py             results.tsv -> 结构化记录
    distance.py           第三项（CHOPCHOP 版）
    locate.py             序列定位（转发 genome_utils）
    cache.py / pipeline.py / jsonio.py / cli.py
```

---

## 踩过的坑（都在代码里注释了）

### 1. 两个网站的"位置"锚点不一样 ⚠️

| 网站 | 报告值指的是 | 正链 guide 的号码是 |
|---|---|---|
| **CRISPOR** | **PAM 区间的左端** | PAM 起点 |
| **CHOPCHOP** | **23 mer 的最左端** | spacer 起点 |

**正链上两者差 20 bp。** 必须都先换算成 PAM 起点再算切点：

```
正链：切点 = PAM起点 − 3.5
反链：切点 = PAM起点 + 5.5
（SpCas9 切在 PAM 5' 侧 3 bp，即 protospacer 第 17/18 位之间）
```

两条链差 11 bp，漏掉 `strand` 会系统性算错。实现见 `genome_utils.py`。

### 2. CHOPCHOP 用序列提交时不给基因组坐标

- 基因名提交 → `chr8:22163178`（基因组坐标）
- 序列提交 → `seq:378`（**序列内坐标**，BED 里染色体名也是 `seq`）

换算公式：`基因组坐标 = genome_offset + N`（实测差 0）。
程序会自动把输入序列在本地定位回基因组，所以用户只需要给碱基。

### 3. CRISPOR 没有官方 API，但有 CGI 接口

```
提交   POST /crispor.py            seq / org / pam / name / submit
取结果 GET  /crispor.py?batchId=X&download=guides&format=tsv
```
比爬 HTML 稳。注意 `ajaxStatus=1` 那个轻量接口在服务器上是坏的（HTTP 500）。

---

## 缓存与速度

| 场景 | 耗时 |
|---|---|
| 全新查询（两个网站都要排队算） | 20 秒 ~ 3 分钟 |
| 同一个序列再跑 | **0.02 秒**（走磁盘缓存） |

缓存目录 `.crispor_cache/` / `.chopchop_cache/` / `.genome_cache/`（已 gitignore）。

**轮询策略**：不刷 265 KB 的结果页，直接反复探测几 KB 的 TSV 下载端点。

---

## 输入限制

- **CRISPOR 上限 2300 bp**（它自己的规矩）。整条 SFTPC 做不了，按外显子截。
- 一条 guide = 20nt spacer + 3nt PAM = **23 bp**，所以 23 bp 输入只能出 1 条候选。
  实测每约 5 bp 出 1 条候选（1600 bp → 312 条）。
- 输入必须是**基因组序列**，不能是 cDNA（跨外显子边界的 guide 没有真实靶点）。

## 测试序列

`SFTPC_测试序列_hg38_2213bp.txt` —— hg38 chr8:22,161,829–22,164,041，
SFTPC (NM_003018) 全部 6 个外显子编码区，含 c.218。直接粘进网页即可。

> 这个文件是**纯碱基、单行、无任何注释**（也没有署名头）—— 因为加了任何字符都会
> 导致粘贴进去报错。它跟 `LICENSE` 是仓库里仅有的两个不带署名头的文件。

预期结果：CRISPOR 434 条 / CHOPCHOP 416 条 / 第一部分 416 / 第二部分 18。

---

## 文档

- [`README_输入格式.md`](README_输入格式.md) —— 输入格式说明（碱基序列 + JSON）
- [`sgweb/README.md`](sgweb/README.md) —— 网页说明
- [`crispor_scraper/README.md`](crispor_scraper/README.md) —— 网站1 接口与坐标约定
- [`chopchop_scraper/README.md`](chopchop_scraper/README.md) —— 网站2 接口与坐标约定
- `四个数据的定义.txt` —— 权重定义

---

## 作者与许可

**本项目由华南师范大学附属中学知识城校区高二一班同学制作**

Made by students of **Class 1, Grade 11**, *The Affiliated High School of SCNU, Knowledge City Campus*
（华南师范大学附属中学知识城校区 高二一班）

- **Copyright (C) 2026 华南师范大学附属中学知识城校区高二一班**
- **许可证：GNU Affero General Public License v3.0 or later（AGPL-3.0-or-later）** —— 全文见 [LICENSE](LICENSE)

### 这意味着什么

| 你可以 | 你必须 |
|---|---|
| ✅ 自由使用、修改、分发 | 📌 保留作者署名和版权声明 |
| ✅ 用于商业用途 | 📌 衍生作品也必须用 AGPL-3.0 开源 |
| ✅ 私有修改（自己用） | 📌 **如果你把它做成网络服务提供给别人用，必须公开完整源码** |

> AGPL 比 GPL 严格的地方就在最后一条：**部署成在线服务也算"分发"**。
> 这是因为本项目的形态就是一个网页服务（`run_web.py`）。

### 引用

如果你在论文、报告或比赛中引用本项目，建议写成：

> 华南师范大学附属中学知识城校区高二一班. *SFTPC 基因 sgRNA 筛选工具*.
> 2026. https://github.com/ZZC-jeizuchin/sftpc-gene-sgRNA-selection
