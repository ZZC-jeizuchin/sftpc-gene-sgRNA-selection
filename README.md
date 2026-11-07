# SFTPC 基因 sgRNA 筛选工具

输入一段碱基序列，自动去 CRISPOR 和 CHOPCHOP 两个网站抓数据，按四个维度加权排序，
用来筛 SFTPC 基因 CRISPR 编辑的 sgRNA。

**两个模型并存**：v1 三维修（7:6:4，只看单条序列）、v2 四维修（35:30:25:10，
多算一维「等位基因特异性」）。网页默认走 v2，开关关掉就退回 v1。

背景是 SFTPC 的 **I73T 突变（c.218T>C）**，这个突变会导致儿童间质性肺病。
挑 sgRNA 主要看三件事：切得动、切得准、切点离突变位点近。

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

不用装任何第三方库，Python 标准库就够（不需要 Flask）。

---

## v1：三个数据与权重（7:6:4）

权重来自 `四个数据的定义.txt`。这一套仍然完整保留，关掉第四维就回到它：

| 权重 | 数据 | 来源 | 处理 |
|---|---|---|---|
| **7** | Doench '16 | CRISPOR | **均匀映射到 0-1**（高考赋分式：最好的 = 1.0，最差的 = 0.0，并列取平均名次） |
| **6** | CFD Spec. score | CRISPOR | **/ 100** |
| **4** | **突变距离评分** | CHOPCHOP 的切割位点 | **分档** + 覆盖奖励，上限 1.0 |

```
总分 = (7 x d16' + 6 x cfd/100 + 4 x 突变距离评分) / 17
```

> **显示约定**：内部一律用 0-1 计算；网页上 **所有分数按 0-100 显示、保留两位小数**
> （比如内部 0.9192 -> 页面显示 `91.92`）。

**突变距离评分**：

```
突变距离评分 = 分档(|切点 - 22163096|) + (spacer 压住目标点 ? 0.1 : 0)      上限 1.0
分档：<20bp=1.0  20-50=0.9  50-100=0.7  100-200=0.4  >200=0.1
```

**切点坐标约定**（`genome_utils.CUT_CONVENTION`，默认 `"int"`）：

| 约定 | 碱基当成 | 切点公式（用 PAM 起点 p） | 距离形态 |
|---|---|---|---|
| `"int"`（默认） | 区间 | 正链 `p-3`，反链 `p+6` | 整数 |
| `"half"` | 点 | 正链 `p-3.5`，反链 `p+5.5` | 半点 |

两种约定指的是**同一个物理切点**，只是坐标系平移了 0.5；实测 434 条 guide 的
**分档结果完全一致**。

目标点 `chr8:22,163,096`（hg38）= SFTPC **c.218** = I73T 突变位点。

---

## v2：第四维「等位基因特异性」（默认开启）

v1 只看一条序列（参考基因组），所以它回答不了「这条 sgRNA 会不会把野生型也切了」。
v2 多算一维：把**只改了一个碱基的突变型序列**再提交一次 CRISPOR，
比较两条序列上的 Doench '16 得分。

```
Δ_raw  = Doench(突变型) − Doench(野生型)          // −100 … +100
ΔScore = 100 × clamp(Δ_raw, 0, 20) / 20           // 20 分差算拉满

位置分：突变落在 guide 的哪一段
    PAM 区(21-23) 100 / Seed 区(1-8) 80 / 中间区(9-14) 40 / PAM 远端(15-20) 10

ASS = 0.6 × ΔScore + 0.4 × 位置分
ΔScore < 20（两条切得差不多）时整体 ×0.3
```

**四维权重**

| 维度 | 权重 | 说明 |
|---|---|---|
| 等位基因特异性 (ASS) | **35** | 能不能区分突变型和野生型 |
| 切割效率 | **30** | Doench '16 |
| 脱靶特异性 | **25** | CFD Spec |
| 突变距离 | **10** | 从 v1 的 23.5% 降下来，只当辅助参考 |

> 只有识别序列（23mer）盖住突变位点的 guide 才可能有非零的 ASS。
> 本案例 434 条里只有 7 条满足，所以综合分结构上上限约 65 分 ——
> 等级因此改用**分位**（前 5% / 20% / 50% / 其余），不用固定分数线。

## 两个网站的 guide 怎么合并（负责人定的规则）

按 PAM 起点 / 20nt spacer 对齐后分成三堆：

| 部分 | 定义 | 处理 |
|---|---|---|
| **第一部分** | 两个网站**都有** | 三项正常算 |
| **第二部分** | CRISPOR 有、CHOPCHOP 没有 | 第三项：能真算就真算，否则用**第一部分第三项的平均值** |
| **第三部分** | CHOPCHOP 有、CRISPOR 没有 | **直接扔掉** |

**滚木(null)**：第一部分为空（两站完全没有共同的 sgRNA）时，第三项全部置 `null`，
权重按 7:6 归一化，页面给出提示。

---

## 目录结构

```
run_web.py                网页入口（默认 1145 端口）
run_crispor.py            CRISPOR 命令行
run_chopchop.py           CHOPCHOP 命令行
genome_utils.py           两个爬虫共用的坐标换算 / 序列定位

sgweb/                    网页
    server.py             HTTP + JSON 接口
    aggregate.py          两站合并 + v1 加权 + v2 四维重排
    allele_specific.py    v2 第四维：ΔScore + 位置分
    index.html            前端（单文件，无外部依赖）

crispor_scraper/          网站1 爬虫
chopchop_scraper/         网站2 爬虫
deploy/                   systemd + nginx + 一键部署

tools/                    取数与复现
    acquire.py            抓野生型数据
    acquire_mutant.py     抓突变型数据（构造"只改一个碱基"的序列）
    acquire_case.py       通用版：任意基因 / 任意突变
    run_all.py            一键复现全流程
    md_to_docx.py         Markdown 转 Word

analysis/                 分析与验证
    score_v2.py           四维评分 + 排名表
    weight_analysis.py    客观赋权 + 消融 + 稳定性
    validate_external.py  用 CFD 独立模型做外部佐证
    simulate_experiment.py 体外切割模拟
    case_compare.py       换基因（HBB 负链）验证通用性

tests/
    test_web_v2.py        起服务 + 打接口 + 校字段

data/                     抓取快照 + CFD 错配矩阵（提交，保证离线可复现）
output/                   分析产物（每次跑都会重写，已 gitignore）
```

---

## 踩过的坑（都在代码里注释了）

### 1. 两个网站的"位置"锚点不一样

| 网站 | 报告值指的是 | 正链 guide 的号码是 |
|---|---|---|
| **CRISPOR** | **PAM 区间的左端** | PAM 起点 |
| **CHOPCHOP** | **23 mer 的最左端** | spacer 起点 |

**正链上这两个数差 20 bp，反链上又正好重合。** 所以必须先换算成同一个基准
（本工具统一用 PAM 起点）再算切点：

```
切点 = PAM起点 - 3（正链）  /  PAM起点 + 6（反链）
SpCas9 切在 PAM 5' 侧 3 bp；反链方向是反的，偏移不一样。
```

漏掉 `strand` 会系统性算错。实现见 `genome_utils.py`。

### 2. CHOPCHOP 用序列提交时不给基因组坐标

- 基因名提交 -> `chr8:22163178`（基因组坐标）
- 序列提交 -> `seq:378`（**序列内坐标**，BED 里染色体名也是 `seq`）

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
| 同一个序列再跑 | **约 0.03 秒**（走磁盘缓存，线上约 0.3 秒） |

缓存目录 `.crispor_cache/` / `.chopchop_cache/` / `.genome_cache/`（已 gitignore）。

**轮询策略**：不刷 265 KB 的结果页，直接反复探测几 KB 的 TSV 下载端点。

---

## 输入限制

- **CRISPOR 上限 2300 bp**（它自己的限制）。整条 SFTPC 放不下，得按外显子截。
- 一条 guide = 20nt spacer + 3nt PAM = **23 bp**，所以只给 23 bp 就只能出 1 条候选。
  实测差不多每 5 bp 出 1 条（1600 bp -> 312 条）。
- 输入必须是**基因组序列**，不能是 cDNA（跨外显子边界的 guide 没有真实靶点）。

## 测试序列

`SFTPC_测试序列_hg38_2213bp.txt` —— hg38 chr8:22,161,829-22,164,041，
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
| 自由使用、修改、分发 | 保留作者署名和版权声明 |
| 用于商业用途 | 衍生作品也必须用 AGPL-3.0 开源 |
| 私有修改（自己用） | **如果你把它做成网络服务提供给别人用，必须公开完整源码** |

> AGPL 比 GPL 严格的地方就在最后一条：**部署成在线服务也算"分发"**。
> 这是因为本项目的形态就是一个网页服务（`run_web.py`）。

### 引用

如果你在论文、报告或比赛中引用本项目，建议写成：

> 华南师范大学附属中学知识城校区高二一班. *SFTPC 基因 sgRNA 筛选工具*.
> 2026. https://github.com/ZZC-jeizuchin/sftpc-gene-sgRNA-selection
