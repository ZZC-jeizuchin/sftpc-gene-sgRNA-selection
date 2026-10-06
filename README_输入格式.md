# 输入格式说明：一段碱基

**两个爬虫（CRISPOR / CHOPCHOP）都只要求一段碱基序列。** 基因组/坐标/偏移全部由程序自己处理。

---

## 1. 最简用法

```bash
# 网站1：Doench '16 + CFD Spec.
python run_crispor.py --seq "CTTCCTTTGTCCCCAATCTG..." --genome hg38 --pam NGG --quiet

# 网站2：到 chr8:22163096 的距离
python run_chopchop.py --seq "CTTCCTTTGTCCCCAATCTG..." --genome hg38 --quiet
```

`--seq` 可以换成长度不限的一段碱基（ACGTN，忽略空白与换行）。
CRISPOR 有 **2300 bp** 硬限制（它的规矩），CHOPCHOP 没明说。

---

## 2. 两个爬虫现在干同样的事

| | CRISPOR | CHOPCHOP |
|---|---|---|
| 输入 | 一段碱基 | 一段碱基 |
| 脱靶搜索 | 它服务器干（≤4 错配） | 它服务器干（≤3 错配） |
| 输出数据 | `doench16`、`cfd_spec`、`distance_to_22163096` | `distance_to_22163096` |
| 联网 | 首次 ~30 s，之后命中缓存 0.1 s | 首次 ~10 s，之后 0.06 s |

CRISPOR 也输出第三个数据（距离），这样可以跟 CHOPCHOP 互相校验。

---

## 3. ⚠️ 唯一容易错的地方：两个网站的"位置"锚点不同

这是本项目最坑的一点，已经统一处理，但要知道：

| 网站 | 报告值指的是 | 正链 guide 的号码是 |
|---|---|---|
| **CRISPOR** | **PAM 区间的左端**（不管链） | PAM 起点 |
| **CHOPCHOP** | **23 mer 的最左端**（不管链） | spacer 起点 |

**正链上两者差 20 bp。**

统一做法：先都换算成 **PAM 起点**，再算切点。

```
正链: 切点 = PAM起点 - 3.5
反链: 切点 = PAM起点 + 5.5
（SpCas9 切在 PAM 5' 侧 3 bp，即 protospacer 第 17/18 位之间）
```

实现都在工作区根目录的 `genome_utils.py`，**两个爬虫共用同一份**，不会各写一套。

JSON 输出里每一步都写了：
- `anchor` —— 该网站用的锚点（`"pam"` / `"23mer"`）
- `reported_position` —— 网站原样给的数字
- `pam_start` —— 归一化后的 PAM 起点（**跨网站对齐用这个**）
- `cut_site` / `distance` / `score` / `band`

---

## 4. 坐标自动定位（所以只需要给碱基）

CHOPCHOP 用"贴序列"方式提交时只回 `seq:378`（序列内坐标），不给基因组坐标。
第三个数据要基因组坐标，所以程序会：

1. 把 SFTPC 基因座（hg38 chr8:22.15–22.17 Mb，20 kb）拉下来缓存到本地
2. 把输入的碱基序列在本地匹配回去（支持 ≤5 个错配、支持反链）
3. 得到 `genome_offset`，然后 `基因组坐标 = genome_offset + 序列内位置`

实测能力：

| 情况 | 结果 |
|---|---|
| 原序列 | ✅ 精确匹配 |
| 任意位置有突变（含第 1 位、最后一位） | ✅ 1 个错配也能定位 |
| 最多 5 个错配 | ✅ |
| 反链 | ✅ |
| 完全无关的序列 | ❌ 正确报错，不瞎猜 |

**不需要下载人类全基因组**（3 GB）—— 只要目标基因座那一小段（20 kb）。
以后换基因，改 `genome_utils.py` 里的 `DEFAULT_LOCUS`，或者传 `--genome-offset` 手动指定。

---

## 5. JSON 输入输出

两个爬虫都支持：

```bash
python run_crispor.py  --input-json q1.json --out-json o1.json --quiet
python run_chopchop.py --input-json q2.json --out-json o2.json --quiet
```

输入（CRISPOR）：

```json
{
  "step1": {"seq": "CTTCCTTT..."},
  "step2": {"genome": "hg38"},
  "step3": {"pam": "NGG"},
  "name": "sftpc",
  "target_pos": 22163096,
  "options": {"cache": true, "auto_locate": true}
}
```

输入（CHOPCHOP）：

```json
{
  "sequence": "CTTCCTTT...",
  "genome": "hg38",
  "nuclease": "CRISPR",
  "application": "knock-out",
  "target_region": "CODING",
  "sequence_name": "sftpc",
  "target_pos": 22163096,
  "options": {"cache": true, "auto_locate": true}
}
```

输出都是：

```json
{
  "ok": true,
  "query": {...},
  "crispor|chopchop": {"batch_id|job_id": "...", "url": "...", "n_guides": 134},
  "datasets": {
    "doench16":  {...},
    "cfd_spec":  {...},
    "distance_to_22163096": {
      "anchor": "pam", "target_pos": 22163096,
      "stats": {"min": 0.5, "within_20bp": 9, ...},
      "values": [
        {"guide_id": "278forw", "reported_position": 278, "pam_start": 22163099,
         "cut_site": 22163095.5, "distance": 0.5, "score": 1.0, "band": "<=20bp"}
      ]
    }
  }
}
```

---

## 6. 交叉验证结果

同一段 600 bp 的 SFTPC 序列，同时喂给两个爬虫：

```
CRISPOR 134 条 / CHOPCHOP 128 条 / 同一 PAM 位点 128 个
切点/距离不一致：0 / 128  → 完全一致 ✓
```

---

## 7. 文件

```
genome_utils.py            两个爬虫共用的坐标/定位逻辑（唯一真相来源）

run_crispor.py             CRISPOR 入口
crispor_scraper/           网站1：Doench'16 / CFD / 距离

run_chopchop.py            CHOPCHOP 入口
chopchop_scraper/          网站2：距离
```

---

> **本项目由华南师范大学附属中学知识城校区高二一班同学制作**
> Made by students of Class 1, Grade 11, The Affiliated High School of SCNU (Knowledge City Campus)
> Copyright (C) 2026 华南师范大学附属中学知识城校区高二一班 · 许可证：[AGPL-3.0-or-later](LICENSE)
