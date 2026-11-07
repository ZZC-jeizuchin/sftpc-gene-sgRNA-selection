========================================================
sgRNA 智能筛选工具 —— 交付说明
========================================================

一、这是什么
--------------------------------------------------------
一个网页工具：把 CRISPOR 和 CHOPCHOP 两个网站的 sgRNA 候选自动采集回来，
统一坐标口径，按四个维度加权排序，帮研究者挑出"最有可能只切突变型"的几条。

配套还有一份参赛报告和一份课题全解（都在本目录下，Word 格式）。


二、最快看到效果（不需要联网、不需要装任何包）
--------------------------------------------------------
环境要求：Python 3.10 或更高版本（Windows / macOS / Linux 都行）

    cd 本项目目录
    python run_web.py 1145

然后用浏览器打开：

    http://127.0.0.1:1145

页面里点"换成 SFTPC 示例序列"→ 点"开始分析"，就能看到完整结果。
网页工具是纯 Python 标准库实现的，零第三方依赖。


三、想复现全部研究结果
--------------------------------------------------------
    pip install -r requirements.txt      # 只有两个包：openpyxl、Pillow

    python tools/run_all.py --skip-net        # 用已下载的数据跑全部分析（推荐）
    python tools/run_all.py                   # 连取数一起跑（要联网，约 6 分钟）

产物都在 output/ 目录下：排名表、CFD 验证、权重分析、模拟凝胶图等。


四、只想知道每个文件是干什么的
--------------------------------------------------------
【网页本体】
  run_web.py                         入口：起网页服务
  sgweb/server.py                    后端（标准库 http.server，零依赖）
  sgweb/index.html                   前端（单文件）
  sgweb/aggregate.py                 ★核心：两站数据合并 + 四维加权 + 分位等级
  sgweb/allele_specific.py           ★核心：等位基因特异性评分（ΔScore + 位置分）

【两个网站的接口反向工程】
  crispor_scraper/                   CRISPOR（无官方 API，靠分析前端请求）
  chopchop_scraper/                  CHOPCHOP（同上）
  genome_utils.py                    两者共用的坐标换算
  run_crispor.py / run_chopchop.py   两站各自的命令行入口

【取数与复现】tools/
  acquire.py                         抓 SFTPC 野生型数据
  acquire_mutant.py                  抓 SFTPC 突变型数据（构造"只改一个碱基"的序列）
  acquire_case.py                    通用版：任意基因 / 任意单碱基突变
  run_all.py                         一键复现：取数 -> 评分 -> 赋权 -> 验证 -> 模拟
  md_to_docx.py                      Markdown 转 Word（生成申报材料用）

【分析与验证】analysis/
  score_v2.py                        四维评分，输出排名表
  weight_analysis.py                 CRITIC 客观赋权 + 消融实验 + 稳定性分析
  validate_external.py               用 CFD 独立模型做外部佐证
  simulate_experiment.py             体外切割实验的计算机模拟 + 凝胶图
  case_compare.py                    多基因案例对比（SFTPC 正链 / HBB 负链）

【测试】tests/
  test_web_v2.py                     起服务、打接口、校验字段

【部署】
  deploy/                            systemd 单元 + nginx 配置 + 一键部署脚本

【文档】
  参赛报告_终稿.docx                  论文体的参赛报告（含四次迭代过程）
  课题全解_思路与方法.docx             从零理解这个课题
  修改方案_v2.docx / P0_执行报告.docx  阶段性工作记录
  实验模拟与外部佐证.docx              实验模拟那一轮的说明
  申报材料/                          查新报告、研究日志模板（需签字盖章）
  docs_src/                          上述 Word 的 Markdown 源文件

【数据与结果】
  data/                              抓取快照 + CFD 错配矩阵（提交，保证离线可复现）
  output/                            分析产物（每次跑都会重写，已 gitignore）


五、数据来源
--------------------------------------------------------
· CRISPOR — https://crispor.gi.ucsc.edu/
· CHOPCHOP — https://chopchop.cbu.uib.no/
· CFD 错配矩阵 — Doench et al. 2016 补充材料，经
  github.com/MichaelLinn/off_target_prediction 公开分发
· 基因组序列与基因结构 — UCSC Genome Browser API、Ensembl REST API

上述数据/工具版权归原团队所有，本项目仅用于科研与竞赛展示。


六、许可证
--------------------------------------------------------
AGPL-3.0-or-later，见 LICENSE。
任何人修改后如果通过网络提供服务，也必须以相同许可证开源。


七、已知限制（请务必先看）
--------------------------------------------------------
1. 没有湿实验验证。报告里所有结论都是计算预测；
   analysis/simulate_experiment.py 产出的是"预期结果"，不是实验数据。
2. 权重 35/30/25/10 是目标驱动的经验判断，不是从数据拟合出来的。
3. 只覆盖 SpCas9 + NGG。
4. 当前版本不独立做全基因组脱靶搜索，脱靶分数直接取 CRISPOR 的 CFD。
