# 部署到 sgss.743.world

本项目由华南师范大学附属中学知识城校区高二一班同学制作 · AGPL-3.0-or-later

---

## 架构

```
用户 → https://sgss.743.world      nginx :443（宝塔管理 + 自动 SSL）
           ↓ proxy_pass
        http://127.0.0.1:1146      本项目 Python 服务（systemd 常驻）
```

**Python 服务只绑 `127.0.0.1`** —— 对外只暴露 nginx 的 80/443，
所以**不需要在阿里云安全组或宝塔防火墙开任何新端口**，最省事也最安全。

> ⚠️ **端口必须是 1146，不能用 1145** —— 服务器上 1145 已被 `veridrop`
> （uvicorn）占用。空闲端口：1146 / 1147 / 8765 / 9000。

---

## 部署前必须确认

| # | 事项 | 谁做 |
|---|---|---|
| 1 | **DNS**：`sgss.743.world` 加 A 记录 → `8.138.245.240` | 域名持有者 |
| 2 | 阿里云安全组放行 80/443（现有站点已放行，通常无需改动） | 域名持有者 |

DNS 没生效前，`certbot` / 宝塔申请 SSL 会失败，nginx 也匹配不到这个域名。

验证 DNS 是否生效：

```bash
dig +short sgss.743.world          # 应该返回 8.138.245.240
```

---

## 部署步骤

### 方式 A：一键脚本（推荐）

```bash
ssh -p 2536 root@8.138.245.240
cd /home/ZZC-jeizuchin/sftpc-sgrna && git pull    # 或先 clone
bash deploy/deploy.sh
```

脚本会做：检查端口 → 拉最新代码 → 装 systemd 单元 → 重启服务 →
本地冒烟测试 → 若 nginx 配置已存在则 `nginx -t && nginx -s reload`。

### 方式 B：宝塔面板（图形界面）

1. **建站**：面板 → 网站 → 添加站点
   - 域名：`sgss.743.world`
   - 根目录：随便（比如 `/www/wwwroot/sgss.743.world`）
   - PHP 版本：**纯静态**（我们不用 PHP）
   - 数据库：不创建

2. **申请 SSL**：站点设置 → SSL → Let's Encrypt → 申请 → 开启「强制 HTTPS」

3. **配反向代理**：站点设置 → 反向代理 → 添加反向代理
   - 代理名称：`sftpc-sgrna`
   - 目标 URL：`http://127.0.0.1:1146`
   - 发送域名：`$host`

4. **⚠️ 改超时**（关键，不改的话长查询会 502/504）：
   编辑 `/www/server/panel/vhost/nginx/sgss.743.world.conf`，
   在 `location` 块里加：
   ```nginx
   proxy_connect_timeout 60s;
   proxy_send_timeout    900s;
   proxy_read_timeout    900s;
   ```
   然后 `nginx -t && nginx -s reload`

---

## 为什么超时要设 900 秒

一次查询要**依次**做这些事：

| 步骤 | 耗时 |
|---|---|
| 提交 CRISPOR + 轮询到算完 | 5 ~ 90 s |
| 提交 CHOPCHOP + 轮询到算完 | 5 ~ 60 s |
| 本地序列定位（20 kb 基因座） | < 1 s（首次会拉一次 UCSC） |

实测：
- 23 bp 短序列 → **14.9 s**
- 2213 bp 的完整 SFTPC CDS → 本地 191 s，服务器 **20.8 s**（网站端有缓存）

nginx 默认 `proxy_read_timeout` 是 **60 秒**，不改必然超时。

> 服务器在阿里云广州，`chopchop.cbu.uib.no` 首次响应要 20 秒，
> 比本机慢，所以超时更要留足。

---

## 运维

```bash
systemctl status sftpc-sgrna          # 状态
systemctl restart sftpc-sgrna         # 重启
journalctl -u sftpc-sgrna -f          # 实时日志
journalctl -u sftpc-sgrna -n 100      # 最近 100 行
ss -ltnp | grep 1146                  # 端口确认（应只有 127.0.0.1）
```

### 缓存

缓存在工作目录下，**删了不影响运行，只是下次要重新抓**：

```bash
du -sh /home/ZZC-jeizuchin/sftpc-sgrna/.?*cache
rm -rf /home/ZZC-jeizuchin/sftpc-sgrna/.crispor_cache \
       /home/ZZC-jeizuchin/sftpc-sgrna/.chopchop_cache
```

模板缓存（`.genome_cache/locus/*.txt`，20 KB）**建议保留** ——
它让"序列 → 基因组坐标"的定位不用每次都请求 UCSC。
换基因时改定位窗口即可，或用 `--locate-window` 参数。

---

## 更新代码

```bash
cd /home/ZZC-jeizuchin/sftpc-sgrna
git pull
systemctl restart sftpc-sgrna
```

---

## 排错

| 症状 | 原因 / 处理 |
|---|---|
| 502 Bad Gateway | 服务没起来 → `systemctl status sftpc-sgrna`；或端口被占 |
| 504 Gateway Timeout | `proxy_read_timeout` 没改 → 按上面改成 900s |
| 域名打不开但 IP 能开 | DNS 没生效，或 nginx 的 `server_name` 不匹配 |
| 查询报「滚动木」/距离缺失 | 输入的序列不在定位窗口内 → 填「序列基因组起点」或改「定位窗口」 |
| 结果为空 / 两站无交集 | 检查输入序列是不是基因组序列（cDNA 不行） |

---

## 文件

```
deploy/
    deploy.sh                  一键部署 / 更新脚本
    sftpc-sgrna.service        systemd 单元
    sgss.743.world.conf        nginx 反向代理配置（手工版）
    README.md                  本文件
```
