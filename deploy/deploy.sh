#!/usr/bin/env bash
# 本项目由华南师范大学附属中学知识城校区高二一班同学制作
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# sgss.743.world 一键部署 / 更新脚本
# 在服务器上以 root 运行：  bash deploy.sh
#
# 先决条件：
#   1. DNS 已把 sgss.743.world 解析到本机
#   2. 宝塔面板里已建好站点 sgss.743.world（或手工放好 nginx 配置）
set -euo pipefail

DOMAIN="sgss.743.world"
PORT=1146                      # ⚠️ 1145 被 veridrop 占了
DIR="/home/ZZC-jeizuchin/sftpc-sgrna"
REPO="https://github.com/ZZC-jeizuchin/sftpc-gene-sgRNA-selection.git"
SERVICE="sftpc-sgrna"
NGINX_CONF="/www/server/panel/vhost/nginx/${DOMAIN}.conf"

say(){ printf '\033[36m==>\033[0m %s\n' "$*"; }
die(){ printf '\033[31m✗ %s\033[0m\n' "$*" >&2; exit 1; }

[ "$(id -u)" = "0" ] || die "需要 root"

# ---------------------------------------------------------------- 1. 端口检查
say "检查端口 $PORT 是否空闲"
if ss -ltn 2>/dev/null | grep -q ":${PORT} "; then
    OWNER=$(ss -ltnp 2>/dev/null | grep ":${PORT} " | grep -oP 'pid=\K[0-9]+' | head -1 || true)
    if [ -n "${OWNER:-}" ] && ps -p "$OWNER" -o cmd= 2>/dev/null | grep -q "run_web.py ${PORT}"; then
        say "端口已被本项目的旧进程占用（PID $OWNER），稍后会重启它"
    else
        die "端口 $PORT 已被别的服务占用：$(ps -p "${OWNER:-0}" -o cmd= 2>/dev/null || echo '未知')"
    fi
fi

# ---------------------------------------------------------------- 2. 拉代码
say "拉取 / 更新代码到 $DIR"
if [ -d "$DIR/.git" ]; then
    git -C "$DIR" fetch --depth 1 origin main
    git -C "$DIR" reset --hard origin/main
else
    rm -rf "$DIR"
    git clone --depth 1 "$REPO" "$DIR"
fi
say "代码版本：$(git -C "$DIR" rev-parse --short HEAD)  文件数：$(find "$DIR" -type f -not -path '*/.git/*' | wc -l)"

# ---------------------------------------------------------------- 3. systemd
say "安装 systemd 单元 /etc/systemd/system/${SERVICE}.service"
if [ -f "$DIR/deploy/sftpc-sgrna.service" ]; then
    cp "$DIR/deploy/sftpc-sgrna.service" "/etc/systemd/system/${SERVICE}.service"
else
    cat > "/etc/systemd/system/${SERVICE}.service" <<EOF
[Unit]
Description=SFTPC sgRNA Designer (CRISPOR + CHOPCHOP)
After=network-online.target
Wants=network-online.target

[Service]
Type=simple
User=root
WorkingDirectory=${DIR}
Environment=PYTHONUNBUFFERED=1
ExecStart=/usr/bin/python3 ${DIR}/run_web.py ${PORT}
Restart=on-failure
RestartSec=5
StandardOutput=journal
StandardError=journal
SyslogIdentifier=${SERVICE}
TimeoutStopSec=30

[Install]
WantedBy=multi-user.target
EOF
fi

# 清掉手工测试进程，交给 systemd 管
pkill -f "run_web.py ${PORT}" 2>/dev/null || true
sleep 1

systemctl daemon-reload
systemctl enable "$SERVICE" >/dev/null 2>&1 || true
systemctl restart "$SERVICE"
sleep 4

# ---------------------------------------------------------------- 4. 验证
say "服务状态"
systemctl is-active "$SERVICE" || die "服务没起来，看日志： journalctl -u ${SERVICE} -n 50"
ss -ltnp 2>/dev/null | grep ":${PORT} " | sed 's/^/    /' || die "端口 $PORT 没在监听"

say "本地冒烟测试"
HEALTH=$(curl -s -m 10 "http://127.0.0.1:${PORT}/api/health" || true)
[ "$HEALTH" = '{"ok": true, "service": "sgRNA Designer"}' ] \
    && say "  health ✓  $HEALTH" \
    || die "  health 返回异常：$HEALTH"
curl -s -m 10 -o /dev/null -w "    index  HTTP %{http_code}  %{size_download} bytes\n" "http://127.0.0.1:${PORT}/"

# ---------------------------------------------------------------- 5. nginx
if [ -f "$NGINX_CONF" ]; then
    say "nginx 配置已存在，测试并重载"
    nginx -t 2>&1 | sed 's/^/    /' && nginx -s reload && say "  已重载 ✓"
else
    say "⚠️  还没建站：$NGINX_CONF 不存在"
    say "   请在宝塔面板建站 $DOMAIN，然后："
    say "     站点设置 → 反向代理 → 目标 URL 填 http://127.0.0.1:${PORT}"
    say "     并把 proxy_read_timeout 改成 900s（长查询要用）"
    say "   或手工： cp $DIR/deploy/${DOMAIN}.conf $NGINX_CONF && nginx -t && nginx -s reload"
fi

say "完成 ✅"
echo
echo "  服务名 : $SERVICE"
echo "  端口   : 127.0.0.1:$PORT（只绑本地，不开防火墙）"
echo "  目录   : $DIR"
echo "  日志   : journalctl -u $SERVICE -f"
echo "  访问   : https://${DOMAIN}/"
