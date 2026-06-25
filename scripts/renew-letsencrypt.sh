#!/usr/bin/env bash
set -euo pipefail

COMPOSE_FILE="${COMPOSE_FILE:-docker-compose.prod.yml}"

docker compose -f "$COMPOSE_FILE" --profile certbot run --rm certbot renew \
  --webroot \
  --webroot-path /var/www/certbot

docker compose -f "$COMPOSE_FILE" exec nginx nginx -t
docker compose -f "$COMPOSE_FILE" exec nginx nginx -s reload
