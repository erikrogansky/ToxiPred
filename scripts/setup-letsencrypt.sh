#!/usr/bin/env bash
set -euo pipefail

DOMAIN="${DOMAIN:-toxipred.fiit.stuba.sk}"
EMAIL="${LE_EMAIL:-}"
COMPOSE_FILE="${COMPOSE_FILE:-docker-compose.prod.yml}"
CERTBOT_CONF_DIR="${CERTBOT_CONF_DIR:-nginx/certbot/conf}"

create_temporary_certificate() {
  mkdir -p "$CERTBOT_CONF_DIR/live/$DOMAIN" nginx/certbot/www

  openssl req -x509 -nodes -newkey rsa:2048 -days 1 \
    -keyout "$CERTBOT_CONF_DIR/live/$DOMAIN/privkey.pem" \
    -out "$CERTBOT_CONF_DIR/live/$DOMAIN/fullchain.pem" \
    -subj "/CN=$DOMAIN"
}

restore_temporary_certificate_on_failure() {
  if [[ ! -f "$CERTBOT_CONF_DIR/live/$DOMAIN/fullchain.pem" || ! -f "$CERTBOT_CONF_DIR/live/$DOMAIN/privkey.pem" ]]; then
    echo "Certificate request failed before a real certificate was installed."
    echo "Restoring a temporary self-signed certificate so nginx can keep starting."
    create_temporary_certificate
    docker compose -f "$COMPOSE_FILE" up -d nginx || true
  fi
}

trap restore_temporary_certificate_on_failure ERR

if [[ -z "$EMAIL" ]]; then
  echo "Set LE_EMAIL before running, for example:"
  echo "  LE_EMAIL=you@example.com $0"
  exit 1
fi

if ! command -v docker >/dev/null 2>&1; then
  echo "docker is required"
  exit 1
fi

if [[ ! -f "$CERTBOT_CONF_DIR/live/$DOMAIN/fullchain.pem" || ! -f "$CERTBOT_CONF_DIR/live/$DOMAIN/privkey.pem" ]]; then
  echo "Creating temporary self-signed certificate for nginx bootstrap..."
  create_temporary_certificate
fi

echo "Starting nginx with the temporary certificate..."
docker compose -f "$COMPOSE_FILE" up -d nginx

echo "Removing temporary certificate files before requesting a real Let's Encrypt certificate..."
rm -rf "$CERTBOT_CONF_DIR/live/$DOMAIN"
rm -rf "$CERTBOT_CONF_DIR/archive/$DOMAIN"
rm -f "$CERTBOT_CONF_DIR/renewal/$DOMAIN.conf"

echo "Requesting Let's Encrypt certificate for $DOMAIN..."
docker compose -f "$COMPOSE_FILE" --profile certbot run --rm certbot certonly \
  --webroot \
  --webroot-path /var/www/certbot \
  --email "$EMAIL" \
  --agree-tos \
  --no-eff-email \
  --force-renewal \
  -d "$DOMAIN"

echo "Reloading nginx with the real certificate..."
trap - ERR
docker compose -f "$COMPOSE_FILE" exec nginx nginx -t
docker compose -f "$COMPOSE_FILE" exec nginx nginx -s reload

echo "Done. Verify with:"
echo "  curl -I https://$DOMAIN/"
