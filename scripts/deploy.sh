#!/bin/bash

###############################################################################
# ToxiPred Deployment Script
#
# Simple deploy: pull code, rebuild changed images, restart services.
###############################################################################

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

COMPOSE_FILE="docker-compose.prod.yml"
HEALTH_CHECK_RETRIES=30
HEALTH_CHECK_INTERVAL=2

log_info()    { echo -e "${BLUE}[INFO]${NC} $1"; }
log_success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
log_warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }
log_error()   { echo -e "${RED}[ERROR]${NC} $1"; }

check_health() {
    log_info "Running health checks..."
    for i in $(seq 1 $HEALTH_CHECK_RETRIES); do
        if docker exec server-toxipred micromamba run -n toxipred python -c "import urllib.request; urllib.request.urlopen('http://localhost:8000/health')" > /dev/null 2>&1; then
            log_success "Health check passed (attempt $i/$HEALTH_CHECK_RETRIES)"
            return 0
        fi
        log_warning "Health check pending (attempt $i/$HEALTH_CHECK_RETRIES)"
        sleep $HEALTH_CHECK_INTERVAL
    done
    log_error "Health check failed after $HEALTH_CHECK_RETRIES attempts"
    return 1
}

deploy() {
    echo ""
    log_info "====== ToxiPred Deployment Started ======"
    echo ""

    # Step 1: Start infrastructure (postgres, redis) if not running
    log_info "Step 1/4: Ensuring infrastructure is up..."
    docker compose -f "$COMPOSE_FILE" up -d toxipred-postgres toxipred-redis
    log_info "Waiting for database to be ready..."
    sleep 5
    log_success "Infrastructure is up"
    echo ""

    # Step 2: Run database migrations
    log_info "Step 2/4: Running database migrations..."
    bash scripts/migrate.sh
    log_success "Database migrations completed"
    echo ""

    # Step 3: Build and restart all services
    log_info "Step 3/4: Building and restarting services..."
    docker compose -f "$COMPOSE_FILE" build
    docker compose -f "$COMPOSE_FILE" up -d --no-build
    # Reload nginx so it resolves fresh upstream IPs after container restarts
    docker compose -f "$COMPOSE_FILE" restart nginx
    echo ""

    # Step 4: Health check
    log_info "Step 4/4: Verifying deployment..."
    sleep 5
    if ! check_health; then
        log_error "Deployment health check failed!"
        log_info "Check logs with: docker compose -f $COMPOSE_FILE logs server-toxipred"
        exit 1
    fi

    echo ""
    log_success "====== Deployment Completed Successfully ======"
    echo ""
}

status() {
    echo ""
    log_info "====== ToxiPred Status ======"
    echo ""

    log_info "Container status:"
    docker compose -f "$COMPOSE_FILE" ps
    echo ""

    echo -n "Backend: "
    if docker exec server-toxipred micromamba run -n toxipred python -c "import urllib.request; urllib.request.urlopen('http://localhost:8000/health')" > /dev/null 2>&1; then
        log_success "Healthy"
    else
        log_error "Unhealthy or stopped"
    fi
    echo ""
}

case "${1:-}" in
    deploy)
        deploy
        ;;
    status)
        status
        ;;
    *)
        echo "ToxiPred Deployment Script"
        echo ""
        echo "Usage: $0 {deploy|status}"
        echo ""
        echo "Commands:"
        echo "  deploy  - Build and restart all services"
        echo "  status  - Show current service status"
        echo ""
        exit 1
        ;;
esac
