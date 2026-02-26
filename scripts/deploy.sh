#!/bin/bash

###############################################################################
# ToxiPred Blue-Green Deployment Script
# 
# This script performs zero-downtime deployments using blue-green strategy.
# It builds new containers, runs health checks, switches traffic, and 
# provides automatic rollback on failure.
###############################################################################

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Configuration
COMPOSE_FILE="docker-compose.prod.yml"
NGINX_UPSTREAM_FILE="nginx/nginx.conf"
NGINX_CONTAINER="toxipred-nginx"
HEALTH_CHECK_RETRIES=30
HEALTH_CHECK_INTERVAL=2

# Logging functions
log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Determine current active deployment
get_active_deployment() {
    if grep -q "server server-toxipred-blue:8000" "$NGINX_UPSTREAM_FILE" && \
       ! grep -q "# server server-toxipred-blue:8000" "$NGINX_UPSTREAM_FILE"; then
        echo "blue"
    else
        echo "green"
    fi
}

# Get inactive deployment (target for new deployment)
get_inactive_deployment() {
    local active=$(get_active_deployment)
    if [ "$active" == "blue" ]; then
        echo "green"
    else
        echo "blue"
    fi
}

# Health check function
check_health() {
    local deployment=$1
    local container="server-toxipred-$deployment"
    
    log_info "Running health checks for $deployment deployment..."
    
    for i in $(seq 1 $HEALTH_CHECK_RETRIES); do
        if docker exec "$container" curl -sf http://localhost:8000/health > /dev/null 2>&1; then
            log_success "Health check passed for $deployment (attempt $i/$HEALTH_CHECK_RETRIES)"
            return 0
        fi
        log_warning "Health check failed for $deployment (attempt $i/$HEALTH_CHECK_RETRIES)"
        sleep $HEALTH_CHECK_INTERVAL
    done
    
    log_error "Health check failed for $deployment after $HEALTH_CHECK_RETRIES attempts"
    return 1
}

# Switch nginx upstream
switch_traffic() {
    local target=$1
    local backup_file="${NGINX_UPSTREAM_FILE}.backup.$(date +%s)"
    
    log_info "Switching traffic to $target deployment..."
    
    # Backup current config
    cp "$NGINX_UPSTREAM_FILE" "$backup_file"
    
    if [ "$target" == "blue" ]; then
        sed -i.tmp 's/# server server-toxipred-blue:8000/server server-toxipred-blue:8000/' "$NGINX_UPSTREAM_FILE"
        sed -i.tmp 's/^[[:space:]]*server server-toxipred-green:8000/        # server server-toxipred-green:8000/' "$NGINX_UPSTREAM_FILE"
    else
        sed -i.tmp 's/# server server-toxipred-green:8000/server server-toxipred-green:8000/' "$NGINX_UPSTREAM_FILE"
        sed -i.tmp 's/^[[:space:]]*server server-toxipred-blue:8000/        # server server-toxipred-blue:8000/' "$NGINX_UPSTREAM_FILE"
    fi
    
    rm -f "${NGINX_UPSTREAM_FILE}.tmp"
    
    # Reload nginx
    docker exec "$NGINX_CONTAINER" nginx -t && \
    docker exec "$NGINX_CONTAINER" nginx -s reload
    
    log_success "Traffic switched to $target deployment"
    log_info "Backup saved to: $backup_file"
}

# Rollback function
rollback() {
    local from_deployment=$1
    local to_deployment=$2
    
    log_warning "Rolling back from $from_deployment to $to_deployment..."
    
    switch_traffic "$to_deployment"
    
    log_info "Stopping failed $from_deployment deployment..."
    if [ "$from_deployment" == "green" ]; then
        docker compose -f "$COMPOSE_FILE" --profile green stop server-toxipred-green
    else
        docker compose -f "$COMPOSE_FILE" stop server-toxipred-blue
    fi
    
    log_success "Rollback completed successfully"
}

# Main deployment function
deploy() {
    echo ""
    log_info "====== ToxiPred Deployment Started ======"
    echo ""
    
    # Get deployment slots
    local active=$(get_active_deployment)
    local target=$(get_inactive_deployment)
    
    log_info "Active deployment: $active"
    log_info "Target deployment: $target"
    echo ""
    
    # Step 1: Run database migrations
    log_info "Step 1/6: Running database migrations..."
    bash scripts/migrate.sh
    log_success "Database migrations completed"
    echo ""
    
    # Step 2: Build new images
    log_info "Step 2/6: Building new Docker images..."
    docker compose -f "$COMPOSE_FILE" build server-toxipred-$target app-toxipred
    log_success "Docker images built successfully"
    echo ""
    
    # Step 3: Start new deployment
    log_info "Step 3/6: Starting $target deployment..."
    if [ "$target" == "green" ]; then
        docker compose -f "$COMPOSE_FILE" --profile green up -d server-toxipred-green
    else
        docker compose -f "$COMPOSE_FILE" up -d server-toxipred-blue
    fi
    log_success "$target deployment started"
    echo ""
    
    # Step 4: Wait for containers to be ready
    log_info "Step 4/6: Waiting for $target deployment to be ready..."
    sleep 10  # Initial wait for container startup
    
    if ! check_health "$target"; then
        log_error "Health check failed for $target deployment"
        rollback "$target" "$active"
        exit 1
    fi
    echo ""
    
    # Step 5: Switch traffic
    log_info "Step 5/6: Switching traffic to $target deployment..."
    switch_traffic "$target"
    
    # Verify traffic switch with another health check
    sleep 5
    if ! check_health "$target"; then
        log_error "Health check failed after traffic switch"
        rollback "$target" "$active"
        exit 1
    fi
    log_success "Traffic successfully switched to $target"
    echo ""
    
    # Step 6: Stop old deployment
    log_info "Step 6/6: Stopping old $active deployment..."
    sleep 10  # Grace period for in-flight requests
    
    if [ "$active" == "green" ]; then
        docker compose -f "$COMPOSE_FILE" --profile green stop server-toxipred-green
    else
        docker compose -f "$COMPOSE_FILE" stop server-toxipred-blue
    fi
    log_success "Old $active deployment stopped"
    echo ""
    
    # Update frontend (rebuilds with production Quasar bundle)
    log_info "Updating frontend (building production Quasar bundle)..."
    log_info "This may take 1-2 minutes for the build..."
    docker compose -f "$COMPOSE_FILE" build app-toxipred
    docker compose -f "$COMPOSE_FILE" up -d app-toxipred
    log_success "Frontend updated (brief ~5s interruption during swap)"
    echo ""
    
    log_success "====== Deployment Completed Successfully ======"
    log_info "Active deployment: $target"
    log_info "Previous deployment ($active) is stopped but can be rolled back if needed"
    echo ""
}

# Quick rollback to previous deployment
quick_rollback() {
    local active=$(get_active_deployment)
    local previous=$(get_inactive_deployment)
    
    echo ""
    log_warning "====== Quick Rollback Initiated ======"
    echo ""
    
    log_info "Current active: $active"
    log_info "Rolling back to: $previous"
    echo ""
    
    # Start previous deployment
    log_info "Starting $previous deployment..."
    if [ "$previous" == "green" ]; then
        docker compose -f "$COMPOSE_FILE" --profile green up -d server-toxipred-green
    else
        docker compose -f "$COMPOSE_FILE" up -d server-toxipred-blue
    fi
    
    sleep 10
    
    if ! check_health "$previous"; then
        log_error "Previous deployment health check failed! Manual intervention required."
        exit 1
    fi
    
    switch_traffic "$previous"
    
    log_success "====== Rollback Completed ======"
    echo ""
}

# Status check
status() {
    local active=$(get_active_deployment)
    
    echo ""
    log_info "====== ToxiPred Deployment Status ======"
    echo ""
    log_info "Active deployment: $active"
    echo ""
    
    log_info "Container health status:"
    docker compose -f "$COMPOSE_FILE" ps
    echo ""
    
    log_info "Backend health checks:"
    echo -n "Blue: "
    if docker exec server-toxipred-blue curl -sf http://localhost:8000/health > /dev/null 2>&1; then
        log_success "Healthy"
    else
        log_error "Unhealthy or stopped"
    fi
    
    echo -n "Green: "
    if docker exec server-toxipred-green curl -sf http://localhost:8000/health > /dev/null 2>&1; then
        log_success "Healthy"
    else
        log_warning "Unhealthy or stopped (this is normal if not active)"
    fi
    echo ""
}

# Command handling
case "${1:-}" in
    deploy)
        deploy
        ;;
    rollback)
        quick_rollback
        ;;
    status)
        status
        ;;
    *)
        echo "ToxiPred Deployment Script"
        echo ""
        echo "Usage: $0 {deploy|rollback|status}"
        echo ""
        echo "Commands:"
        echo "  deploy    - Deploy new version with zero downtime"
        echo "  rollback  - Quickly rollback to previous deployment"
        echo "  status    - Show current deployment status"
        echo ""
        exit 1
        ;;
esac
