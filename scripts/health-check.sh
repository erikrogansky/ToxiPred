#!/bin/bash

###############################################################################
# System Health Check Script
# 
# Comprehensive health check for ToxiPred deployment
###############################################################################

set -e

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

COMPOSE_FILE="docker-compose.prod.yml"

echo ""
echo "======================================"
echo "  ToxiPred System Health Check"
echo "======================================"
echo ""

# Check disk space
echo -e "${BLUE}[DISK SPACE]${NC}"
df -h / | tail -1 | awk '{
    usage = substr($5, 1, length($5)-1);
    if (usage > 80) {
        printf "\033[0;31m⚠ WARNING: Disk usage at %s%%\033[0m\n", usage;
        exit 1;
    } else {
        printf "\033[0;32m✓ Disk usage: %s%%\033[0m\n", usage;
    }
}'
echo ""

# Check Docker disk
echo -e "${BLUE}[DOCKER DISK]${NC}"
docker system df | tail -n +2
echo ""

# Check containers
echo -e "${BLUE}[CONTAINERS]${NC}"
docker compose -f "$COMPOSE_FILE" ps --format "table {{.Name}}\t{{.Status}}\t{{.Health}}"
echo ""

# Check backend health
echo -e "${BLUE}[BACKEND HEALTH]${NC}"
for deployment in blue green; do
    container="server-toxipred-$deployment"
    echo -n "$deployment: "
    if docker ps --format '{{.Names}}' | grep -q "^${container}$"; then
        if docker exec "$container" curl -sf http://localhost:8000/health > /dev/null 2>&1; then
            echo -e "${GREEN}✓ Healthy${NC}"
        else
            echo -e "${RED}✗ Unhealthy${NC}"
        fi
    else
        echo -e "${YELLOW}○ Stopped${NC}"
    fi
done
echo ""

# Check database
echo -e "${BLUE}[DATABASE]${NC}"
if docker exec toxipred-postgres pg_isready -U toxipred > /dev/null 2>&1; then
    echo -e "${GREEN}✓ PostgreSQL is ready${NC}"
    
    # Check database size
    db_size=$(docker exec toxipred-postgres psql -U toxipred -d toxipred -t -c "SELECT pg_size_pretty(pg_database_size('toxipred'));")
    echo "  Database size: $db_size"
else
    echo -e "${RED}✗ PostgreSQL is not ready${NC}"
fi
echo ""

# Check Redis
echo -e "${BLUE}[REDIS]${NC}"
if docker exec toxipred-redis redis-cli ping > /dev/null 2>&1; then
    echo -e "${GREEN}✓ Redis is ready${NC}"
    
    # Check memory
    redis_mem=$(docker exec toxipred-redis redis-cli info memory | grep used_memory_human | cut -d: -f2)
    echo "  Memory used: $redis_mem"
else
    echo -e "${RED}✗ Redis is not ready${NC}"
fi
echo ""

# Check workers
echo -e "${BLUE}[CELERY WORKERS]${NC}"
if docker exec toxipred-worker bash -c "micromamba run -n toxipred celery -A app.workers.celery_app:celery_app inspect ping" > /dev/null 2>&1; then
    echo -e "${GREEN}✓ Workers are responding${NC}"
    
    # Active tasks
    active=$(docker exec toxipred-worker bash -c "micromamba run -n toxipred celery -A app.workers.celery_app:celery_app inspect active" 2>/dev/null | grep -c "task" || echo "0")
    echo "  Active tasks: $active"
else
    echo -e "${RED}✗ Workers are not responding${NC}"
fi
echo ""

# Check recent logs for errors
echo -e "${BLUE}[RECENT ERRORS]${NC}"
error_count=$(docker compose -f "$COMPOSE_FILE" logs --since 1h 2>/dev/null | grep -i error | wc -l)
if [ "$error_count" -gt 0 ]; then
    echo -e "${YELLOW}⚠ Found $error_count errors in last hour${NC}"
    echo "  Run 'make prod-logs' to investigate"
else
    echo -e "${GREEN}✓ No errors in last hour${NC}"
fi
echo ""

# Check backups
echo -e "${BLUE}[BACKUPS]${NC}"
if [ -d "backups" ]; then
    backup_count=$(ls -1 backups/toxipred_backup_*.sql 2>/dev/null | wc -l)
    if [ "$backup_count" -gt 0 ]; then
        latest_backup=$(ls -t backups/toxipred_backup_*.sql 2>/dev/null | head -1)
        backup_age=$(( ($(date +%s) - $(stat -f %m "$latest_backup" 2>/dev/null || stat -c %Y "$latest_backup")) / 86400 ))
        echo -e "${GREEN}✓ $backup_count backups found${NC}"
        echo "  Latest backup: $(basename "$latest_backup") ($backup_age days old)"
    else
        echo -e "${YELLOW}⚠ No backups found${NC}"
    fi
else
    echo -e "${YELLOW}⚠ Backup directory not found${NC}"
fi
echo ""

# Resource usage
echo -e "${BLUE}[RESOURCE USAGE]${NC}"
docker stats --no-stream --format "table {{.Name}}\t{{.CPUPerc}}\t{{.MemUsage}}" | grep toxipred
echo ""

echo "======================================"
echo "  Health Check Complete"
echo "======================================"
echo ""
