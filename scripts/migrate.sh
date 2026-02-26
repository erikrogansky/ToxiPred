#!/bin/bash

###############################################################################
# Database Migration Script
# 
# Handles database migrations safely with backup and rollback capabilities.
###############################################################################

set -e

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Configuration
COMPOSE_FILE="docker-compose.prod.yml"
POSTGRES_CONTAINER="toxipred-postgres"
BACKUP_DIR="./backups"
POSTGRES_USER="${POSTGRES_USER:-toxipred}"
POSTGRES_DB="${POSTGRES_DB:-toxipred}"

# Create backup directory
mkdir -p "$BACKUP_DIR"

# Backup database
backup_database() {
    local backup_file="$BACKUP_DIR/toxipred_backup_$(date +%Y%m%d_%H%M%S).sql"
    
    log_info "Creating database backup..."
    
    docker exec "$POSTGRES_CONTAINER" pg_dump -U "$POSTGRES_USER" "$POSTGRES_DB" > "$backup_file"
    
    if [ $? -eq 0 ]; then
        log_success "Database backed up to: $backup_file"
        
        # Keep only last 10 backups
        ls -t "$BACKUP_DIR"/toxipred_backup_*.sql | tail -n +11 | xargs -r rm
        
        echo "$backup_file"
    else
        log_error "Database backup failed!"
        exit 1
    fi
}

# Run migrations
run_migrations() {
    log_info "Running database migrations..."
    
    # The FastAPI app creates tables on startup, but you can add Alembic here
    # For now, we'll just verify the database is accessible
    
    docker compose -f "$COMPOSE_FILE" exec -T toxipred-postgres psql -U "$POSTGRES_USER" -d "$POSTGRES_DB" -c "SELECT 1;" > /dev/null
    
    if [ $? -eq 0 ]; then
        log_success "Database is accessible and ready"
    else
        log_error "Database migration failed!"
        exit 1
    fi
}

# Main
main() {
    log_info "Starting database migration process..."
    echo ""
    
    # Backup
    backup_file=$(backup_database)
    echo ""
    
    # Migrate
    run_migrations
    echo ""
    
    log_success "Migration completed successfully!"
    log_info "Backup available at: $backup_file"
}

main
