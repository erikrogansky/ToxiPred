# 🔧 Advanced Deployment Guide

**Detailed technical information, advanced troubleshooting, and customization options.**

For basic setup, see: `SETUP.md`

---

## 📐 Architecture Overview

### System Components

```
┌─────────────────────────────────────────────────────────────────┐
│                        Internet / Users                          │
└────────────────────────────┬────────────────────────────────────┘
                             │
                             ▼
                    ┌────────────────┐
                    │  Nginx (Port 80)│
                    │  Load Balancer  │
                    └────────┬───────┘
                             │
           ┌─────────────────┼─────────────────┐
           │                 │                 │
           ▼                 ▼                 ▼
    ┌───────────┐    ┌───────────┐    ┌───────────┐
    │  Backend  │    │  Backend  │    │ Frontend  │
    │   BLUE    │    │   GREEN   │    │  (Quasar) │
    │  (Active) │    │ (Standby) │    │           │
    └─────┬─────┘    └─────┬─────┘    └───────────┘
          │                │
          └────────┬───────┘
                   │
    ┌──────────────┼──────────────────────────┐
    │              │                          │
    ▼              ▼                          ▼
┌──────────┐  ┌──────────┐            ┌──────────┐
│PostgreSQL│  │  Redis   │            │  Celery  │
│Database  │  │  Queue   │            │  Workers │
└──────────┘  └──────────┘            └──────────┘
     │             │                         │
     ▼             ▼                         ▼
┌──────────┐  ┌──────────┐            ┌──────────┐
│ Volume:  │  │ Volume:  │            │  Flower  │
│ pgdata   │  │  redis   │            │Monitor UI│
└──────────┘  └──────────┘            └──────────┘
```

### Network Configuration

All services run on an internal Docker network (`toxipred-network`):
- Only Nginx exposes ports to the host
- Backend services communicate via internal DNS
- Persistent volumes for data safety

**Exposed Ports:**
- `80` - HTTP (Nginx)
- `443` - HTTPS (Nginx, if configured)
- `5555` - Flower monitoring UI

---

## 🔄 Blue-Green Deployment Deep Dive

### How It Works

The system maintains two identical backend environments (Blue and Green). At any time, one is active and serving traffic, while the other is standby.

**Deployment Flow:**

```
State 1: Blue Active
├─ Nginx → Blue:8000 (serving users)
└─ Green (stopped)

State 2: Deploy to Green
├─ Nginx → Blue:8000 (still serving users)
├─ Green starting...
├─ Green health checks running...
└─ Green ready ✓

State 3: Traffic Switch
├─ Nginx configuration updated
├─ Nginx reload (0 downtime)
├─ Nginx → Green:8000 (now serving users)
└─ Blue draining connections...

State 4: Cleanup
├─ Nginx → Green:8000 (serving users)
├─ Blue stopped
└─ Ready for next deployment (will use Blue)
```

### Health Check Mechanism

Before switching traffic, the deployment script verifies:

1. **Container Running:** `docker ps | grep server-toxipred-green`
2. **HTTP Endpoint:** `curl http://localhost:8000/health`
3. **Database Connectivity:** `curl http://localhost:8000/ready`
4. **Retry Logic:** 30 attempts, 2-second intervals

If ANY check fails → automatic rollback.

### Nginx Traffic Switching

The deployment script modifies `nginx/nginx.conf`:

```nginx
# Before (Blue active)
upstream backend {
    server server-toxipred-blue:8000;
    # server server-toxipred-green:8000;
}

# After (Green active)
upstream backend {
    # server server-toxipred-blue:8000;
    server server-toxipred-green:8000;
}
```

Then executes: `nginx -t && nginx -s reload` (zero downtime reload).

---

## 🗄️ Database Management

### Automatic Backups

Every deployment triggers:
```bash
docker exec toxipred-postgres pg_dump -U toxipred toxipred > \
  backups/toxipred_backup_$(date +%Y%m%d_%H%M%S).sql
```

**Retention:** Last 10 backups kept automatically.

### Manual Backup & Restore

**Backup:**
```bash
# Full database dump
docker exec toxipred-postgres pg_dump -U toxipred toxipred > backup.sql

# Compressed backup
docker exec toxipred-postgres pg_dump -U toxipred toxipred | gzip > backup.sql.gz

# Schema only
docker exec toxipred-postgres pg_dump -U toxipred --schema-only toxipred > schema.sql

# Data only
docker exec toxipred-postgres pg_dump -U toxipred --data-only toxipred > data.sql
```

**Restore:**
```bash
# Stop backend to prevent conflicts
docker compose -f docker-compose.prod.yml stop server-toxipred-blue server-toxipred-green

# Restore database
docker exec -i toxipred-postgres psql -U toxipred toxipred < backup.sql

# Or from compressed
gunzip -c backup.sql.gz | docker exec -i toxipred-postgres psql -U toxipred toxipred

# Restart backend
make prod-up
```

### Database Migrations

Currently using SQLAlchemy's `Base.metadata.create_all()` on startup.

**For production-grade migrations (optional):**

Install Alembic:
```bash
cd server-toxipred
pip install alembic
alembic init alembic
```

Configure `alembic/env.py` to use your models, then:
```bash
# Create migration
alembic revision --autogenerate -m "Add new column"

# Apply migration
alembic upgrade head

# Rollback
alembic downgrade -1
```

---

## 🔧 Advanced Configuration

### Environment Variables Reference

**Database:**
- `POSTGRES_USER` - Database username (default: `toxipred`)
- `POSTGRES_PASSWORD` - Database password (⚠️ change in production!)
- `POSTGRES_DB` - Database name (default: `toxipred`)
- `DATABASE_URL` - Full connection string (auto-generated)

**Celery/Redis:**
- `CELERY_BROKER_URL` - Redis URL for task queue
- `CELERY_RESULT_BACKEND` - Redis URL for results
- `CELERY_DEFAULT_QUEUE` - Queue name (default: `default`)

**Application:**
- `PUBCHEM_LOOKUP` - Enable PubChem lookups (1=on, 0=off)
- `API_URL` - Frontend API endpoint URL

### Nginx Configuration

**Rate Limiting:**

Edit `nginx/nginx.conf`:
```nginx
# Current: 10 requests/second for API
limit_req_zone $binary_remote_addr zone=api_limit:10m rate=10r/s;

# Increase to 30/s:
limit_req_zone $binary_remote_addr zone=api_limit:10m rate=30r/s;
```

**HTTPS/SSL:**

1. Place certificates in `nginx/ssl/`:
   - `cert.pem`
   - `key.pem`

2. Edit `nginx/conf.d/toxipred.conf`:
```nginx
server {
    listen 443 ssl;
    server_name your-domain.com;

    ssl_certificate /etc/nginx/ssl/cert.pem;
    ssl_certificate_key /etc/nginx/ssl/key.pem;
    
    # Rest of configuration...
}

# Redirect HTTP to HTTPS
server {
    listen 80;
    server_name your-domain.com;
    return 301 https://$server_name$request_uri;
}
```

3. Reload: `docker exec toxipred-nginx nginx -s reload`

**Custom Headers:**

Add to `nginx/conf.d/toxipred.conf`:
```nginx
location / {
    add_header X-Custom-Header "value" always;
    # Rest of configuration...
}
```

### Scaling Workers

**Horizontal Scaling:**
```bash
# Scale to 5 workers
docker compose -f docker-compose.prod.yml up -d --scale toxipred-worker=5

# Scale back to 1
docker compose -f docker-compose.prod.yml up -d --scale toxipred-worker=1
```

**Worker Configuration:**

Edit `server-toxipred/app/workers/celery_app.py`:
```python
celery_app.conf.update(
    worker_prefetch_multiplier=4,  # Tasks per worker
    worker_max_tasks_per_child=1000,  # Restart after N tasks
    task_time_limit=300,  # Task timeout (5 minutes)
    task_soft_time_limit=270,  # Soft timeout (4.5 minutes)
)
```

### Load Balancing Multiple Backends

To run both Blue and Green simultaneously (2x capacity):

Edit `nginx/nginx.conf`:
```nginx
upstream backend {
    least_conn;  # Distribute based on connections
    server server-toxipred-blue:8000 max_fails=3 fail_timeout=30s;
    server server-toxipred-green:8000 max_fails=3 fail_timeout=30s;
}
```

Start both:
```bash
docker compose -f docker-compose.prod.yml up -d server-toxipred-blue
docker compose -f docker-compose.prod.yml --profile green up -d server-toxipred-green
```

---

## 🐛 Advanced Troubleshooting

### Container Logs Analysis

**Follow logs with filters:**
```bash
# Only errors
docker compose -f docker-compose.prod.yml logs | grep -i error

# Specific timeframe (last hour)
docker compose -f docker-compose.prod.yml logs --since 1h

# JSON format for parsing
docker compose -f docker-compose.prod.yml logs --json

# Export logs
docker compose -f docker-compose.prod.yml logs > full_logs.txt
```

### Performance Debugging

**Check resource usage:**
```bash
# Container stats
docker stats --no-stream

# Specific container
docker stats server-toxipred-blue --no-stream

# Continuous monitoring
watch -n 2 'docker stats --no-stream'
```

**Database query performance:**
```bash
# Connect to database
docker exec -it toxipred-postgres psql -U toxipred toxipred

# Show slow queries
SELECT * FROM pg_stat_statements ORDER BY mean_time DESC LIMIT 10;

# Show active connections
SELECT * FROM pg_stat_activity;

# Kill a query
SELECT pg_terminate_backend(pid) FROM pg_stat_activity WHERE pid = 12345;
```

**Redis monitoring:**
```bash
# Connect to Redis
docker exec -it toxipred-redis redis-cli

# Show stats
INFO

# Monitor commands
MONITOR

# Check queue length
LLEN celery

# Show memory usage
MEMORY STATS
```

### Network Issues

**Check container connectivity:**
```bash
# Test DNS resolution
docker exec server-toxipred-blue ping toxipred-postgres

# Test port connectivity
docker exec server-toxipred-blue nc -zv toxipred-postgres 5432

# Check network
docker network inspect toxipred-network
```

**Nginx debugging:**
```bash
# Test configuration
docker exec toxipred-nginx nginx -t

# Check upstream status
docker exec toxipred-nginx cat /etc/nginx/nginx.conf

# Access logs
docker exec toxipred-nginx tail -f /var/log/nginx/access.log

# Error logs
docker exec toxipred-nginx tail -f /var/log/nginx/error.log
```

### Deployment Script Debugging

Enable verbose output:
```bash
# Edit scripts/deploy.sh
set -x  # Add after 'set -e'

# Run deployment
./scripts/deploy.sh deploy 2>&1 | tee deployment.log
```

Check deployment state:
```bash
# Current active deployment
grep "server server-toxipred-" nginx/nginx.conf

# Last deployment backup
ls -lt backups/ | head -5

# Container uptime
docker ps --format "table {{.Names}}\t{{.Status}}"
```

### Recovery Scenarios

**Scenario 1: Both Blue and Green Down**

```bash
# Check what's running
docker compose -f docker-compose.prod.yml ps

# Restart infrastructure
docker compose -f docker-compose.prod.yml restart toxipred-postgres toxipred-redis

# Start backend (whichever was active)
docker compose -f docker-compose.prod.yml up -d server-toxipred-blue

# Verify
make health
```

**Scenario 2: Database Corrupted**

```bash
# Stop backends
docker compose -f docker-compose.prod.yml stop server-toxipred-blue server-toxipred-green

# Backup current (corrupted) database
docker exec toxipred-postgres pg_dump -U toxipred toxipred > corrupted_backup.sql

# Drop and recreate database
docker exec toxipred-postgres psql -U toxipred -c "DROP DATABASE toxipred;"
docker exec toxipred-postgres psql -U toxipred -c "CREATE DATABASE toxipred;"

# Restore from backup
docker exec -i toxipred-postgres psql -U toxipred toxipred < backups/latest_good_backup.sql

# Restart
make prod-up
make health
```

**Scenario 3: Disk Full**

```bash
# Check usage
df -h
docker system df

# Emergency cleanup
docker system prune -a -f
docker volume prune -f

# Clean old backups
find backups/ -name "*.sql" -mtime +30 -delete

# Clean logs
sudo truncate -s 0 /var/lib/docker/containers/*/*-json.log
```

**Scenario 4: Nginx Not Routing**

```bash
# Check Nginx is running
docker ps | grep nginx

# Test config
docker exec toxipred-nginx nginx -t

# Check upstream config
docker exec toxipred-nginx cat /etc/nginx/nginx.conf | grep "server server-toxipred"

# Reload config
docker exec toxipred-nginx nginx -s reload

# Restart if needed
docker compose -f docker-compose.prod.yml restart nginx
```

---

## 🔐 Security Hardening

### Server-Level Security

**Firewall (UFW):**
```bash
sudo ufw default deny incoming
sudo ufw default allow outgoing
sudo ufw allow 22/tcp    # SSH
sudo ufw allow 80/tcp    # HTTP
sudo ufw allow 443/tcp   # HTTPS
sudo ufw enable
```

**SSH Hardening:**

Edit `/etc/ssh/sshd_config`:
```
PermitRootLogin no
PasswordAuthentication no
PubkeyAuthentication yes
Port 22  # Consider changing
```

Restart SSH: `sudo systemctl restart sshd`

**Fail2Ban:**
```bash
sudo apt install fail2ban
sudo systemctl enable fail2ban
sudo systemctl start fail2ban
```

### Application-Level Security

**Secure Flower (Celery Monitor):**

Edit `nginx/conf.d/toxipred.conf`:
```nginx
server {
    listen 5555;
    location / {
        auth_basic "Flower Monitoring";
        auth_basic_user_file /etc/nginx/.htpasswd;
        proxy_pass http://flower/;
    }
}
```

Create password file:
```bash
# Install htpasswd
sudo apt install apache2-utils

# Create password
htpasswd -c nginx/.htpasswd admin

# Rebuild nginx
docker compose -f docker-compose.prod.yml restart nginx
```

**Database Access:**

Restrict PostgreSQL to internal network only:
```yaml
# docker-compose.prod.yml
toxipred-postgres:
  # Remove ports section to prevent external access
  # ports:
  #   - "5432:5432"
```

**Environment Secrets:**

Never commit `.env`:
```bash
# Verify it's ignored
git status .env  # Should not appear

# If accidentally committed
git rm --cached .env
git commit -m "Remove .env from tracking"
```

---

## 📊 Monitoring & Observability

### Health Check Endpoints

**Backend Health:**
- `GET /health` - Basic health check
- `GET /ready` - Readiness check (includes DB connectivity)

**Response format:**
```json
{
  "status": "healthy",
  "service": "toxipred-api"
}
```

### Custom Health Monitoring

Create a monitoring script:
```bash
#!/bin/bash
# monitoring.sh

ENDPOINTS=(
    "http://your-server.com/api/health"
    "http://your-server.com/api/ready"
)

for endpoint in "${ENDPOINTS[@]}"; do
    if curl -sf "$endpoint" > /dev/null; then
        echo "✓ $endpoint OK"
    else
        echo "✗ $endpoint FAILED"
        # Send alert (email, Slack, etc.)
    fi
done
```

Run via cron every 5 minutes:
```bash
*/5 * * * * /path/to/monitoring.sh >> /var/log/health-monitor.log 2>&1
```

### Centralized Logging (Optional)

**Using Docker logging driver:**

Edit `docker-compose.prod.yml`:
```yaml
services:
  server-toxipred-blue:
    logging:
      driver: "json-file"
      options:
        max-size: "10m"
        max-file: "3"
```

**Export logs to external service:**
```bash
# Example: Send to Papertrail
docker logs -f server-toxipred-blue | \
  ncat --ssl logs.papertrailapp.com 12345
```

---

## 🚀 Performance Optimization

### Database Optimization

**Connection pooling:**

Edit `server-toxipred/app/db.py`:
```python
engine = create_engine(
    DATABASE_URL,
    pool_size=20,           # Max connections
    max_overflow=10,        # Extra connections
    pool_pre_ping=True,     # Verify connections
    pool_recycle=3600,      # Recycle after 1 hour
)
```

**Indexes:**
```sql
-- Connect to database
docker exec -it toxipred-postgres psql -U toxipred toxipred

-- Add indexes for common queries
CREATE INDEX idx_jobs_created_at ON jobs(created_at);
CREATE INDEX idx_jobs_status ON jobs(status);
```

### Redis Optimization

**Increase memory limit:**

Edit `docker-compose.prod.yml`:
```yaml
toxipred-redis:
  command: redis-server --maxmemory 512mb --maxmemory-policy allkeys-lru
```

### Nginx Caching

**Static asset caching:**

Edit `nginx/conf.d/toxipred.conf`:
```nginx
location ~* \.(jpg|jpeg|png|gif|ico|css|js)$ {
    expires 1y;
    add_header Cache-Control "public, immutable";
}
```

**API response caching:**
```nginx
proxy_cache_path /var/cache/nginx levels=1:2 keys_zone=api_cache:10m;

location /api/ {
    proxy_cache api_cache;
    proxy_cache_valid 200 5m;
    proxy_cache_key "$request_uri";
}
```

---

## 🔄 CI/CD Customization

### Custom Deploy Conditions

Edit `.github/workflows/deploy.yml`:

**Deploy only on tags:**
```yaml
on:
  push:
    tags:
      - 'v*'
```

**Deploy on schedule:**
```yaml
on:
  schedule:
    - cron: '0 2 * * *'  # Daily at 2 AM
```

**Deploy to staging first:**
```yaml
jobs:
  deploy-staging:
    runs-on: ubuntu-latest
    steps:
      # Deploy to staging server
      
  deploy-production:
    needs: deploy-staging
    runs-on: ubuntu-latest
    steps:
      # Deploy to production after staging succeeds
```

### Notifications

**Slack notifications:**
```yaml
- name: Notify Slack
  if: always()
  uses: slackapi/slack-github-action@v1.24.0
  with:
    webhook-url: ${{ secrets.SLACK_WEBHOOK }}
    payload: |
      {
        "text": "Deployment ${{ job.status }}: ${{ github.event.head_commit.message }}"
      }
```

**Discord notifications:**
```yaml
- name: Notify Discord
  if: always()
  run: |
    curl -X POST "${{ secrets.DISCORD_WEBHOOK }}" \
      -H "Content-Type: application/json" \
      -d '{"content": "Deployment ${{ job.status }}"}'
```

---

## 📝 Maintenance Tasks

### Regular Maintenance Checklist

**Weekly:**
- [ ] Check disk space: `df -h`
- [ ] Review logs for errors: `make prod-logs`
- [ ] Verify backups exist: `ls -lh backups/`
- [ ] Check application health: `make health`

**Monthly:**
- [ ] Clean old Docker images: `make clean`
- [ ] Archive old backups
- [ ] Review and optimize database
- [ ] Check security updates: `sudo apt update && sudo apt list --upgradable`

**Quarterly:**
- [ ] Review and update dependencies
- [ ] Performance audit
- [ ] Security audit
- [ ] Disaster recovery test

### Automated Maintenance

**Auto-cleanup script:**
```bash
#!/bin/bash
# cleanup.sh

# Remove old Docker images
docker system prune -a -f

# Remove old backups (keep last 30)
ls -t backups/*.sql | tail -n +31 | xargs -r rm

# Remove old logs (older than 30 days)
find /var/log/toxipred -name "*.log" -mtime +30 -delete

echo "Cleanup completed: $(date)"
```

Add to cron:
```bash
0 3 * * 0 /path/to/cleanup.sh >> /var/log/cleanup.log 2>&1
```

---

## 📚 Additional Resources

### Useful Docker Commands

```bash
# View container resource usage
docker stats

# Inspect container
docker inspect server-toxipred-blue

# Execute command in container
docker exec -it server-toxipred-blue bash

# Copy files from container
docker cp server-toxipred-blue:/app/file.txt ./

# View container filesystem
docker exec server-toxipred-blue ls -la /app

# Check container IP
docker inspect server-toxipred-blue | grep IPAddress
```

### PostgreSQL Commands

```bash
# Database size
docker exec toxipred-postgres psql -U toxipred -c "SELECT pg_size_pretty(pg_database_size('toxipred'));"

# Table sizes
docker exec toxipred-postgres psql -U toxipred -d toxipred -c "SELECT relname, pg_size_pretty(pg_total_relation_size(relid)) FROM pg_catalog.pg_statio_user_tables ORDER BY pg_total_relation_size(relid) DESC;"

# Connection count
docker exec toxipred-postgres psql -U toxipred -c "SELECT count(*) FROM pg_stat_activity;"

# Vacuum database
docker exec toxipred-postgres psql -U toxipred -d toxipred -c "VACUUM ANALYZE;"
```

---

**For basic setup and usage, return to:** `SETUP.md`
