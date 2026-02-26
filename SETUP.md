# 🚀 ToxiPred Production Deployment - Setup Guide

**Everything you need to deploy ToxiPred to production with zero-downtime and auto-deploy from GitHub.**

---

## 📋 What You'll Get

- ✅ **Zero-downtime deployments** (blue-green strategy for backend, optimized build for frontend)
- ✅ **Auto-deploy on git push** to main branch
- ✅ **Automatic rollback** if health checks fail
- ✅ **Database backups** before each deployment
- ✅ **Data persistence** (no data loss ever)
- ✅ **Production-optimized frontend** (built Quasar bundle, not dev server)
- ✅ **Simple commands:** `git push` to deploy, `make rollback` to undo

**Time to setup:** 10 minutes  
**Deployment time:** Automatic (~3 minutes after push, includes Quasar build)

---

## 🎯 Step-by-Step Setup

### Part 1: Prepare Your Local Machine (2 minutes)

On your **Mac/PC**, generate an SSH key for GitHub Actions:

```bash
# Generate SSH key (press Enter for no passphrase)
ssh-keygen -t ed25519 -C "github-actions" -f ~/.ssh/toxipred_deploy

# Display private key (you'll need this for GitHub)
cat ~/.ssh/toxipred_deploy
# Copy the ENTIRE output (including -----BEGIN/END-----)

# Display public key (you'll need this for server)
cat ~/.ssh/toxipred_deploy.pub
# Copy this line
```

**Save these two keys somewhere - you'll need them in the next steps!**

---

### Part 2: Setup Your Server (5 minutes)

SSH into your server and run these commands:

```bash
# 1. Clone the repository
git clone https://github.com/erikrogansky/ToxiPred.git
cd ToxiPred

# 2. Run the setup script
./scripts/server-setup.sh

# 3. Add GitHub Actions public key (from Part 1)
nano ~/.ssh/authorized_keys
# Paste the public key on a new line, save (Ctrl+X, Y, Enter)

# 4. Configure environment
cp .env.prod.example .env
nano .env
# Change POSTGRES_PASSWORD to something secure
# Save (Ctrl+X, Y, Enter)

# 5. Start the application
make prod-up

# 6. Wait 30 seconds for everything to start
sleep 30

# 7. Verify it's working
make health
# You should see green checkmarks ✓
```

**Note the following from your server:**
- Server hostname/IP: `hostname -f` (or your public IP)
- SSH username: `whoami`
- Deploy path: `pwd` (run this inside the ToxiPred directory)

---

### Part 3: Configure GitHub Secrets (2 minutes)

1. Go to your GitHub repository: `https://github.com/erikrogansky/ToxiPred`

2. Click: **Settings** → **Secrets and variables** → **Actions**

3. Click: **New repository secret**

4. Add these 4 secrets (one at a time):

| Secret Name | Value | Example |
|-------------|-------|---------|
| `SSH_PRIVATE_KEY` | Private key from Part 1 | `-----BEGIN OPENSSH PRIVATE KEY-----...` |
| `SERVER_HOST` | Your server hostname or IP | `toxipred.example.com` or `192.168.1.100` |
| `SERVER_USER` | Your SSH username | `ubuntu` or `your-username` |
| `DEPLOY_PATH` | Full path to ToxiPred directory | `/home/ubuntu/ToxiPred` |

---

### Part 4: Test Auto-Deploy (1 minute)

```bash
# On your local machine, make a test commit
cd /path/to/your/local/ToxiPred
echo "# Test auto-deploy" >> README.md
git add README.md
git commit -m "Test automatic deployment"
git push origin main

# Watch the deployment happen!
# Go to: https://github.com/erikrogansky/ToxiPred/actions
```

You should see:
- ✅ Green checkmark when deployment succeeds
- Your server automatically updated with the latest code

---

## 🎮 Daily Usage

### Deploy New Code

```bash
# Make your changes
git add .
git commit -m "Add awesome new feature"
git push origin main

# That's it! GitHub Actions deploys automatically
# Watch progress at: github.com/erikrogansky/ToxiPred/actions
```

### Skip Deployment (for docs/minor changes)

```bash
git commit -m "Update README [skip-deploy]"
git push origin main
# Pushes to GitHub but doesn't deploy to server
```

### Manually Trigger Deployment

1. Go to: `github.com/erikrogansky/ToxiPred/actions`
2. Click: **Deploy to Production**
3. Click: **Run workflow** → **Run workflow**

### Check Status on Server

```bash
ssh your-user@your-server.com
cd /path/to/ToxiPred

make health        # Comprehensive health check
make prod-status   # Current deployment status
make prod-logs     # View live logs
```

### Rollback if Needed

```bash
ssh your-user@your-server.com
cd /path/to/ToxiPred
make rollback
# Switches back to previous version in ~10 seconds
```

---

## 🔧 Common Tasks

### Update Environment Variables

```bash
ssh your-user@your-server.com
cd /path/to/ToxiPred
nano .env
# Make your changes
make prod-down
make prod-up
```

### View Application Logs

```bash
ssh your-user@your-server.com
cd /path/to/ToxiPred

# All logs
make prod-logs

# Specific service
docker compose -f docker-compose.prod.yml logs -f server-toxipred-blue

# Last 100 lines
docker compose -f docker-compose.prod.yml logs --tail=100
```

### Manual Backup

```bash
ssh your-user@your-server.com
cd /path/to/ToxiPred
make backup
# Backup saved in ./backups/
```

### Check Disk Space

```bash
ssh your-user@your-server.com
df -h                    # Overall disk usage
docker system df         # Docker disk usage
make clean              # Clean up old Docker images
```

### Scale Workers (if you need more processing power)

```bash
ssh your-user@your-server.com
cd /path/to/ToxiPred
docker compose -f docker-compose.prod.yml up -d --scale toxipred-worker=3
# Now you have 3 workers instead of 1
```

---

## 🆘 Quick Troubleshooting

### GitHub Actions Deployment Failing

**Check the logs:**
1. Go to: `github.com/erikrogansky/ToxiPred/actions`
2. Click on the failed workflow
3. Read the error message

**Common issues:**

- **"Permission denied (publickey)"**
  ```bash
  # Test SSH connection from local machine:
  ssh -i ~/.ssh/toxipred_deploy your-user@your-server.com
  
  # If fails, check:
  # 1. Public key is in server's ~/.ssh/authorized_keys
  # 2. Private key is correct in GitHub Secret SSH_PRIVATE_KEY
  ```

- **"Directory not found"**
  ```bash
  # Verify DEPLOY_PATH matches actual path on server
  ssh your-user@your-server.com
  cd /path/to/ToxiPred && pwd
  # Update DEPLOY_PATH secret if different
  ```

### Application Not Working on Server

```bash
# SSH to server and run health check
ssh your-user@your-server.com
cd /path/to/ToxiPred
make health

# Check which services are down
docker compose -f docker-compose.prod.yml ps

# View logs of failing service
docker compose -f docker-compose.prod.yml logs service-name

# Restart everything if needed
make prod-down
make prod-up
```

### Health Checks Failing

```bash
ssh your-user@your-server.com
cd /path/to/ToxiPred

# Run detailed health check
./scripts/health-check.sh

# Check specific issues:
docker compose -f docker-compose.prod.yml ps          # Container status
docker compose -f docker-compose.prod.yml logs -f     # Live logs
make prod-status                                      # Deployment status
```

### Out of Disk Space

```bash
ssh your-user@your-server.com
df -h                          # Check disk usage
docker system df               # Check Docker disk usage

# Clean up
make clean                     # Remove old Docker images
docker system prune -a         # More aggressive cleanup

# Remove old backups (keeps last 10 automatically)
ls -t backups/*.sql | tail -n +11 | xargs rm
```

### Need to Rollback

```bash
# If deployment just failed, GitHub Actions already rolled back
# But if you need to manually rollback:
ssh your-user@your-server.com
cd /path/to/ToxiPred
make rollback
```

### Database Issues

```bash
ssh your-user@your-server.com
cd /path/to/ToxiPred

# Check database is running
docker compose -f docker-compose.prod.yml ps toxipred-postgres

# Check database logs
docker compose -f docker-compose.prod.yml logs toxipred-postgres

# Connect to database
docker exec -it toxipred-postgres psql -U toxipred toxipred

# Restore from backup (if needed)
docker exec -i toxipred-postgres psql -U toxipred toxipred < backups/toxipred_backup_YYYYMMDD_HHMMSS.sql
```

---

## 🔒 Security Checklist

Before going live:

- [ ] Changed `POSTGRES_PASSWORD` in `.env` to a strong password
- [ ] `.env` file is NOT committed to git (it's in .gitignore)
- [ ] GitHub Secrets are properly configured
- [ ] Firewall allows only ports 22 (SSH), 80 (HTTP), 443 (HTTPS)
- [ ] SSH key authentication is working
- [ ] Using HTTPS in production (optional but recommended)

---

## 📊 Monitoring Your Application

### Application URLs

- **Frontend:** `http://your-server.com`
- **API Health:** `http://your-server.com/api/health`
- **Flower (Celery monitoring):** `http://your-server.com:5555`

### Regular Health Checks

```bash
# Automated health check
ssh your-user@your-server.com 'cd /path/to/ToxiPred && make health'

# Or setup a cron job to check daily
# Add to crontab: crontab -e
0 9 * * * cd /path/to/ToxiPred && ./scripts/health-check.sh > /tmp/health.log 2>&1
```

### View Recent Deployments

- **GitHub:** `github.com/erikrogansky/ToxiPred/actions`
- **Server:** `ssh your-user@your-server.com 'cd /path/to/ToxiPred && git log --oneline -10'`

---

## 🎓 Understanding How It Works

### Blue-Green Deployment

```
┌─────────┐
│  Nginx  │  ← Load balancer (always running)
└────┬────┘
     │
     ├─→ [Blue Server]   ← Currently active (serving users)
     │
     └─→ [Green Server]  ← Standby (deployment target)
```

**Deployment process:**
1. Build new code on Green (standby)
2. Run health checks on Green
3. Switch Nginx to point to Green
4. Green is now active, Blue becomes standby
5. Next deployment goes to Blue

**Result:** Users never see downtime!

### Frontend Deployment

The frontend (Quasar) is handled differently but also with minimal interruption:

1. **Build Phase:** Quasar builds production bundle (happens during Docker build, ~1-2 minutes)
2. **Deploy Phase:** New nginx container with built files replaces old one (~5 seconds)
3. **Nginx serves static files:** No backend needed, just HTML/CSS/JS
4. **Browser caching:** Users on site continue using cached version

**Frontend downtime:** ~5 seconds (only during container swap)
**Backend downtime:** 0 seconds (blue-green deployment)

### What Persists Across Deployments

- ✅ PostgreSQL database (all jobs, results, shared data)
- ✅ Redis queue (tasks being processed)
- ✅ Docker volumes (mounted at `/var/lib/docker/volumes/`)

**These NEVER get deleted during deployments or container restarts!**

### Automatic Rollback

If health checks fail during deployment:
1. Deployment script detects failure
2. Traffic stays on old version
3. New containers are stopped
4. Deployment marked as failed
5. You get notification in GitHub Actions

---

## 💡 Pro Tips

1. **Always test locally first:**
   ```bash
   make up          # Test in development mode
   make build       # Ensure it builds correctly
   ```

2. **Use descriptive commit messages:**
   ```bash
   git commit -m "Fix user authentication bug"
   # Better than: git commit -m "fix"
   ```

3. **Skip deployment for non-code changes:**
   ```bash
   git commit -m "Update documentation [skip-deploy]"
   ```

4. **Monitor deployments:**
   - Keep GitHub Actions tab open during deployment
   - Watch for green checkmark ✅

5. **Backup before major changes:**
   ```bash
   ssh your-user@your-server.com 'cd /path/to/ToxiPred && make backup'
   ```

---

## 📞 Quick Command Reference

```bash
# Deployment (happens automatically on push)
git push origin main                    # Auto-deploy
git commit -m "message [skip-deploy]"  # Don't deploy

# Server operations
make health              # Health check
make prod-status         # Deployment status
make prod-logs           # View logs
make rollback            # Undo deployment
make backup              # Manual backup
make clean               # Clean up Docker

# Monitoring
make health              # Quick health check
./scripts/health-check.sh  # Detailed health report
```

---

## ✅ You're All Set!

Your deployment pipeline is now:

1. **Write code** on your machine
2. **Commit and push** to GitHub
3. **Automatic deployment** to your server
4. **Zero downtime** for users
5. **Automatic rollback** if anything fails

That's it! Professional DevOps without the complexity. 🚀

---

**For advanced troubleshooting and detailed information, see:** `ADVANCED-DEPLOYMENT.md`
