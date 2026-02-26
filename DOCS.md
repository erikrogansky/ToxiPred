# 📚 Documentation Index

## Quick Navigation

### 🎯 Getting Started
- **[README.md](README.md)** - Project overview & development setup
- **[SETUP.md](SETUP.md)** - ⭐ Production deployment setup (start here!)

### 🔧 Advanced
- **[ADVANCED.md](ADVANCED.md)** - Detailed troubleshooting, architecture, optimization

---

## What Each File Does

### SETUP.md
**Use this to:** Deploy ToxiPred to production

Contains:
- Step-by-step server setup
- GitHub Actions configuration
- Daily usage commands
- Basic troubleshooting
- Quick command reference

**Time needed:** 10 minutes for setup, then automatic deployments

---

### ADVANCED.md
**Use this when:** You need to debug issues or customize the deployment

Contains:
- Architecture diagrams
- Database management
- Performance optimization
- Security hardening
- Advanced troubleshooting
- CI/CD customization
- Monitoring setup

**Reference this for:** Problems not solved in SETUP.md

---

## Quick Commands

```bash
# Setup production (first time)
# Follow: SETUP.md

# Daily deployment
git push origin main           # Auto-deploys
git commit -m "msg [skip-deploy]"  # Don't deploy

# Check status
ssh user@server
make health

# Rollback if needed
make rollback

# Advanced debugging
# See: ADVANCED.md
```

---

## File Structure

```
ToxiPred/
├── README.md          ← Project overview, development setup
├── SETUP.md           ← Production deployment guide ⭐
├── ADVANCED.md        ← Advanced troubleshooting & customization
├── DOCS.md            ← This file (documentation index)
│
├── docker-compose.yml              ← Development
├── docker-compose.prod.yml         ← Production
├── .env.example                    ← Environment template
│
├── scripts/
│   ├── deploy.sh                   ← Main deployment script
│   ├── migrate.sh                  ← Database migrations
│   ├── health-check.sh             ← Health monitoring
│   └── server-setup.sh             ← Server preparation
│
├── nginx/                          ← Load balancer config
└── .github/workflows/deploy.yml    ← Auto-deployment
```

---

**Need help?**
1. Check **SETUP.md** for common tasks
2. Check **ADVANCED.md** for detailed troubleshooting
3. Run `make health` on your server
