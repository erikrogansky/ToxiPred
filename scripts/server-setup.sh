#!/bin/bash

###############################################################################
# Server Setup Script for GitHub Actions Deployment
# 
# Run this on your server to prepare it for GitHub Actions deployments
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

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

echo ""
log_info "====== ToxiPred Server Setup for GitHub Actions ======"
echo ""

# Check if Docker is installed
log_info "Checking Docker installation..."
if ! command -v docker &> /dev/null; then
    log_error "Docker is not installed!"
    log_info "Install Docker first: https://docs.docker.com/engine/install/"
    exit 1
fi
log_success "Docker is installed"

# Check if Docker Compose is installed
log_info "Checking Docker Compose installation..."
if ! command -v docker compose &> /dev/null; then
    log_error "Docker Compose is not installed!"
    log_info "Install Docker Compose: https://docs.docker.com/compose/install/"
    exit 1
fi
log_success "Docker Compose is installed"

# Check if git is installed
log_info "Checking Git installation..."
if ! command -v git &> /dev/null; then
    log_error "Git is not installed!"
    log_info "Install git: sudo apt-get install git"
    exit 1
fi
log_success "Git is installed"

# Print current directory info
echo ""
log_info "Current directory: $(pwd)"
log_info "Current user: $(whoami)"
echo ""

# Ask for repository URL if not already cloned
if [ ! -d ".git" ]; then
    log_warning "This doesn't appear to be a git repository"
    echo ""
    echo "Please clone the repository first:"
    echo "  git clone https://github.com/erikrogansky/ToxiPred.git"
    echo "  cd ToxiPred"
    echo "  ./scripts/server-setup.sh"
    echo ""
    exit 1
fi

log_success "Git repository detected"

# Make scripts executable
log_info "Making scripts executable..."
chmod +x scripts/*.sh
log_success "Scripts are now executable"

# Create .env if it doesn't exist
if [ ! -f ".env" ]; then
    log_info "Creating .env from template..."
    cp .env.prod.example .env
    log_warning "Please edit .env and set secure passwords!"
    log_info "Edit with: nano .env"
else
    log_success ".env file already exists"
fi

# Create backups directory
log_info "Creating backups directory..."
mkdir -p backups
log_success "Backups directory ready"

# Check SSH setup
echo ""
log_info "Checking SSH setup..."
if [ -f "$HOME/.ssh/authorized_keys" ]; then
    log_success "SSH authorized_keys file exists"
    log_info "Add your GitHub Actions deploy key to: $HOME/.ssh/authorized_keys"
else
    log_info "Creating SSH directory..."
    mkdir -p "$HOME/.ssh"
    chmod 700 "$HOME/.ssh"
    touch "$HOME/.ssh/authorized_keys"
    chmod 600 "$HOME/.ssh/authorized_keys"
    log_success "SSH directory created"
    log_warning "Add your GitHub Actions deploy key to: $HOME/.ssh/authorized_keys"
fi

# Check firewall
echo ""
log_info "Checking firewall recommendations..."
if command -v ufw &> /dev/null; then
    log_info "UFW firewall detected"
    echo "Recommended firewall rules:"
    echo "  sudo ufw allow 22    # SSH"
    echo "  sudo ufw allow 80    # HTTP"
    echo "  sudo ufw allow 443   # HTTPS"
    echo "  sudo ufw enable"
else
    log_info "UFW not detected, ensure your firewall allows:"
    echo "  - Port 22 (SSH)"
    echo "  - Port 80 (HTTP)"
    echo "  - Port 443 (HTTPS)"
fi

# Print GitHub Actions secrets info
echo ""
log_info "====== GitHub Secrets Configuration ======"
echo ""
echo "Add these secrets to your GitHub repository:"
echo ""
echo "1. SSH_PRIVATE_KEY"
echo "   Generate on your LOCAL machine:"
echo "   ssh-keygen -t ed25519 -C 'github-actions' -f ~/.ssh/toxipred_deploy"
echo "   cat ~/.ssh/toxipred_deploy  # Copy this to GitHub Secret"
echo ""
echo "2. SERVER_HOST"
echo "   Value: $(hostname -f 2>/dev/null || hostname)"
echo "   Or your server's public IP/domain"
echo ""
echo "3. SERVER_USER"
echo "   Value: $(whoami)"
echo ""
echo "4. DEPLOY_PATH"
echo "   Value: $(pwd)"
echo ""

# Print next steps
log_info "====== Next Steps ======"
echo ""
echo "1. Edit .env file and set secure passwords:"
echo "   nano .env"
echo ""
echo "2. Add GitHub Actions deploy key (public key) to authorized_keys:"
echo "   echo 'your-public-key' >> ~/.ssh/authorized_keys"
echo ""
echo "3. Test the deployment locally:"
echo "   make prod-up"
echo "   make health"
echo ""
echo "4. Configure GitHub Secrets (see above)"
echo ""
echo "5. Push to main branch to trigger auto-deployment!"
echo ""

log_success "====== Server Setup Complete! ======"
echo ""
log_info "Your deploy path: $(pwd)"
log_info "Your SSH user: $(whoami)"
echo ""
