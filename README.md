# 🧪 TOXIPRED – Toxicity Prediction Platform

Welcome to **TOXIPRED**, a public research project for predicting chemical toxicity using machine learning.
This repository contains both the frontend (**app-toxipred**) and backend (**server-toxipred**) components,
connected via Docker for easy development and deployment.

---

## 🚀 Quick Start

### 1. Prerequisites
- [Docker Desktop](https://www.docker.com/products/docker-desktop/) installed and running
- Git clone of this repository

```bash
git clone https://github.com/erikrogansky/toxipred.git
cd toxipred
```

---

### 2. Start the Development Environment

Build and run everything (frontend + backend):

```bash
docker compose up --build
```

After the first successful build, you can simply start it later with:

```bash
docker compose up
```

Then open:

- **Frontend (app-toxipred):** http://localhost:9000  
- **Backend (server-toxipred):** http://localhost:8000/health

---

### 3. Project Structure

```
toxipred/
├─ app-toxipred/           # Quasar (Vue 3 + TypeScript) frontend
│  ├─ src/
│  ├─ quasar.config.ts
│  ├─ Dockerfile
│  └─ .env.development
│
├─ server-toxipred/        # FastAPI backend
│  ├─ app/main.py
│  ├─ pyproject.toml
│  ├─ requirements.txt (auto-generated)
│  └─ Dockerfile
│
├─ docker-compose.yml      # Orchestrates frontend + backend
└─ README.md               # You are here ✨
```

---

### 4. Useful Commands

```bash
# Rebuild containers if dependencies changed
docker compose build

# Run in background (detached mode)
docker compose up -d

# Stop everything
docker compose down

# View logs (live)
docker compose logs -f

# Restart only one service (example: backend)
docker compose restart server-toxipred
```

---

### 5. Notes for Contributors

- Both containers use **live volume mounts** — any code changes reload automatically:
  - Frontend → Hot Module Reload (HMR) via Quasar
  - Backend → Auto-reload via Uvicorn

- Environment variables for the frontend live in `app-toxipred/.env.development`:
  ```env
  VITE_API_BASE=http://localhost:8000
  ```

- Backend endpoints (FastAPI) are defined in `server-toxipred/app/main.py`.

---

### 6. Production Build (optional)

To build the static frontend bundle for deployment:

```bash
docker compose exec app-toxipred npm run build
```

This will generate optimized static files inside `/app/dist/spa` within the container.

---

### 🧠 Summary

| Action | Command |
|--------|----------|
| Build + Run (first time) | `docker compose up --build` |
| Start again (afterwards) | `docker compose up` |
| Stop all containers | `docker compose down` |
| Open frontend | [http://localhost:9000](http://localhost:9000) |
| Open backend | [http://localhost:8000/health](http://localhost:8000/health) |

---

Developed with ❤️ by the **TOXIPRED** team  