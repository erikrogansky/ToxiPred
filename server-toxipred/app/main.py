from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from sqlalchemy import inspect, text
from app.db import Base, engine
from app.api.routes_meta import router as meta_router
from app.api.routes_jobs import router as jobs_router
from app.api.routes_sharing import router as sharing_router
from app.api.routes_demos import router as demos_router

Base.metadata.create_all(bind=engine)

# Safe migration: add demo columns if they don't exist yet
_inspector = inspect(engine)
if _inspector.has_table("job_results"):
    _existing = {c["name"] for c in _inspector.get_columns("job_results")}
    with engine.begin() as conn:
        if "is_demo" not in _existing:
            conn.execute(text(
                "ALTER TABLE job_results ADD COLUMN is_demo BOOLEAN NOT NULL DEFAULT false"
            ))
        if "demo_title" not in _existing:
            conn.execute(text(
                "ALTER TABLE job_results ADD COLUMN demo_title VARCHAR"
            ))
        if "demo_description" not in _existing:
            conn.execute(text(
                "ALTER TABLE job_results ADD COLUMN demo_description VARCHAR"
            ))

app = FastAPI(title="TOXIPRED API", version="0.1.0")
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_credentials=True,
                   allow_methods=["*"], allow_headers=["*"])

app.include_router(meta_router)
app.include_router(jobs_router)
app.include_router(sharing_router)
app.include_router(demos_router)


@app.get("/health")
async def health_check():
    """Health check endpoint for load balancer and deployment monitoring."""
    return {"status": "healthy", "service": "toxipred-api"}


@app.get("/ready")
async def readiness_check():
    """Readiness check - verifies database connectivity."""
    try:
        from app.db import SessionLocal
        db = SessionLocal()
        db.execute("SELECT 1")
        db.close()
        return {"status": "ready", "database": "connected"}
    except Exception as e:
        from fastapi import HTTPException
        raise HTTPException(status_code=503, detail=f"Database not ready: {str(e)}")
