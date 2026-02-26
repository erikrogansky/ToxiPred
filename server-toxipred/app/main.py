from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from app.db import Base, engine
from app.api.routes_meta import router as meta_router
from app.api.routes_jobs import router as jobs_router
from app.api.routes_sharing import router as sharing_router

Base.metadata.create_all(bind=engine)

app = FastAPI(title="TOXIPRED API", version="0.1.0")
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_credentials=True,
                   allow_methods=["*"], allow_headers=["*"])

app.include_router(meta_router)
app.include_router(jobs_router)
app.include_router(sharing_router)


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
