from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from app.db import Base, engine
from app.api.routes_meta import router as meta_router
from app.api.routes_jobs import router as jobs_router

Base.metadata.create_all(bind=engine)

app = FastAPI(title="TOXIPRED API", version="0.1.0")
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_credentials=True,
                   allow_methods=["*"], allow_headers=["*"])

app.include_router(meta_router)
app.include_router(jobs_router)
