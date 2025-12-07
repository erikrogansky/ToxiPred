from fastapi import APIRouter, HTTPException
from datetime import datetime
from pydantic import BaseModel
from typing import Optional

from app.db import SessionLocal
from app.models.job_result import JobResult
from app.models.shared_job import SharedJob

router = APIRouter(tags=["sharing"])


class ShareRequest(BaseModel):
    """Request to create a share link."""
    custom_password: Optional[str] = None  # Optional custom password


class ShareResponse(BaseModel):
    """Response containing share link details."""
    share_token: str
    password: str
    expires_at: str


class AccessSharedRequest(BaseModel):
    """Request to access a shared job."""
    password: str


class AccessSharedResponse(BaseModel):
    """Response containing shared job data."""
    job_data: dict
    can_import: bool = True


class ImportSharedRequest(BaseModel):
    """Request to import a shared job into workspace."""
    password: str
    existing_job_ids: Optional[list[str]] = None  # Job IDs already in user's workspace


class ImportSharedResponse(BaseModel):
    """Response after importing a shared job."""
    job_id: str
    message: str
    already_exists: bool = False  # True if job was already in workspace


@router.post("/jobs/share/{job_id}", response_model=ShareResponse)
def create_share(job_id: str, req: Optional[ShareRequest] = None):
    """
    Create a shareable link for a job result.
    Returns a token and password that can be used to access the job.
    Optionally accepts a custom password in the request body.
    """
    db = SessionLocal()
    try:
        # Verify the job exists
        job = db.get(JobResult, job_id)
        if job is None:
            raise HTTPException(status_code=404, detail="Job not found")
        
        # Check if job is expired
        now = datetime.utcnow()
        if job.expires_at < now:
            raise HTTPException(status_code=404, detail="Job has expired")
        
        # Get custom password if provided
        custom_password = req.custom_password if req else None
        
        # Create the share
        shared, password = SharedJob.create_share(job_id, custom_password=custom_password, ttl_days=30)
        db.add(shared)
        db.commit()
        db.refresh(shared)
        
        return ShareResponse(
            share_token=shared.share_token,
            password=password,
            expires_at=shared.expires_at.isoformat(),
        )
    finally:
        db.close()


@router.post("/jobs/shared/{token}", response_model=AccessSharedResponse)
def access_shared(token: str, req: AccessSharedRequest):
    """
    Access a shared job result using token and password.
    Returns the job data if password is correct.
    """
    db = SessionLocal()
    try:
        # Find the shared job
        shared = db.query(SharedJob).filter(SharedJob.share_token == token).first()
        
        if shared is None:
            raise HTTPException(status_code=404, detail="Share link not found or has been deleted")
        
        # Check if share is expired
        now = datetime.utcnow()
        if shared.expires_at < now:
            raise HTTPException(status_code=410, detail="Share link has expired")
        
        # Verify password
        if not shared.verify(req.password):
            raise HTTPException(status_code=401, detail="Invalid password")
        
        # Get the source job
        job = db.get(JobResult, shared.source_job_id)
        if job is None:
            raise HTTPException(
                status_code=410, 
                detail="The shared prediction has been deleted by all users"
            )
        
        # Check if source job is expired
        if job.expires_at < now:
            raise HTTPException(status_code=410, detail="The prediction has expired")
        
        return AccessSharedResponse(
            job_data=job.payload,
            can_import=True,
        )
    finally:
        db.close()


@router.post("/jobs/import/{token}", response_model=ImportSharedResponse)
def import_shared(token: str, req: ImportSharedRequest):
    """
    Import a shared job into the user's workspace.
    This increments the reference count on the source job.
    Returns already_exists=true if the job is already in the user's workspace.
    """
    db = SessionLocal()
    try:
        # Find the shared job
        shared = db.query(SharedJob).filter(SharedJob.share_token == token).first()
        
        if shared is None:
            raise HTTPException(status_code=404, detail="Share link not found or has been deleted")
        
        # Check if share is expired
        now = datetime.utcnow()
        if shared.expires_at < now:
            raise HTTPException(status_code=410, detail="Share link has expired")
        
        # Verify password
        if not shared.verify(req.password):
            raise HTTPException(status_code=401, detail="Invalid password")
        
        # Get the source job
        job = db.get(JobResult, shared.source_job_id)
        if job is None:
            raise HTTPException(
                status_code=410, 
                detail="The shared prediction has been deleted by all users"
            )
        
        # Check if source job is expired
        if job.expires_at < now:
            raise HTTPException(status_code=410, detail="The prediction has expired")
        
        # Check if job already exists in user's workspace
        existing_ids = req.existing_job_ids or []
        if shared.source_job_id in existing_ids:
            return ImportSharedResponse(
                job_id=shared.source_job_id,
                message="This prediction is already in your workspace",
                already_exists=True,
            )
        
        # Increment reference count only for new imports
        shared.increment_refs()
        db.commit()
        
        return ImportSharedResponse(
            job_id=shared.source_job_id,
            message="Successfully imported prediction to workspace",
            already_exists=False,
        )
    finally:
        db.close()


@router.get("/jobs/shared/{token}/check")
def check_shared(token: str):
    """
    Check if a share link exists and is valid (without requiring password).
    Used to show appropriate UI state before password entry.
    """
    db = SessionLocal()
    try:
        shared = db.query(SharedJob).filter(SharedJob.share_token == token).first()
        
        if shared is None:
            return {"valid": False, "reason": "not_found"}
        
        now = datetime.utcnow()
        if shared.expires_at < now:
            return {"valid": False, "reason": "expired"}
        
        # Check if source job still exists
        job = db.get(JobResult, shared.source_job_id)
        if job is None:
            return {"valid": False, "reason": "deleted"}
        
        if job.expires_at < now:
            return {"valid": False, "reason": "job_expired"}
        
        return {"valid": True}
    finally:
        db.close()
