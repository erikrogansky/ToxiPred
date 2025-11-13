from fastapi import APIRouter, HTTPException
from datetime import datetime
from pydantic import BaseModel
from app.services.chem_resolver import resolve_chemical
import redis, os
import json, asyncio
from fastapi.responses import StreamingResponse
from typing import List

from app.db import SessionLocal
from app.models.job_result import JobResult

class JobsValidateRequest(BaseModel):
    job_ids: List[str]

class JobsValidateResponse(BaseModel):
    invalid_ids: List[str]

TERMINAL = {"SUCCESS", "FAILURE", "REVOKED"}

from app.workers.celery_app import celery_app

router = APIRouter(tags=["jobs"])

PROGRESS_DB = int(os.getenv("PROGRESS_DB", "2"))
r = redis.Redis(host="toxipred-redis", port=6379, db=PROGRESS_DB)

class EnqueueRequest(BaseModel):
    query: str

@router.post("/jobs/predict/{model_name}")
def enqueue(model_name: str, req: EnqueueRequest):
    chem = resolve_chemical(req.query)

    if chem is None:
        return {"error": "Unable to interpret input as SMILES, CAS number, name, or InChI."}

    task = celery_app.send_task(
        "toxipred.predict",
        args=[model_name, chem.smiles],
        kwargs={
            "display_name": chem.name,
            "formula": chem.formula,
            "input_query": chem.query,
            "input_type": chem.source,
        },
    )

    return {"job_id": task.id}

@router.get("/jobs/status/{job_id}")
def status(job_id: str):
    ar = celery_app.AsyncResult(job_id)

    try:
        prog = {k.decode(): v.decode() for k, v in r.hgetall(f"progress:{job_id}").items()}
    except Exception:
        prog = {}

    payload = {
        "state": ar.state,
        "progress_pct": (float(prog["pct"]) if "pct" in prog else None),
        "msg": prog.get("msg"),
    }

    if ar.successful():
        res = ar.result or {}
        payload["result"] = {
            "name": res.get("name"),
            "formula": res.get("formula"),
            "model": res.get("model"),
            "prediction": res.get("prediction"),
        }
    elif ar.failed():
        payload["error"] = str(ar.result)

    return payload

@router.get("/jobs/stream/{job_id}")
async def stream(job_id: str):
    async def event_gen():
        last_pct = None
        last_state = None

        while True:
            ar = celery_app.AsyncResult(job_id)

            try:
                prog = {k.decode(): v.decode() for k, v in r.hgetall(f"progress:{job_id}").items()}
            except Exception:
                prog = {}

            payload = {
                "task_id": job_id,
                "state": ar.state,
                "progress_pct": (float(prog["pct"]) if "pct" in prog else None),
                "msg": prog.get("msg"),
            }

            if ar.successful():
                res = ar.result or {}
                payload["result"] = {
                    "name": res.get("name"),
                    "formula": res.get("formula"),
                    "model": res.get("model"),
                    "prediction": res.get("prediction"),
                }
            elif ar.failed():
                payload["error"] = str(ar.result)

            pct = payload["progress_pct"]
            if payload["state"] != last_state or pct != last_pct or payload.get("error") or payload.get("result"):
                yield f"data: {json.dumps(payload)}\n\n"
                last_state = payload["state"]
                last_pct = pct

            if payload["state"] in TERMINAL:
                break

            await asyncio.sleep(0.7)

    return StreamingResponse(
        event_gen(),
        media_type="text/event-stream",
        headers={
            "Cache-Control": "no-cache",
            "Connection": "keep-alive",
        },
    )


@router.get("/jobs/result/{job_id}")
def job_result(job_id: str):
    db = SessionLocal()
    try:
        obj = db.get(JobResult, job_id)
        if obj is None:
            ar = celery_app.AsyncResult(job_id)
            if not ar.ready():
                return {"error": "Job not finished", "state": ar.state}
            if ar.failed():
                return {"error": str(ar.result), "state": ar.state}
            return ar.result

        now = datetime.utcnow()
        if obj.expires_at < now:
            raise HTTPException(status_code=404, detail="Result expired")

        return obj.payload
    finally:
        db.close()

@router.post("/jobs/validate", response_model=JobsValidateResponse)
def validate_jobs(body: JobsValidateRequest):
    job_ids = list(set(body.job_ids))
    if not job_ids:
        return JobsValidateResponse(invalid_ids=[])

    db = SessionLocal()
    try:
        now = datetime.utcnow()
        rows = (
            db.query(JobResult)
            .filter(JobResult.id.in_(job_ids))
            .all()
        )

        found_ids = {row.id for row in rows}
        invalid_ids = set(job_ids) - found_ids

        expired_ids = {row.id for row in rows if row.expires_at < now}
        invalid_ids |= expired_ids

        if expired_ids:
            (
                db.query(JobResult)
                .filter(JobResult.id.in_(list(expired_ids)))
                .delete(synchronize_session=False)
            )
            db.commit()

        return JobsValidateResponse(invalid_ids=sorted(invalid_ids))
    finally:
        db.close()


@router.delete("/jobs/delete/{job_id}")
def delete_job(job_id: str):
    db = SessionLocal()
    try:
        deleted_count = (
            db.query(JobResult)
            .filter(JobResult.id == job_id)
            .delete(synchronize_session=False)
        )
        db.commit()

        if deleted_count == 0:
            raise HTTPException(status_code=404, detail="Job not found")

        return {"message": "Job deleted successfully", "job_id": job_id}
    finally:
        db.close()
