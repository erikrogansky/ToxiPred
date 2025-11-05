from fastapi import APIRouter
from pydantic import BaseModel
import redis, os

# ✅ use your Celery app instance
from app.workers.celery_app import celery_app

router = APIRouter(prefix="/api", tags=["jobs"])

PROGRESS_DB = int(os.getenv("PROGRESS_DB", "2"))
r = redis.Redis(host="toxipred-redis", port=6379, db=PROGRESS_DB)

class EnqueueRequest(BaseModel):
    query: str

@router.post("/jobs/predict/{model_name}")
def enqueue(model_name: str, req: EnqueueRequest):
    from rdkit import Chem
    m = Chem.MolFromSmiles(req.query)
    if not m:
        return {"error": "Invalid SMILES."}
    smiles = Chem.MolToSmiles(m, canonical=True)

    # send task via the same Celery app
    task = celery_app.send_task("toxipred.predict", args=[model_name, smiles])
    return {"job_id": task.id}

@router.get("/jobs/status/{job_id}")
def status(job_id: str):
    # ✅ bind AsyncResult to your Celery app
    ar = celery_app.AsyncResult(job_id)

    # progress hash from Redis (optional nicety)
    try:
        prog = {k.decode(): v.decode() for k, v in r.hgetall(f"progress:{job_id}").items()}
    except Exception:
        prog = {}

    payload = {
        "state": ar.state,  # PENDING / STARTED / SUCCESS / FAILURE / RETRY
        "progress_pct": (float(prog["pct"]) if "pct" in prog else None),
        "msg": prog.get("msg"),
    }

    if ar.successful():
        payload["result"] = ar.result
    elif ar.failed():
        payload["error"] = str(ar.result)

    return payload