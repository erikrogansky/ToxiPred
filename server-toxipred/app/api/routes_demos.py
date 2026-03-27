from fastapi import APIRouter
from app.db import SessionLocal
from app.models.job_result import JobResult

router = APIRouter(tags=["demos"])


@router.get("/demos")
def list_demos():
    """Return all demo predictions."""
    db = SessionLocal()
    try:
        demos = (
            db.query(JobResult)
            .filter(JobResult.is_demo == True)
            .order_by(JobResult.created_at.asc())
            .all()
        )
        return [
            {
                "job_id": d.id,
                "model": d.model,
                "name": d.name,
                "trivial_name": d.trivial_name,
                "formula": d.formula,
                "canonical_smiles": d.canonical_smiles,
                "demo_title": d.demo_title,
                "demo_description": d.demo_description,
                "prediction": d.payload.get("prediction") if d.payload else None,
            }
            for d in demos
        ]
    finally:
        db.close()
