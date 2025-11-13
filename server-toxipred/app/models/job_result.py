from datetime import datetime, timedelta
from sqlalchemy import Column, String, DateTime, JSON
from app.db import Base

class JobResult(Base):
  __tablename__ = "job_results"

  id = Column(String, primary_key=True)
  model = Column(String, nullable=False)
  name = Column(String, nullable=True)
  formula = Column(String, nullable=True)
  canonical_smiles = Column(String, nullable=True)

  payload = Column(JSON, nullable=False)

  created_at = Column(DateTime, nullable=False, default=datetime.utcnow, index=True)
  expires_at = Column(DateTime, nullable=False, index=True)

  @staticmethod
  def from_payload(job_id: str, payload: dict, ttl_days: int = 14) -> "JobResult":
    now = datetime.utcnow()
    return JobResult(
      id=job_id,
      model=payload.get("model"),
      name=payload.get("name"),
      formula=payload.get("formula"),
      canonical_smiles=payload.get("canonical_smiles"),
      payload=payload,
      created_at=now,
      expires_at=now + timedelta(days=ttl_days),
    )
