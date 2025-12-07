from datetime import datetime, timedelta
from sqlalchemy import Column, String, DateTime, JSON, ARRAY
from app.db import Base

class JobResult(Base):
  __tablename__ = "job_results"

  id = Column(String, primary_key=True)
  model = Column(String, nullable=False)
  name = Column(String, nullable=True)
  trivial_name = Column(String, nullable=True)
  formula = Column(String, nullable=True)
  canonical_smiles = Column(String, nullable=True)
  other_names = Column(ARRAY(String), nullable=True)

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
      trivial_name=payload.get("trivial_name"),
      formula=payload.get("formula"),
      canonical_smiles=payload.get("canonical_smiles"),
      other_names=payload.get("other_names"),
      payload=payload,
      created_at=now,
      expires_at=now + timedelta(days=ttl_days),
    )
