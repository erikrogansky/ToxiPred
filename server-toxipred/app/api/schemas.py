from pydantic import BaseModel, Field
from typing import Literal, List, Optional, Union, Dict

class PredictRequest(BaseModel):
    query: str = Field(..., description="SMILES or InChI")

class PredictResponse(BaseModel):
    job_id: str

class JobStatus(BaseModel):
    state: str
    progress_pct: Optional[float] = None
    msg: Optional[str] = None
    result: Optional[Dict] = None
