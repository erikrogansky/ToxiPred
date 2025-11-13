from pydantic import BaseModel
from typing import Optional, Dict, List

class ModelDetail(BaseModel):
  file: str
  path: str
  kind: str
  features_in_spec: Optional[int] = None
  note: Optional[str] = None

class ModelsResponse(BaseModel):
  available_models: List[str]
  details: Dict[str, ModelDetail]

class IntrospectResponse(BaseModel):
  model: str
  file: str
  used_features: List[str]
  note: Optional[str] = None
