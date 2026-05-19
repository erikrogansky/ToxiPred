from pathlib import Path
import time
from typing import Dict, Any, Tuple
from fastapi import Request
from fastapi.responses import JSONResponse, Response

from app.models.models import MODELS
from app.services.model_loader import load_sidecar_features, sidecar_paths
from app.utils.etag import make_etag

def feature_count(v) -> int | None:
  if v.features:
    return len(v.features)
  sidecar = load_sidecar_features(v)
  return len(sidecar) if sidecar else None

def models_payload() -> Dict[str, Any]:
  return {
    "available_models": list(MODELS.keys()),
    "details": {
      k: {
        "file": v.path.name,
        "path": str(v.path),
        "kind": v.kind,
        "features_in_spec": feature_count(v),
        "note": v.note,
        "test_type": v.test_type,
        "prediction_target": v.prediction_target,
        "positive_label": v.positive_label,
        "negative_label": v.negative_label,
        "classification_threshold": v.classification_threshold,
        "dataset": v.dataset,
      }
      for k, v in MODELS.items()
    },
  }

def last_models_mtime() -> float:
  mtimes = []
  for v in MODELS.values():
    try:
      mtimes.append(Path(v.path).stat().st_mtime)
    except FileNotFoundError:
      pass
    for sidecar in sidecar_paths(v):
      try:
        mtimes.append(sidecar.stat().st_mtime)
      except FileNotFoundError:
        pass
  return max(mtimes) if mtimes else time.time()

def list_models_response(request: Request) -> Response:
  payload = models_payload()
  last_mod = last_models_mtime()
  etag = make_etag(payload, last_mod)

  if request.headers.get("if-none-match") == etag:
    return Response(
      status_code=304,
      headers={
        "ETag": etag,
        "Cache-Control": "no-cache, must-revalidate",
        "Last-Modified": time.strftime("%a, %d %b %Y %H:%M:%S GMT", time.gmtime(last_mod)),
      },
    )

  return JSONResponse(
    content=payload,
    headers={
      "ETag": etag,
      "Cache-Control": "no-cache, must-revalidate",
      "Last-Modified": time.strftime("%a, %d %b %Y %H:%M:%S GMT", time.gmtime(last_mod)),
    },
  )
