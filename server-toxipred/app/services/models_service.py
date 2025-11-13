from pathlib import Path
import time
from typing import Dict, Any, Tuple
from fastapi import Request
from fastapi.responses import JSONResponse, Response

from app.models.models import MODELS
from app.utils.etag import make_etag

def models_payload() -> Dict[str, Any]:
  return {
    "available_models": list(MODELS.keys()),
    "details": {
      k: {
        "file": v.path.name,
        "path": str(v.path),
        "kind": v.kind,
        "features_in_spec": (len(v.features) if v.features else None),
        "note": v.note,
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
