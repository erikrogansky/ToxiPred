from typing import Any, List, Optional, Tuple
from fastapi import HTTPException
import json
from pathlib import Path
from app.models.models import ModelSpec

_cache: dict[str, Any] = {}

def load_pickle(spec: ModelSpec):
    if not spec.path.exists():
        raise HTTPException(500, f"Model file not found: {spec.path}")
    try:
        import joblib; return joblib.load(spec.path)
    except Exception: ...
    try:
        import pickle; return pickle.load(open(spec.path, "rb"))
    except Exception: ...
    try:
        import dill;   return dill.load(open(spec.path, "rb"))
    except Exception as e:
        raise HTTPException(500, f"Failed to load pickle '{spec.path.name}': {e}")

def get_model(name: str, specs: dict[str, ModelSpec]):
    spec = specs.get(name)
    if not spec: raise HTTPException(404, f"Unknown model '{name}'")
    if name not in _cache: _cache[name] = load_pickle(spec)
    return _cache[name]

def sidecar_path(spec: ModelSpec) -> Path:
    return spec.path.with_suffix(spec.path.suffix + ".features.json")

def load_sidecar_features(spec: ModelSpec) -> Optional[List[str]]:
    p = sidecar_path(spec)
    if p.exists():
        try:
            data = json.load(open(p, "r"))
            if isinstance(data, list) and all(isinstance(x, str) for x in data):
                return data
        except Exception: ...
    return None

def infer_feature_names(model: Any) -> Tuple[Optional[List[str]], Optional[int]]:
    names = None; n_in = getattr(model, "n_features_in_", None)
    if hasattr(model, "feature_names_in_"):
        try: names = [str(x) for x in model.feature_names_in_]
        except Exception: ...
    return names, (int(n_in) if n_in is not None else None)

def resolve_feature_names(spec: ModelSpec, model: Any) -> List[str]:
    if spec.features: return spec.features
    side = load_sidecar_features(spec)
    if side: return side
    names, n_in = infer_feature_names(model)
    if names: return names
    if n_in:
        raise HTTPException(500, f"Model reports n_features_in_={n_in} but no names. "
                                 f"Provide {sidecar_path(spec)} or set spec.features.")
    raise HTTPException(500, "Could not resolve feature names for the model.")
