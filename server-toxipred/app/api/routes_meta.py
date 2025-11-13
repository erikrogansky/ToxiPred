from fastapi import APIRouter, HTTPException, Query, Request
from rdkit import Chem

from app.models.models import MODELS
from app.services.model_loader import get_model, resolve_feature_names
from app.services.models_service import list_models_response
from app.domain.descriptors import descriptor_map
from app.schemas.meta import ModelsResponse, IntrospectResponse

router = APIRouter(tags=["meta"])

@router.get("/models", response_model=ModelsResponse)
def list_models(request: Request):
  return list_models_response(request)

@router.get("/descriptors")
def describe_smiles(smiles: str = Query(...)):
    m = Chem.MolFromSmiles(smiles)
    if not m: raise HTTPException(400, "Invalid SMILES.")
    can = Chem.MolToSmiles(m, canonical=True)
    feats = descriptor_map(m)
    return {"canonical_smiles": can, "computed_features": sorted(feats.keys()), "values": feats}

@router.get("/models/{model_name}/introspect")
def introspect_model(model_name: str):
    spec = MODELS.get(model_name)
    if not spec: raise HTTPException(404, f"Unknown model '{model_name}'")
    model = get_model(model_name, MODELS)
    names = resolve_feature_names(spec, model)
    return {"model": model_name, "file": str(spec.path), "used_features": names, "note": spec.note or None}
