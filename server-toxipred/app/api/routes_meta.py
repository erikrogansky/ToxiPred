from fastapi import APIRouter, HTTPException, Query
from rdkit import Chem
from app.models.models import MODELS
from app.services.model_loader import get_model, resolve_feature_names
from app.domain.descriptors import descriptor_map

router = APIRouter(prefix="/api", tags=["meta"])

@router.get("/models")
def list_models():
    return {
        "available_models": list(MODELS.keys()),
        "details": {k: {"file": v.path.name, "path": str(v.path),
                        "kind": v.kind, "features_in_spec": (len(v.features) if v.features else None),
                        "note": v.note} for k, v in MODELS.items()}
    }

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
