# app/workers/tasks.py
import os
import redis
from celery import shared_task
from rdkit import Chem

from app.models.models import MODELS
from app.services.model_loader import get_model, resolve_feature_names
from app.domain.descriptors import descriptor_map
from app.services.predict import predict_vectorized

# Use the docker service name (not "redis")
PROGRESS_DB = int(os.getenv("PROGRESS_DB", "2"))
r = redis.Redis(host="toxipred-redis", port=6379, db=PROGRESS_DB)

def set_progress(job_id, pct, msg=""):
    r.hset(f"progress:{job_id}", mapping={"pct": pct, "msg": msg})

@shared_task(bind=True, name="toxipred.predict")
def predict_task(self, model_name: str, smiles: str):
    job_id = self.request.id
    set_progress(job_id, 1, "queued")

    m   = Chem.MolFromSmiles(smiles)
    can = Chem.MolToSmiles(m, canonical=True)
    feats = descriptor_map(m)

    spec  = MODELS[model_name]
    model = get_model(model_name, MODELS)
    names = resolve_feature_names(spec, model)

    import numpy as np
    X = np.array([feats[f] for f in names], dtype=float).reshape(1, -1)

    set_progress(job_id, 30, "featurized")
    y, conf = predict_vectorized(model, X)
    set_progress(job_id, 100, "done")

    return {
        "canonical_smiles": can,
        "prediction": y,
        "confidence": conf,
        "features_used": names,
        "model": model_name,
    }
