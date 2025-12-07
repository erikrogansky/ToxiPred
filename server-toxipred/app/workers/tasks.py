import os
import redis
import numpy as np
from celery import shared_task
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from datetime import datetime

from app.db import SessionLocal
from app.models.job_result import JobResult
from app.models.models import MODELS
from app.services.model_loader import get_model, resolve_feature_names
from app.domain.descriptors import descriptor_map
from app.services.predict import predict_vectorized
from app.services.explain import explain_atom_importance, explain_feature_importance

PROGRESS_DB = int(os.getenv("PROGRESS_DB", "2"))
r = redis.Redis(host="toxipred-redis", port=6379, db=PROGRESS_DB)

def set_progress(job_id, pct, msg=""):
    r.hset(f"progress:{job_id}", mapping={"pct": pct, "msg": msg})


def make_descriptor_fn(feature_names):
    def _fn(mol):
        feats = descriptor_map(mol)
        vec = np.array([feats.get(f, np.nan) for f in feature_names], dtype=float)
        if np.isnan(vec).any():
            vec = np.nan_to_num(vec, nan=0.0)
        return vec
    return _fn


@shared_task(bind=True, name="toxipred.predict")
def predict_task(
    self,
    model_name: str,
    smiles: str,
    display_name: str | None = None,
    formula: str | None = None,
    input_query: str | None = None,
    input_type: str | None = None,
    trivial_name: str | None = None,
    other_names: list[str] | None = None,
):
    job_id = self.request.id
    set_progress(job_id, 1, "queued")
    self.update_state(state="STARTED")

    m = Chem.MolFromSmiles(smiles)
    if m is None:
        set_progress(job_id, 100, "invalid_smiles")
        return {"error": "Invalid SMILES."}
    can = Chem.MolToSmiles(m, canonical=True)

    if formula is None:
        formula = rdMolDescriptors.CalcMolFormula(m)

    if not display_name:
        display_name = input_query or can

    spec  = MODELS[model_name]
    model = get_model(model_name, MODELS)
    names = resolve_feature_names(spec, model)

    feats_dict = descriptor_map(m)
    X = np.array([feats_dict.get(f, np.nan) for f in names], dtype=float).reshape(1, -1)
    if np.isnan(X).any():
        X = np.nan_to_num(X, nan=0.0)

    set_progress(job_id, 30, "featurized")
    self.update_state(state="PROGRESS", meta={"pct": 30, "msg": "featurized"})

    y, conf = predict_vectorized(model, X)

    set_progress(job_id, 55, "explaining")
    self.update_state(state="PROGRESS", meta={"pct": 55, "msg": "explaining"})

    feature_scores = None
    try:
        feature_scores_np, base_pred = explain_feature_importance(
            model,
            X,
            baseline=None,
            normalize=False,
        )
        feature_scores = feature_scores_np.tolist()
    except Exception as e:
        print("Error computing feature_scores:", e)
        feature_scores = None

    try:
        descriptor_fn = make_descriptor_fn(names)
        atom_scores_np, base_pred = explain_atom_importance(model, m, descriptor_fn)
        atom_scores = atom_scores_np.tolist()
    except Exception:
        atom_scores = None

    set_progress(job_id, 95, "assembling_result")

    # Extract feature values for the summary table
    feature_values = [feats_dict.get(f, None) for f in names]

    payload = {
        "input_query": input_query,
        "input_type": input_type,
        "name": display_name,
        "trivial_name": trivial_name,
        "formula": formula,
        "canonical_smiles": can,
        "other_names": other_names,

        "prediction": y,
        "confidence": conf,
        "features_used": names,
        "feature_values": feature_values,
        "model": model_name,

        "feature_scores": feature_scores,
        "atom_scores": atom_scores,
    }

    try:
        db = SessionLocal()
        db_obj = JobResult.from_payload(job_id, payload, ttl_days=14)
        existing = db.get(JobResult, job_id)
        if existing:
            db.delete(existing)
            db.flush()
        db.add(db_obj)
        db.commit()
    except Exception as e:
        print("Error saving JobResult:", e)
    finally:
        try:
            db.close()
        except:
            pass

    set_progress(job_id, 100, "done")
    self.update_state(state="PROGRESS", meta={"pct": 100, "msg": "done"})
    return payload


@shared_task
def prune_old_results():
    db = SessionLocal()
    try:
        now = datetime.utcnow()
        db.query(JobResult).filter(JobResult.expires_at < now).delete()
        db.commit()
    finally:
        db.close()
