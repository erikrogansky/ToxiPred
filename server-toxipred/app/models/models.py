# app/models/models.py
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, Optional
import os

BASE_DIR = Path(__file__).resolve().parent
MODELS_DIR = Path(os.getenv("TOXIPRED_MODELS_DIR", BASE_DIR))

ModelKind = Literal["pickle"]

@dataclass(frozen=True)
class ModelSpec:
    path: Path
    kind: ModelKind
    features: Optional[list[str]] = None  # None = infer/sidecar; list = explicit
    note: str = ""

# XGB (already working for you; keep as-is if you want)
XGB_FEATURES = [
    "BertzCT", "PEOE_VSA1", "HallKierAlpha", "SlogP_VSA2", "MinAbsEStateIndex",
    "EState_VSA4", "PEOE_VSA6", "EState_VSA3", "PEOE_VSA7", "BalabanJ", "PEOE_VSA9",
]

MODELS: dict[str, ModelSpec] = {
    "XGB": ModelSpec(
        path=MODELS_DIR / "xgb_model.pkl",
        kind="pickle",
        features=XGB_FEATURES,
        note="XGBoost (sklearn wrapper) pickled – 11 descriptors",
    ),
    "Ensamble": ModelSpec(
        path=MODELS_DIR / "ensamble_model.pkl",
        kind="pickle",
        features=[
            "MaxEStateIndex",
            "qed",
            "AvgIpc",
            "BertzCT",
            "Chi4n",
            "PEOE_VSA1",
            "PEOE_VSA11",
            "PEOE_VSA3",
            "PEOE_VSA6",
            "PEOE_VSA9",
            "SMR_VSA6",
            "EState_VSA2",
            "EState_VSA3",
            "EState_VSA4",
            "VSA_EState2",
            "VSA_EState4",
            "VSA_EState6",
            "NumHeteroatoms",
            "fr_NH1",
            "fr_NH2",
            "fr_amide",
        ],
        note="Sklearn ensemble pickled – 21 descriptors (explicit order)",
    ),
}