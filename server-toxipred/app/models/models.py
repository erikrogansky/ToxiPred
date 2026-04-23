from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal, Optional
import os

BASE_DIR = Path(__file__).resolve().parent
MODELS_DIR = Path(os.getenv("TOXIPRED_MODELS_DIR", BASE_DIR))

ModelKind = Literal["pickle"]
ExplainerKind = Literal["tree", "linear", "deep", "kernel"]
FeatureSource = Literal["descriptors", "mixed"]
TestType = Literal["in_vitro", "in_vivo", "in_chemico"]
PredictionTarget = Literal["corrosion", "photo_irritation", "photo_toxicity"]

@dataclass(frozen=True)
class ModelSpec:
    path: Path
    kind: ModelKind
    features: Optional[list[str]] = None
    note: str = ""
    explainer: Optional[ExplainerKind] = None
    feature_source: FeatureSource = "descriptors"
    test_type: Optional[TestType] = None
    prediction_target: Optional[PredictionTarget] = None
    positive_label: str = "Positive"
    negative_label: str = "Negative"
    classification_threshold: Optional[float] = None
    dataset: Optional[str] = None

# ── Feature lists (from training CSV headers) ─────────────────────────

_FEATURES_CORR = [
    "EState_VSA6", "PEOE_VSA9", "MolWt", "MaxAbsEStateIndex",
    "PEOE_VSA10", "Ipc", "FpDensityMorgan1", "SlogP_VSA1",
    "MaxPartialCharge", "fr_Al_OH", "BCUT2D_MWLOW", "Kappa3",
    "MinEStateIndex", "EState_VSA8", "VSA_EState5", "PEOE_VSA1",
    "FractionCSP3", "EState_VSA7", "MinPartialCharge", "SMR_VSA10",
    "SPS", "pka", "EState_VSA9",
]

_FEATURES_CHEMICO: Optional[list[str]] = None
_FEATURES_3T3: Optional[list[str]] = None
_FEATURES_GB_IN_VIVO_CORR = [
    "MinEStateIndex",
    "MolWt",
    "BCUT2D_MWLOW",
    "BertzCT",
    "Ipc",
    "Kappa2",
]

def _load_json_features(name: str) -> list[str]:
    import json
    candidates = [
        MODELS_DIR / "XGB" / f"{name}.features.json",
        MODELS_DIR / f"{name}.features.json",  # legacy location
    ]
    for p in candidates:
        if p.exists():
            with open(p) as f:
                return json.load(f)
    raise FileNotFoundError(f"Missing features sidecar for '{name}': tried {', '.join(str(c) for c in candidates)}")

def _chemico_features() -> list[str]:
    global _FEATURES_CHEMICO
    if _FEATURES_CHEMICO is None:
        _FEATURES_CHEMICO = _load_json_features("xgb_in_chemico_photo")
    return _FEATURES_CHEMICO

def _3t3_features() -> list[str]:
    global _FEATURES_3T3
    if _FEATURES_3T3 is None:
        _FEATURES_3T3 = _load_json_features("xgb_in_vitro_3T3_photo")
    return _FEATURES_3T3

# ── Model registry ────────────────────────────────────────────────────

MODELS: dict[str, ModelSpec] = {
    "XGB Corrosion": ModelSpec(
        path=MODELS_DIR / "XGB/xgb_in_vitro_corr.pkl",
        kind="pickle",
        features=_FEATURES_CORR,
        note="XGBoost – 23 descriptors, in vitro corrosion",
        explainer="tree",
        feature_source="descriptors",
        test_type="in_vitro",
        prediction_target="corrosion",
        positive_label="Corrosive",
        negative_label="Non-corrosive",
        classification_threshold=0.5,
        dataset="In Vitro (151 compounds)",
    ),
    "XGB Phototox Chemico": ModelSpec(
        path=MODELS_DIR / "XGB/xgb_in_chemico_photo.pkl",
        kind="pickle",
        note="XGBoost – 48 descriptors + MACCS, in chemico phototoxicity",
        explainer="tree",
        feature_source="mixed",
        test_type="in_chemico",
        prediction_target="photo_toxicity",
        positive_label="Phototoxic",
        negative_label="Non-phototoxic",
        classification_threshold=0.45,
        dataset="In Chemico (162 compounds)",
    ),
    "XGB Phototox 3T3": ModelSpec(
        path=MODELS_DIR / "XGB/xgb_in_vitro_3T3_photo.pkl",
        kind="pickle",
        note="XGBoost – 52 descriptors + AtomPair FP, in vitro 3T3 phototoxicity",
        explainer="tree",
        feature_source="mixed",
        test_type="in_vitro",
        prediction_target="photo_toxicity",
        positive_label="Phototoxic",
        negative_label="Non-phototoxic",
        classification_threshold=0.5,
        dataset="In Vitro 3T3 (396 compounds)",
    ),
    "GB Corrosion (in vivo)": ModelSpec(
        path=MODELS_DIR / "GB/gb_in_vivo_corrosion.pkl",
        kind="pickle",
        features=_FEATURES_GB_IN_VIVO_CORR,
        note="GradientBoosting – 6 descriptors, in vivo corrosion",
        explainer="tree",
        feature_source="descriptors",
        test_type="in_vivo",
        prediction_target="corrosion",
        positive_label="Corrosive",
        negative_label="Non-corrosive",
        classification_threshold=0.5,
        dataset="In Vivo (189 compounds)",
    ),
}
