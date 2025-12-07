from dataclasses import dataclass
from pathlib import Path
from typing import Literal, Optional
import os

BASE_DIR = Path(__file__).resolve().parent
MODELS_DIR = Path(os.getenv("TOXIPRED_MODELS_DIR", BASE_DIR))

ModelKind = Literal["pickle"]
ExplainerKind = Literal["tree", "linear", "deep", "kernel"]
FeatureSource = Literal["descriptors", "morgan", "tokens"]
TestType = Literal["in_vitro", "in_vivo", "in_chemico"]
PredictionTarget = Literal["photo_irritation", "photo_toxicity"]

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

XGB_FEATURES = [
    "BertzCT", "PEOE_VSA1", "HallKierAlpha", "SlogP_VSA2", "MinAbsEStateIndex",
    "EState_VSA4", "PEOE_VSA6", "EState_VSA3", "PEOE_VSA7", "BalabanJ", "PEOE_VSA9",
]

MODELS: dict[str, ModelSpec] = {
    "XGBoost": ModelSpec(
        path=MODELS_DIR / "xgb_model.pkl",
        kind="pickle",
        features=XGB_FEATURES,
        note="XGBoost (sklearn wrapper) pickled – 11 descriptors",
        explainer="tree",
        feature_source="descriptors",
        test_type="in_vitro",
        prediction_target="photo_toxicity",
        positive_label="Phototoxic",
        negative_label="Non-phototoxic",
    ),
    "Ensamble": ModelSpec(
        path=MODELS_DIR / "ensamble_model.pkl",
        kind="pickle",
        features=[
            "MaxEStateIndex","qed","AvgIpc","BertzCT","Chi4n",
            "PEOE_VSA1","PEOE_VSA11","PEOE_VSA3","PEOE_VSA6","PEOE_VSA9",
            "SMR_VSA6","EState_VSA2","EState_VSA3","EState_VSA4",
            "VSA_EState2","VSA_EState4","VSA_EState6",
            "NumHeteroatoms","fr_NH1","fr_NH2","fr_amide",
        ],
        note="Sklearn ensemble pickled – 21 descriptors (explicit order)",
        explainer=None,
        feature_source="descriptors",
        test_type="in_chemico",
        prediction_target="photo_irritation",
        positive_label="Irritant",
        negative_label="Non-irritant",
    ),
}
