from typing import Tuple, Any, Optional
import numpy as np


def _sigmoid(x: float) -> float:
    x = float(np.clip(x, -709, 709))
    return float(1.0 / (1.0 + np.exp(-x)))


def _decision_confidence(model: Any, X) -> float | None:
    if not hasattr(model, "decision_function"):
        return None
    try:
        margin = np.asarray(model.decision_function(X), dtype=float)
    except Exception:
        return None

    if margin.ndim == 1 or margin.shape[-1] == 1:
        p1 = _sigmoid(margin.reshape(-1)[0])
        return max(p1, 1.0 - p1)

    if margin.ndim == 2 and margin.shape[0] > 0:
        row = margin[0]
        row = row - np.max(row)
        exp = np.exp(row)
        denom = float(exp.sum())
        if denom > 0:
            return float(exp.max() / denom)

    return None

def predict_vectorized(
    model: Any,
    X,
    threshold: Optional[float] = None,
) -> Tuple[float | int | str, float | None]:
    if not hasattr(model, "predict"):
        raise ValueError("Object has no predict()")

    # Custom threshold: use predict_proba and threshold manually
    if threshold is not None and hasattr(model, "predict_proba"):
        proba = model.predict_proba(X)[0]
        conf = float(max(proba)) if hasattr(proba, "__iter__") else float(proba)
        y0 = int(proba[1] >= threshold) if len(proba) > 1 else int(proba[0] >= threshold)
        return y0, conf

    y = model.predict(X)
    y0 = y[0].item() if hasattr(y[0], "item") else y[0]
    conf = None
    if hasattr(model, "predict_proba"):
        try:
            p0 = model.predict_proba(X)[0]
            conf = float(max(p0)) if hasattr(p0, "__iter__") else float(p0)
        except Exception:
            pass
    if conf is None:
        conf = _decision_confidence(model, X)
    return y0, conf
