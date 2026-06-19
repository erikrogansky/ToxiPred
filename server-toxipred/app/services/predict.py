from typing import Tuple, Any, Optional

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
    return y0, conf
