from typing import Tuple, Any

def predict_vectorized(model: Any, X) -> Tuple[float | int | str, float | None]:
    if not hasattr(model, "predict"):
        raise ValueError("Object has no predict()")
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
