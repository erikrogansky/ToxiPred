from typing import Any, List, Optional, Tuple
from fastapi import HTTPException
import json
import pickle as _pickle
from pathlib import Path
from app.models.models import ModelSpec

_cache: dict[str, Any] = {}

# ── Cross-version compatibility unpickler ───────────────────────────────
#
# numpy 2.0 renamed the internal submodule ``numpy.core`` → ``numpy._core``.
# Pickles saved with numpy >= 2 therefore reference ``numpy._core.multiarray``
# (and similar), which does not exist in numpy 1.x. When such a pickle is
# loaded under numpy 1.26 the stdlib unpickler raises ``ModuleNotFoundError``
# which joblib propagates; the ``dill`` fallback then surfaces the misleading
# "STACK_GLOBAL requires str" error.
#
# This remap transparently rewrites the module path so the pickle can load
# while we stay pinned to numpy 1.26.x. It is a no-op under numpy >= 2.

_NUMPY_MODULE_REMAP: tuple[tuple[str, str], ...] = (
    ("numpy._core", "numpy.core"),
    ("numpy.core",  "numpy._core"),
)


class _CompatUnpickler(_pickle.Unpickler):
    def find_class(self, module: str, name: str):
        try:
            return super().find_class(module, name)
        except (ModuleNotFoundError, AttributeError):
            for src, dst in _NUMPY_MODULE_REMAP:
                if module == src or module.startswith(src + "."):
                    remapped = dst + module[len(src):]
                    try:
                        return super().find_class(remapped, name)
                    except (ModuleNotFoundError, AttributeError):
                        continue
            raise


def _load_with_compat_unpickler(path: Path):
    with open(path, "rb") as f:
        return _CompatUnpickler(f).load()


def load_pickle(spec: ModelSpec):
    if not spec.path.exists():
        raise HTTPException(500, f"Model file not found: {spec.path}")

    errors: list[str] = []

    try:
        import joblib
        return joblib.load(spec.path)
    except Exception as e:
        errors.append(f"joblib: {e!r}")

    try:
        return _load_with_compat_unpickler(spec.path)
    except Exception as e:
        errors.append(f"pickle(compat): {e!r}")

    try:
        import dill
        with open(spec.path, "rb") as f:
            return dill.load(f)
    except Exception as e:
        errors.append(f"dill: {e!r}")

    raise HTTPException(
        500,
        f"Failed to load pickle '{spec.path.name}': " + " | ".join(errors),
    )


def get_model(name: str, specs: dict[str, ModelSpec]):
    spec = specs.get(name)
    if not spec: raise HTTPException(404, f"Unknown model '{name}'")
    if name not in _cache: _cache[name] = load_pickle(spec)
    return _cache[name]

def sidecar_path(spec: ModelSpec) -> Path:
    return spec.path.with_suffix(spec.path.suffix + ".features.json")

def load_sidecar_features(spec: ModelSpec) -> Optional[List[str]]:
    p = sidecar_path(spec)
    if p.exists():
        try:
            data = json.load(open(p, "r"))
            if isinstance(data, list) and all(isinstance(x, str) for x in data):
                return data
        except Exception: ...
    return None

def infer_feature_names(model: Any) -> Tuple[Optional[List[str]], Optional[int]]:
    names = None; n_in = getattr(model, "n_features_in_", None)
    if hasattr(model, "feature_names_in_"):
        try: names = [str(x) for x in model.feature_names_in_]
        except Exception: ...
    return names, (int(n_in) if n_in is not None else None)

def resolve_feature_names(spec: ModelSpec, model: Any) -> List[str]:
    if spec.features: return spec.features
    side = load_sidecar_features(spec)
    if side: return side
    names, n_in = infer_feature_names(model)
    if names: return names
    if n_in:
        raise HTTPException(500, f"Model reports n_features_in_={n_in} but no names. "
                                 f"Provide {sidecar_path(spec)} or set spec.features.")
    raise HTTPException(500, "Could not resolve feature names for the model.")
