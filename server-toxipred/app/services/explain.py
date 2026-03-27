from __future__ import annotations
import logging
from typing import Callable, Iterable, List, Tuple, Optional
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdmolops

logger = logging.getLogger(__name__)

def _score_fn(model) -> Callable[[np.ndarray], np.ndarray]:
    if hasattr(model, "predict_proba"):
        return lambda X: np.asarray(model.predict_proba(X))[:, 1]
    if hasattr(model, "decision_function"):
        def f(X):
            z = np.asarray(model.decision_function(X)).reshape(-1)
            return 1.0 / (1.0 + np.exp(-z))
        return f
    return lambda X: np.asarray(model.predict(X)).reshape(-1).astype(float)

def _sanitize(m: Chem.Mol) -> Chem.Mol:
    try:
        Chem.SanitizeMol(m)
    except Exception:
        pass
    return m


def _mask_atoms_with_dummy(mol: Chem.Mol, atom_idxs: Iterable[int]) -> Chem.Mol:
    rw = Chem.RWMol(mol)
    atom_idxs = set(atom_idxs)
    for idx in atom_idxs:
        a = rw.GetAtomWithIdx(idx)
        a.SetAtomicNum(0)
        a.SetFormalCharge(0)
        a.SetIsAromatic(False)
        a.SetNoImplicit(True)
    return _sanitize(rw.GetMol())

def _ring_atom_groups(mol: Chem.Mol) -> List[List[int]]:
    rings = mol.GetRingInfo().AtomRings()
    return [list(r) for r in rings]

def _featurize(descriptor_fn: Callable[[Chem.Mol], np.ndarray],
               mol: Chem.Mol) -> np.ndarray:
    """Assumes your descriptor_fn returns 1D array for a single mol."""
    X = descriptor_fn(mol)
    X = np.asarray(X).reshape(1, -1)
    return X

def _explain_by_groups(
    model,
    mol: Chem.Mol,
    descriptor_fn: Callable[[Chem.Mol], np.ndarray],
    groups: List[List[int]],
    batch_size: int = 64,
) -> Tuple[np.ndarray, float]:
    scorer = _score_fn(model)
    X0 = _featurize(descriptor_fn, mol)
    base = float(scorer(X0)[0])

    masked_preds: List[float] = []
    for g in groups:
        try:
            mg = _mask_atoms_with_dummy(mol, g)
            Xg = _featurize(descriptor_fn, mg)
            pred = float(scorer(Xg)[0])
        except Exception:
            logger.debug("masking failed for atom group %s, scoring as base", g)
            pred = base
        masked_preds.append(pred)

    masked_preds = np.asarray(masked_preds, dtype=float)
    atts = base - masked_preds
    return atts, base


def explain_atom_importance(
    model,
    mol: Chem.Mol,
    descriptor_fn: Callable[[Chem.Mol], np.ndarray],
    max_atoms: int = 120,
    normalize: bool = False,
) -> Tuple[np.ndarray, float]:
    n = mol.GetNumAtoms()
    if n == 0:
        return np.array([]), float("nan")
    if n > max_atoms:
        idxs = list(range(n))
    else:
        idxs = list(range(n))
    groups = [[i] for i in idxs]
    atts, base = _explain_by_groups(model, mol, descriptor_fn, groups)
    atom_atts = np.zeros(n, dtype=float)
    for i, aidx in enumerate(idxs):
        atom_atts[aidx] = atts[i]
    if normalize and atom_atts.std() > 0:
        atom_atts = (atom_atts - atom_atts.mean()) / atom_atts.std()
    return atom_atts, base

def explain_ring_importance(
    model,
    mol: Chem.Mol,
    descriptor_fn: Callable[[Chem.Mol], np.ndarray],
    normalize: bool = False,
) -> Tuple[np.ndarray, float, List[List[int]]]:
    rings = _ring_atom_groups(mol)
    if not rings:
        return np.array([]), float("nan"), []
    atts, base = _explain_by_groups(model, mol, descriptor_fn, rings)
    if normalize and atts.std() > 0:
        atts = (atts - atts.mean()) / atts.std()
    return atts, base, rings

def explain_smarts_groups(
    model,
    mol: Chem.Mol,
    descriptor_fn: Callable[[Chem.Mol], np.ndarray],
    smarts_list: List[str],
    normalize: bool = False,
) -> Tuple[np.ndarray, float, List[List[int]]]:
    groups: List[List[int]] = []
    for s in smarts_list:
        pat = Chem.MolFromSmarts(s)
        if not pat:
            continue
        for match in mol.GetSubstructMatches(pat):
            groups.append(list(match))
    if not groups:
        return np.array([]), float("nan"), []
    atts, base = _explain_by_groups(model, mol, descriptor_fn, groups)
    if normalize and atts.std() > 0:
        atts = (atts - atts.mean()) / atts.std()
    return atts, base, groups

def explain_feature_importance(
    model,
    x: np.ndarray,
    baseline: Optional[np.ndarray] = None,
    normalize: bool = False,
) -> Tuple[np.ndarray, float]:
    x = np.asarray(x, dtype=float).reshape(1, -1)
    n_features = x.shape[1]

    scorer = _score_fn(model)
    base = float(scorer(x)[0])

    if baseline is not None:
        baseline_vec = np.asarray(baseline, dtype=float).reshape(1, -1)
        if baseline_vec.shape[1] != n_features:
            raise ValueError("baseline must have same number of features as x")
    else:
        baseline_vec = None

    scores = np.zeros(n_features, dtype=float)

    for j in range(n_features):
        x_pert = x.copy()
        if baseline_vec is not None:
            x_pert[0, j] = baseline_vec[0, j]
        else:
            x_pert[0, j] = 0.0

        pred_pert = float(scorer(x_pert)[0])
        scores[j] = base - pred_pert

    if normalize and scores.std() > 0:
        scores = (scores - scores.mean()) / scores.std()

    return scores, base


# ── SHAP-based explainers ────────────────────────────────────────────


def _is_xgboost(model) -> bool:
    """Check if the model is an XGBoost model (sklearn API or raw Booster)."""
    try:
        import xgboost
        return isinstance(model, (xgboost.XGBModel, xgboost.Booster))
    except ImportError:
        return False


def _xgboost_shap(model, x: np.ndarray) -> Tuple[np.ndarray, float]:
    """Use XGBoost's built-in TreeSHAP (pred_contribs).

    Returns SHAP values in log-odds (margin) space.
    The last column of contribs is the bias term (base_score).
    """
    import xgboost

    booster = model.get_booster() if hasattr(model, "get_booster") else model
    dm = xgboost.DMatrix(x, feature_names=booster.feature_names)
    contribs = booster.predict(dm, pred_contribs=True)  # shape (1, n_features+1)
    sv = contribs[0, :-1]
    base_value = float(contribs[0, -1])
    return sv.astype(float), base_value


def explain_feature_importance_shap(
    model,
    x: np.ndarray,
    explainer_kind: str = "tree",
    normalize: bool = False,
) -> Tuple[np.ndarray, float]:
    """Compute per-feature SHAP values.

    *explainer_kind* selects the SHAP explainer:
        "tree"   → XGBoost built-in or shap.TreeExplainer
        "linear" → shap.LinearExplainer
        "kernel" → shap.KernelExplainer (model-agnostic, slow)
    """
    x = np.asarray(x, dtype=float).reshape(1, -1)

    if explainer_kind == "tree" and _is_xgboost(model):
        sv, base_value = _xgboost_shap(model, x)
    elif explainer_kind == "tree":
        import shap
        explainer = shap.TreeExplainer(model)
        shap_values = explainer.shap_values(x)
        if isinstance(shap_values, list):
            sv = np.asarray(shap_values[1]).reshape(-1)
        elif shap_values.ndim == 3:
            sv = shap_values[0, :, 1]
        elif shap_values.ndim == 2:
            sv = shap_values[0]
        else:
            sv = shap_values.reshape(-1)
        base_value = float(
            explainer.expected_value[1]
            if isinstance(explainer.expected_value, (list, np.ndarray))
            and len(explainer.expected_value) > 1
            else explainer.expected_value
        )
    elif explainer_kind == "linear":
        import shap
        explainer = shap.LinearExplainer(model, x)
        shap_values = explainer.shap_values(x)
        sv = np.asarray(shap_values).reshape(-1)
        base_value = float(explainer.expected_value)
    elif explainer_kind == "kernel":
        import shap
        scorer = _score_fn(model)
        explainer = shap.KernelExplainer(scorer, x)
        shap_values = explainer.shap_values(x)
        sv = np.asarray(shap_values).reshape(-1)
        base_value = float(explainer.expected_value)
    else:
        raise ValueError(f"Unknown explainer kind: {explainer_kind}")

    if normalize and sv.std() > 0:
        sv = (sv - sv.mean()) / sv.std()

    return sv.astype(float), base_value
