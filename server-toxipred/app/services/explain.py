# app/services/explain.py
from __future__ import annotations
from typing import Callable, Iterable, List, Tuple, Optional
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdmolops

# --- Score helper: returns P(positive) if available, else a squashed score ---
def _score_fn(model) -> Callable[[np.ndarray], np.ndarray]:
    if hasattr(model, "predict_proba"):
        return lambda X: np.asarray(model.predict_proba(X))[:, 1]
    if hasattr(model, "decision_function"):
        def f(X):
            z = np.asarray(model.decision_function(X)).reshape(-1)
            return 1.0 / (1.0 + np.exp(-z))
        return f
    return lambda X: np.asarray(model.predict(X)).reshape(-1).astype(float)

# --- Safe sanitize (don’t bomb out if masking breaks valence heuristics) ---
def _sanitize(m: Chem.Mol) -> Chem.Mol:
    try:
        Chem.SanitizeMol(m)
    except Exception:
        pass
    return m

# --- Masking primitives ------------------------------------------------------

def _mask_atoms_with_dummy(mol: Chem.Mol, atom_idxs: Iterable[int]) -> Chem.Mol:
    """Replace selected atoms with a dummy [*] to keep topology mostly intact."""
    rw = Chem.RWMol(mol)
    atom_idxs = set(atom_idxs)
    for idx in atom_idxs:
        a = rw.GetAtomWithIdx(idx)
        a.SetAtomicNum(0)      # 0 => dummy
        a.SetFormalCharge(0)
        a.SetIsAromatic(False)
        a.SetNoImplicit(True)  # avoid auto H completion
    return _sanitize(rw.GetMol())

def _ring_atom_groups(mol: Chem.Mol) -> List[List[int]]:
    rings = mol.GetRingInfo().AtomRings()
    return [list(r) for r in rings]

# --- Descriptor wrapper -------------------------------------------------------

def _featurize(descriptor_fn: Callable[[Chem.Mol], np.ndarray],
               mol: Chem.Mol) -> np.ndarray:
    """Assumes your descriptor_fn returns 1D array for a single mol."""
    X = descriptor_fn(mol)
    X = np.asarray(X).reshape(1, -1)
    return X

# --- Core perturbation engine -------------------------------------------------

def _explain_by_groups(
    model,
    mol: Chem.Mol,
    descriptor_fn: Callable[[Chem.Mol], np.ndarray],
    groups: List[List[int]],
    batch_size: int = 64,
) -> Tuple[np.ndarray, float]:
    """
    Returns (attr_per_group, base_pred). Attribution = base_pred - masked_pred.
    """
    scorer = _score_fn(model)
    # Base prediction on the unmasked molecule
    X0 = _featurize(descriptor_fn, mol)
    base = float(scorer(X0)[0])

    # Build masked molecules in batches
    masked_preds: List[float] = []
    for i in range(0, len(groups), batch_size):
        batch = groups[i:i+batch_size]
        masked = []
        for g in batch:
            mg = _mask_atoms_with_dummy(mol, g)
            masked.append(_featurize(descriptor_fn, mg))
        if not masked:
            continue
        Xb = np.vstack(masked)  # (B, D)
        preds = scorer(Xb)      # (B,)
        masked_preds.extend([float(p) for p in preds])

    masked_preds = np.asarray(masked_preds, dtype=float)
    atts = base - masked_preds  # positive => group increases risk
    return atts, base

# --- Public APIs --------------------------------------------------------------

def explain_atom_importance(
    model,
    mol: Chem.Mol,
    descriptor_fn: Callable[[Chem.Mol], np.ndarray],
    max_atoms: int = 120,
    normalize: bool = False,
) -> Tuple[np.ndarray, float]:
    """
    Per-atom attribution via single-atom masking.
    Returns (atom_importances, base_pred).
    """
    n = mol.GetNumAtoms()
    if n == 0:
        return np.array([]), float("nan")
    if n > max_atoms:
        # Light sampling for very large molecules to keep runtime sane
        idxs = list(range(n))
    else:
        idxs = list(range(n))
    groups = [[i] for i in idxs]
    atts, base = _explain_by_groups(model, mol, descriptor_fn, groups)
    # If we sampled, stitch back (here we didn’t subsample, but hook left in)
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
    """
    Per-ring attribution via masking all atoms in each ring.
    Returns (ring_importances, base_pred, rings_as_atom_idx_lists).
    """
    rings = _ring_atom_groups(mol)
    if not rings:
        return np.array([]), float("nan"), []
    atts, base = _explain_by_groups(model, mol, descriptor_fn, rings)
    if normalize and atts.std() > 0:
        atts = (atts - atts.mean()) / atts.std()
    return atts, base, rings

# Optional: group atoms by functional groups or SMARTS matches
def explain_smarts_groups(
    model,
    mol: Chem.Mol,
    descriptor_fn: Callable[[Chem.Mol], np.ndarray],
    smarts_list: List[str],
    normalize: bool = False,
) -> Tuple[np.ndarray, float, List[List[int]]]:
    """
    Attribution for arbitrary SMARTS-defined groups.
    """
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
