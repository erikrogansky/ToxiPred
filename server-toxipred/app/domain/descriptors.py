"""Unified descriptor computation for all toxicity models.

Feature groups computed on demand:
  - RDKit standard descriptors (EState, PEOE_VSA, BCUT2D, fr_*, etc.)
  - Mordred autocorrelation (GATS, MATS, ATS) and Chi connectivity
  - Gasteiger charge summaries (Q-descriptors)
  - HOMO / LUMO via RDKit Extended Hückel Theory
  - pKa rough estimate (functional-group based)
  - SPS (Wiener index approximation)
  - MACCS fingerprint bits
  - Hashed AtomPair fingerprint bits
"""

from __future__ import annotations

import logging
import re
from typing import Dict, List, Optional, Set

import numpy as np
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    MACCSkeys,
    rdEHTTools,
    rdMolDescriptors,
    rdPartialCharges,
    rdmolops,
)

logger = logging.getLogger(__name__)

# ── mordred (lazy singleton) ──────────────────────────────────────────

_mordred_calc = None


def _get_mordred_calc():
    global _mordred_calc
    if _mordred_calc is None:
        from mordred import Calculator
        from mordred.Autocorrelation import ATS, AATS, ATSC, AATSC, MATS, GATS
        from mordred.Chi import Chi

        _mordred_calc = Calculator()
        for D in (ATS, AATS, ATSC, AATSC, MATS, GATS, Chi):
            _mordred_calc.register(D)
    return _mordred_calc


# ── name mapping: mordred → CSV convention ────────────────────────────

_AC_RE = re.compile(r"^(GATS|MATS|ATS|AATS|ATSC|AATSC)(\d+)([a-z]+)$")


def _mordred_to_csv(name: str) -> str:
    """GATS4e → GATSe4, MATS7p → MATSp7, etc."""
    m = _AC_RE.match(name)
    if m:
        prefix, lag, prop = m.groups()
        return f"{prefix}{prop}{lag}"
    return name


# ── feature-group helpers ─────────────────────────────────────────────


def _rdkit_descriptors(mol: Chem.Mol) -> Dict[str, float]:
    d = {}
    for k, v in Descriptors.CalcMolDescriptors(mol).items():
        try:
            d[k] = float(v) if v is not None else float("nan")
        except (TypeError, ValueError):
            d[k] = float("nan")
    # Extras that may not be registered in CalcMolDescriptors
    for name, fn in (
        ("NumAmideBonds", rdMolDescriptors.CalcNumAmideBonds),
        ("NumSaturatedRings", rdMolDescriptors.CalcNumSaturatedRings),
    ):
        if name not in d:
            try:
                d[name] = float(fn(mol))
            except Exception:
                d[name] = float("nan")
    return d


def _has_dummy_atoms(mol: Chem.Mol) -> bool:
    return any(a.GetAtomicNum() == 0 for a in mol.GetAtoms())


def _mordred_descriptors(mol: Chem.Mol) -> Dict[str, float]:
    if _has_dummy_atoms(mol):
        return {}
    calc = _get_mordred_calc()
    result = calc(mol)
    d: Dict[str, float] = {}
    for desc, val in zip(calc.descriptors, result):
        name = str(desc)
        try:
            fval = float(val)
        except (ValueError, TypeError):
            fval = float("nan")
        d[name] = fval
        csv_name = _mordred_to_csv(name)
        if csv_name != name:
            d[csv_name] = fval
    return d


def _charge_descriptors(mol: Chem.Mol) -> Dict[str, float]:
    """Gasteiger-charge Q-descriptors (QNmax, QCmin, QHss, …)."""
    try:
        rdPartialCharges.ComputeGasteigerCharges(mol)
    except Exception:
        return {}

    by_elem: Dict[str, list] = {}
    all_q: list[float] = []
    for atom in mol.GetAtoms():
        try:
            q = atom.GetDoubleProp("_GasteigerCharge")
            if np.isnan(q) or np.isinf(q):
                continue
        except Exception:
            continue
        all_q.append(q)
        by_elem.setdefault(atom.GetSymbol(), []).append(q)

    d: Dict[str, float] = {}
    if all_q:
        d["Qmax"] = max(all_q)
        d["Qmin"] = min(all_q)
        d["Qass"] = sum(abs(v) for v in all_q)

    for sym, prefix in (("N", "QN"), ("C", "QC"), ("H", "QH"), ("O", "QO"), ("S", "QS")):
        vals = by_elem.get(sym, [])
        if vals:
            d[f"{prefix}max"] = max(vals)
            d[f"{prefix}min"] = min(vals)
            d[f"{prefix}ss"] = sum(v ** 2 for v in vals)
        else:
            d[f"{prefix}max"] = 0.0
            d[f"{prefix}min"] = 0.0
            d[f"{prefix}ss"] = 0.0
    return d


def _homo_lumo(mol: Chem.Mol) -> Dict[str, float]:
    """HOMO / LUMO energies via Extended Hückel Theory (RDKit built-in)."""
    d: Dict[str, float] = {
        "HOMO_eV": float("nan"),
        "LUMO_eV": float("nan"),
        "HL_Gap_eV": float("nan"),
    }
    if _has_dummy_atoms(mol):
        return d
    try:
        mol3d = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        if AllChem.EmbedMolecule(mol3d, params) != 0:
            if AllChem.EmbedMolecule(mol3d, AllChem.ETKDG()) != 0:
                return d
        try:
            AllChem.MMFFOptimizeMolecule(mol3d, maxIters=500)
        except Exception:
            try:
                AllChem.UFFOptimizeMolecule(mol3d, maxIters=500)
            except Exception:
                pass
        ok, results = rdEHTTools.RunMol(mol3d)
        if ok:
            orb = results.GetOrbitalEnergies()
            n_e = sum(a.GetAtomicNum() for a in mol3d.GetAtoms())
            homo_i = n_e // 2 - 1
            lumo_i = homo_i + 1
            if 0 <= homo_i < len(orb) and lumo_i < len(orb):
                d["HOMO_eV"] = float(orb[homo_i])
                d["LUMO_eV"] = float(orb[lumo_i])
                d["HL_Gap_eV"] = d["LUMO_eV"] - d["HOMO_eV"]
    except Exception:
        logger.debug("EHT computation failed", exc_info=True)
    return d


def _sps(mol: Chem.Mol) -> float:
    """Sum of Path Scores ≈ Wiener index (half-sum of distance matrix)."""
    try:
        dm = rdmolops.GetDistanceMatrix(mol)
        return float(dm.sum() / 2.0)
    except Exception:
        return float("nan")


_PKA_SMARTS: list[tuple[str, float]] = [
    ("S(=O)(=O)O", 1.5),
    ("P(=O)(O)O",  2.0),
    ("C(=O)O",     4.0),
    ("c1ccccc1O", 10.0),
    ("CO",        16.0),
]


def _pka_estimate(mol: Chem.Mol) -> float:
    """Very rough pKa of the most acidic functional group."""
    best = float("nan")
    for sma, pka_val in _PKA_SMARTS:
        pat = Chem.MolFromSmarts(sma)
        if pat and mol.HasSubstructMatch(pat):
            if np.isnan(best) or pka_val < best:
                best = pka_val
    return best


def _maccs_fingerprint(mol: Chem.Mol) -> Dict[str, float]:
    fp = MACCSkeys.GenMACCSKeys(mol)
    return {f"MACCS_{i}": float(fp[i]) for i in range(167)}


def _atompair_fingerprint(mol: Chem.Mol, n_bits: int = 512) -> Dict[str, float]:
    fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol, nBits=n_bits)
    return {f"AtomPair_{i}": float(fp[i]) for i in range(n_bits)}


# ── need-detection helpers ────────────────────────────────────────────

_MORDRED_PREFIXES = ("GATS", "MATS", "ATS", "AATS", "ATSC", "AATSC")


def _needs_mordred(wanted: Set[str]) -> bool:
    for n in wanted:
        if any(n.startswith(p) for p in _MORDRED_PREFIXES):
            return True
        if re.match(r"^Chiv?\d+(pc|ch|c)$", n):
            return True
    return False


# ── public API ────────────────────────────────────────────────────────


def descriptor_map(
    mol: Chem.Mol,
    feature_names: Optional[List[str]] = None,
) -> Dict[str, float]:
    """Compute all molecular features a model may need, keyed by name.

    When *feature_names* is given, expensive groups (mordred, EHT, FP)
    are only computed if at least one requested name belongs to that group.
    """
    Chem.SanitizeMol(mol)
    try:
        rdPartialCharges.ComputeGasteigerCharges(mol)
    except Exception:
        pass

    wanted: Optional[Set[str]] = set(feature_names) if feature_names else None
    d: Dict[str, float] = {}

    # 1. RDKit standard descriptors (fast – always run)
    d.update(_rdkit_descriptors(mol))

    # 2. Mordred autocorrelation + Chi variants (only if needed)
    if wanted is None or _needs_mordred(wanted):
        d.update(_mordred_descriptors(mol))

    # 3. Q-charge descriptors
    if wanted is None or any(n.startswith("Q") and n[1:2].isupper() for n in (wanted or set())):
        d.update(_charge_descriptors(mol))

    # 4. HOMO / LUMO
    if wanted is None or (wanted & {"LUMO_eV", "HL_Gap_eV", "HOMO_eV"}):
        d.update(_homo_lumo(mol))

    # 5. SPS
    if wanted is None or "SPS" in (wanted or set()):
        d["SPS"] = _sps(mol)

    # 6. pKa
    if wanted is None or "pka" in (wanted or set()):
        d["pka"] = _pka_estimate(mol)

    # 7. MACCS fingerprint
    if wanted is None or any(n.startswith("MACCS_") for n in (wanted or set())):
        d.update(_maccs_fingerprint(mol))

    # 8. AtomPair fingerprint
    if wanted is None or any(n.startswith("AtomPair_") for n in (wanted or set())):
        d.update(_atompair_fingerprint(mol))

    return d