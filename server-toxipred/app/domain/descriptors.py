from typing import Dict
from fastapi import HTTPException
from rdkit import Chem
from rdkit.Chem import QED, Descriptors, rdMolDescriptors, Fragments, rdPartialCharges

def descriptor_map(mol: Chem.Mol) -> Dict[str, float]:
    Chem.SanitizeMol(mol)
    try: rdPartialCharges.ComputeGasteigerCharges(mol)
    except Exception: pass
    d: Dict[str, float] = {}

    def add(name, fn):
        try: d[name] = float(fn(mol))
        except Exception: d[name] = float("nan")

    # include the ones your models need (and a few common ones)
    add("MaxEStateIndex", Descriptors.MaxEStateIndex)
    add("qed", QED.qed)
    add("AvgIpc", Descriptors.AvgIpc)
    add("BertzCT", Descriptors.BertzCT)
    add("Chi4n", Descriptors.Chi4n)
    for i in (1,11,3,6,9):
        add(f"PEOE_VSA{i}", getattr(Descriptors, f"PEOE_VSA{i}"))
    add("SMR_VSA6", Descriptors.SMR_VSA6)
    for i in (2,3,4):
        add(f"EState_VSA{i}", getattr(Descriptors, f"EState_VSA{i}"))
    for i in (2,4,6):
        add(f"VSA_EState{i}", getattr(Descriptors, f"VSA_EState{i}"))
    add("NumHeteroatoms", rdMolDescriptors.CalcNumHeteroatoms)
    add("fr_NH1", Fragments.fr_NH1)
    add("fr_NH2", Fragments.fr_NH2)
    add("fr_amide", Fragments.fr_amide)

    # extras your XGB uses
    add("HallKierAlpha", Descriptors.HallKierAlpha)
    add("SlogP_VSA2", Descriptors.SlogP_VSA2)
    add("MinAbsEStateIndex", Descriptors.MinAbsEStateIndex)
    add("PEOE_VSA7", Descriptors.PEOE_VSA7)
    add("BalabanJ", Descriptors.BalabanJ)

    return d