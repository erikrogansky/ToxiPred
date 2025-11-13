import os
import re
import logging
from dataclasses import dataclass
from typing import Optional

import requests
from requests import RequestException
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

log = logging.getLogger(__name__)

CAS_RE = re.compile(r"^\d{2,7}-\d{2}-\d$")
PUBCHEM_ENABLED = os.getenv("PUBCHEM_LOOKUP", "1") == "1"

@dataclass
class ResolvedChem:
    smiles: str
    mol: Chem.Mol
    name: Optional[str]
    formula: str
    source: str 
    query: str


def _from_smiles(query: str) -> Optional[ResolvedChem]:
    m = Chem.MolFromSmiles(query)
    if m is None:
        return None
    can = Chem.MolToSmiles(m, canonical=True)
    formula = rdMolDescriptors.CalcMolFormula(m)
    return ResolvedChem(
        smiles=can,
        mol=m,
        name=None,
        formula=formula,
        source="smiles",
        query=query,
    )


def _from_inchi(query: str) -> Optional[ResolvedChem]:
    if not query.startswith("InChI="):
        return None
    m = Chem.MolFromInchi(query)
    if m is None:
        return None
    can = Chem.MolToSmiles(m, canonical=True)
    formula = rdMolDescriptors.CalcMolFormula(m)
    return ResolvedChem(
        smiles=can,
        mol=m,
        name=None,
        formula=formula,
        source="inchi",
        query=query,
    )


def _from_pubchem_name_or_cas(query: str) -> Optional[ResolvedChem]:
    if not PUBCHEM_ENABLED:
        return None

    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
        f"{requests.utils.quote(query)}/property/CanonicalSMILES,"
        "MolecularFormula,IUPACName,IsomericSMILES,SMILES,ConnectivitySMILES/JSON"
    )

    try:
        resp = requests.get(url, timeout=5)
        resp.raise_for_status()
    except RequestException as e:
        log.warning("PubChem lookup failed for %r: %s", query, e)
        return None

    try:
        props = resp.json()["PropertyTable"]["Properties"][0]

        smiles = (
            props.get("CanonicalSMILES")
            or props.get("IsomericSMILES")
            or props.get("SMILES")
            or props.get("ConnectivitySMILES")
        )

        if not smiles:
            log.warning(
                "PubChem response for %r has no usable SMILES field: %r",
                query,
                props,
            )
            return None

        formula = props.get("MolecularFormula")
        name = props.get("IUPACName") or query

        m = Chem.MolFromSmiles(smiles)
        if m is None:
            log.warning("RDKit failed to parse SMILES from PubChem for %r: %r", query, smiles)
            return None

        can = Chem.MolToSmiles(m, canonical=True)
        if not formula:
            formula = rdMolDescriptors.CalcMolFormula(m)

        return ResolvedChem(
            smiles=can,
            mol=m,
            name=name,
            formula=formula,
            source="cas" if CAS_RE.match(query) else "name",
            query=query,
        )
    except Exception as e:
        log.warning("PubChem JSON parsing failed for %r: %s", query, e)
        return None



def resolve_chemical(query: str) -> Optional[ResolvedChem]:
    query = query.strip()

    if query.startswith("InChI="):
        res = _from_inchi(query)
        if res:
            return res

    res = _from_smiles(query)
    if res:
        return res

    if CAS_RE.match(query) or query.replace(" ", "").isalpha():
        res = _from_pubchem_name_or_cas(query)
        if res:
            return res

    res = _from_pubchem_name_or_cas(query)
    if res:
        return res

    return None


