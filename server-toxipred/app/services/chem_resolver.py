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
    trivial_name: Optional[str] = None
    other_names: Optional[list[str]] = None


def _get_cid_from_smiles(smiles: str) -> Optional[int]:
    """Get PubChem CID from SMILES string."""
    if not PUBCHEM_ENABLED:
        return None
    
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{requests.utils.quote(smiles)}/cids/JSON"
        resp = requests.get(url, timeout=5)
        resp.raise_for_status()
        data = resp.json()
        cids = data.get("IdentifierList", {}).get("CID", [])
        return cids[0] if cids else None
    except Exception as e:
        log.warning("Failed to get CID from SMILES %r: %s", smiles, e)
        return None


def _fetch_all_synonyms(cid: int) -> Optional[list[str]]:
    """Fetch all synonyms for a given PubChem CID."""
    if not PUBCHEM_ENABLED:
        return None
    
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/JSON"
        resp = requests.get(url, timeout=5)
        resp.raise_for_status()
        data = resp.json()
        synonyms = data.get("InformationList", {}).get("Information", [{}])[0].get("Synonym", [])
        return synonyms if synonyms else None
    except Exception as e:
        log.warning("Failed to fetch synonyms for CID %d: %s", cid, e)
        return None


def _filter_synonyms(synonyms: list[str], trivial_name: Optional[str] = None, max_synonyms: int = 10) -> Optional[list[str]]:
    """Filter synonyms to exclude duplicates and trivial name."""
    if not synonyms:
        return None
    
    # Track unique names (case-insensitive)
    seen_lower = set()
    if trivial_name:
        seen_lower.add(trivial_name.lower())
    
    unique_synonyms = []
    for s in synonyms:
        if len(s) >= 100:  # Skip very long names
            continue
        s_lower = s.lower()
        if s_lower not in seen_lower:
            seen_lower.add(s_lower)
            unique_synonyms.append(s)
            if len(unique_synonyms) >= max_synonyms:
                break
    
    return unique_synonyms if unique_synonyms else None


def _fetch_synonyms(cid: int, trivial_name: Optional[str] = None, max_synonyms: int = 10) -> Optional[list[str]]:
    """Fetch synonyms for a given PubChem CID, excluding duplicates and trivial name."""
    all_synonyms = _fetch_all_synonyms(cid)
    if not all_synonyms:
        return None
    return _filter_synonyms(all_synonyms, trivial_name=trivial_name, max_synonyms=max_synonyms)


def _from_smiles(query: str) -> Optional[ResolvedChem]:
    m = Chem.MolFromSmiles(query)
    if m is None:
        return None
    can = Chem.MolToSmiles(m, canonical=True)
    formula = rdMolDescriptors.CalcMolFormula(m)
    
    # Try to enrich with PubChem data
    trivial_name = None
    other_names = None
    name = None
    
    cid = _get_cid_from_smiles(can)
    if cid:
        all_synonyms = _fetch_all_synonyms(cid)
        if all_synonyms:
            # Use first synonym as trivial name if available
            trivial_name = all_synonyms[0] if all_synonyms else None
            # Get filtered other_names
            other_names = _filter_synonyms(all_synonyms, trivial_name=trivial_name, max_synonyms=10)
            # Use IUPAC name if available, otherwise first synonym
            name = trivial_name
    
    return ResolvedChem(
        smiles=can,
        mol=m,
        name=name,
        formula=formula,
        source="smiles",
        query=query,
        trivial_name=trivial_name,
        other_names=other_names,
    )


def _from_inchi(query: str) -> Optional[ResolvedChem]:
    if not query.startswith("InChI="):
        return None
    m = Chem.MolFromInchi(query)
    if m is None:
        return None
    can = Chem.MolToSmiles(m, canonical=True)
    formula = rdMolDescriptors.CalcMolFormula(m)
    
    # Try to enrich with PubChem data
    trivial_name = None
    other_names = None
    name = None
    
    cid = _get_cid_from_smiles(can)
    if cid:
        all_synonyms = _fetch_all_synonyms(cid)
        if all_synonyms:
            trivial_name = all_synonyms[0] if all_synonyms else None
            other_names = _filter_synonyms(all_synonyms, trivial_name=trivial_name, max_synonyms=10)
            name = trivial_name
    
    return ResolvedChem(
        smiles=can,
        mol=m,
        name=name,
        formula=formula,
        source="inchi",
        query=query,
        trivial_name=trivial_name,
        other_names=other_names,
    )


def _from_pubchem_name_or_cas(query: str) -> Optional[ResolvedChem]:
    if not PUBCHEM_ENABLED:
        return None

    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
        f"{requests.utils.quote(query)}/property/CanonicalSMILES,"
        "MolecularFormula,IUPACName,IsomericSMILES,SMILES,ConnectivitySMILES,Title/JSON"
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
        trivial_name = props.get("Title")

        m = Chem.MolFromSmiles(smiles)
        if m is None:
            log.warning("RDKit failed to parse SMILES from PubChem for %r: %r", query, smiles)
            return None

        can = Chem.MolToSmiles(m, canonical=True)
        if not formula:
            formula = rdMolDescriptors.CalcMolFormula(m)

        # Fetch synonyms (other names) for the compound
        cid = _get_cid_from_smiles(can)
        all_synonyms = None
        if cid:
            all_synonyms = _fetch_all_synonyms(cid)
        
        # If user's query matches a synonym (case-insensitive) and it's not a long IUPAC name,
        # use the properly cased synonym as the trivial name
        if all_synonyms and len(query) < 50 and not query.startswith("InChI"):
            query_lower = query.lower()
            for syn in all_synonyms:
                if syn.lower() == query_lower:
                    trivial_name = syn
                    break
        
        # Get filtered other_names (excluding trivial name and duplicates)
        other_names = None
        if all_synonyms:
            other_names = _filter_synonyms(all_synonyms, trivial_name=trivial_name, max_synonyms=10)

        return ResolvedChem(
            smiles=can,
            mol=m,
            name=name,
            formula=formula,
            source="cas" if CAS_RE.match(query) else "name",
            query=query,
            trivial_name=trivial_name,
            other_names=other_names,
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


