#!/usr/bin/env python3
"""Fill pKa values for datasets in ``new_datasets`` using local-first prediction.

This script mirrors ``fill_pka_temp_datasets.py`` structure but prefers
on-device inference and only falls back to PubChem when local methods cannot
resolve a value.

Backends:
- ``auto`` (default): try ``pkasolver`` first, then ``dimorphite_dl``, then
    PubChem.
- ``pkasolver``: require pkasolver and use it directly.
- ``dimorphite``: estimate a representative pKa by detecting protonation-state
    transition pH using ``dimorphite_dl``.

Output schema is normalized to: ``smiles,pka,target``.
"""

from __future__ import annotations

import argparse
import json
import re
import time
from collections import Counter
from pathlib import Path
from typing import Any
from urllib import error, parse, request

import pandas as pd
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize


DEFAULT_INPUT_DIR = Path("new_datasets")
DEFAULT_CACHE = Path("new_datasets/pka_local_cache.json")
DEFAULT_BACKEND = "auto"
USER_AGENT = "ToxiPred-pKa-Filler/1.0 (+https://pubchem.ncbi.nlm.nih.gov/)"
DEFAULT_API_ERROR_TTL_HOURS = 24.0


CacheEntry = dict[str, Any]

_PKASOLVER_QUERY_MODEL: Any | None = None


def discover_default_inputs(output_suffix: str) -> list[Path]:
    candidates = sorted(DEFAULT_INPUT_DIR.glob("*.csv"))
    return [path for path in candidates if not path.stem.endswith(output_suffix)]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Fill pka column for CSV files under new_datasets using local methods "
            "with PubChem fallback in auto mode."
        )
    )
    parser.add_argument(
        "--inputs",
        nargs="+",
        type=Path,
        default=None,
        help="Input CSV files. Defaults to all CSVs under new_datasets/.",
    )
    parser.add_argument(
        "--suffix",
        default="_with_pka",
        help="Suffix appended before .csv for output files.",
    )
    parser.add_argument(
        "--cache-file",
        type=Path,
        default=DEFAULT_CACHE,
        help="JSON cache file for resolved pKa values.",
    )
    parser.add_argument(
        "--backend",
        choices=["auto", "pkasolver", "dimorphite"],
        default=DEFAULT_BACKEND,
        help="Local pKa backend selection.",
    )
    parser.add_argument(
        "--max-rows",
        type=int,
        default=None,
        help="Optional max rows to process per input file for quick tests.",
    )
    parser.add_argument(
        "--ph-min",
        type=float,
        default=0.0,
        help="Minimum pH for dimorphite transition scan.",
    )
    parser.add_argument(
        "--ph-max",
        type=float,
        default=14.0,
        help="Maximum pH for dimorphite transition scan.",
    )
    parser.add_argument(
        "--ph-step",
        type=float,
        default=0.5,
        help="pH step size for dimorphite transition scan.",
    )
    parser.add_argument(
        "--max-variants",
        type=int,
        default=64,
        help="Maximum protonation variants per pH point for dimorphite.",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=15.0,
        help="HTTP timeout in seconds for PubChem requests.",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=3,
        help="Number of HTTP retries per PubChem request.",
    )
    parser.add_argument(
        "--sleep",
        type=float,
        default=0.2,
        help="Sleep time between PubChem API calls in seconds.",
    )
    parser.add_argument(
        "--retry-api-errors",
        action="store_true",
        default=True,
        help=(
            "Retry compounds cached as PubChem API errors from previous runs. "
            "Enabled by default."
        ),
    )
    parser.add_argument(
        "--no-retry-api-errors",
        dest="retry_api_errors",
        action="store_false",
        help="Do not retry compounds previously cached as PubChem API errors.",
    )
    parser.add_argument(
        "--api-error-ttl-hours",
        type=float,
        default=DEFAULT_API_ERROR_TTL_HOURS,
        help=(
            "Retry cached api_error entries only when older than this many hours. "
            "Use 0 to retry all api_error entries every run."
        ),
    )
    return parser.parse_args()


def load_cache(cache_path: Path) -> dict[str, CacheEntry]:
    if not cache_path.exists():
        return {}

    try:
        raw = json.loads(cache_path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return {}

    cache: dict[str, CacheEntry] = {}
    if not isinstance(raw, dict):
        return cache

    for key, value in raw.items():
        if not isinstance(key, str):
            continue

        # Backward-compatible numeric cache values.
        if isinstance(value, (int, float)):
            cache[key] = {"status": "resolved", "value": float(value)}
            continue

        if not isinstance(value, dict):
            continue

        status = value.get("status")
        if status == "resolved" and isinstance(value.get("value"), (int, float)):
            cache[key] = {"status": "resolved", "value": float(value["value"])}
        elif status in {"not_found", "invalid_smiles", "predict_error"}:
            cache[key] = {
                "status": status,
                "error": str(value.get("error", "")),
            }
        elif status == "api_error":
            err_at = value.get("error_at")
            cache[key] = {
                "status": "api_error",
                "error": str(value.get("error", "")),
                "error_at": float(err_at) if isinstance(err_at, (int, float)) else None,
            }

    return cache


def should_retry_api_error(
    entry: CacheEntry,
    retry_api_errors: bool,
    api_error_ttl_hours: float,
    now_epoch: float,
) -> bool:
    if not retry_api_errors:
        return False

    if api_error_ttl_hours <= 0:
        return True

    error_at = entry.get("error_at")
    if not isinstance(error_at, (int, float)):
        return True

    age_seconds = now_epoch - float(error_at)
    ttl_seconds = api_error_ttl_hours * 3600.0
    return age_seconds >= ttl_seconds


def save_cache(cache_path: Path, cache: dict[str, CacheEntry]) -> None:
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(
        json.dumps(cache, indent=2, sort_keys=True, ensure_ascii=True),
        encoding="utf-8",
    )


def canonicalize_smiles(smiles: str) -> str | None:
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
        if mol is None:
            return None

        chooser = rdMolStandardize.LargestFragmentChooser()
        mol = chooser.choose(mol)
        Chem.SanitizeMol(mol)
        return Chem.MolToSmiles(mol, canonical=True)
    except Exception:
        return None


def fetch_json(url: str, timeout: float, retries: int) -> Any:
    last_error: Exception | None = None
    for attempt in range(retries):
        req = request.Request(url, headers={"User-Agent": USER_AGENT})
        try:
            with request.urlopen(req, timeout=timeout) as response:
                return json.loads(response.read().decode("utf-8"))
        except (
            error.URLError,
            error.HTTPError,
            TimeoutError,
            json.JSONDecodeError,
        ) as exc:
            last_error = exc
            if attempt < retries - 1:
                time.sleep(0.6 * (attempt + 1))
    if last_error is not None:
        raise last_error
    raise RuntimeError("Unexpected network error without exception")


def _extract_numbers_from_text(text: str) -> list[float]:
    matches = re.findall(r"[+-]?\d+(?:\.\d+)?", text)
    numbers: list[float] = []
    for token in matches:
        try:
            numbers.append(float(token))
        except ValueError:
            continue
    return numbers


def _collect_candidate_numbers(node: Any) -> list[float]:
    numbers: list[float] = []

    if isinstance(node, dict):
        if "Number" in node and isinstance(node["Number"], list):
            for item in node["Number"]:
                if isinstance(item, (int, float)):
                    numbers.append(float(item))

        if "String" in node and isinstance(node["String"], str):
            text = node["String"]
            if "pka" in text.lower():
                numbers.extend(_extract_numbers_from_text(text))

        if "StringWithMarkup" in node and isinstance(node["StringWithMarkup"], list):
            for item in node["StringWithMarkup"]:
                if not isinstance(item, dict):
                    continue
                text = item.get("String")
                if isinstance(text, str) and "pka" in text.lower():
                    numbers.extend(_extract_numbers_from_text(text))

        for value in node.values():
            numbers.extend(_collect_candidate_numbers(value))

    elif isinstance(node, list):
        for item in node:
            numbers.extend(_collect_candidate_numbers(item))

    return numbers


def _extract_pka_from_pug_view(data: Any) -> float | None:
    candidates = _collect_candidate_numbers(data)
    if not candidates:
        return None
    return candidates[0]


def _resolve_cid_for_smiles(smiles: str, timeout: float, retries: int) -> int | None:
    quoted = parse.quote(smiles, safe="")
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/"
        f"{quoted}/cids/JSON"
    )
    payload = fetch_json(url, timeout=timeout, retries=retries)
    cids = (
        payload.get("IdentifierList", {}).get("CID", [])
        if isinstance(payload, dict)
        else []
    )
    if cids and isinstance(cids[0], int):
        return cids[0]
    return None


def resolve_pka_from_pubchem(smiles: str, timeout: float, retries: int) -> float | None:
    cid = _resolve_cid_for_smiles(smiles, timeout=timeout, retries=retries)
    if cid is None:
        return None

    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/"
        f"{cid}/JSON?heading={parse.quote('Dissociation Constants', safe='')}"
    )
    payload = fetch_json(url, timeout=timeout, retries=retries)
    return _extract_pka_from_pug_view(payload)


def _extract_first_float(value: Any) -> float | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        return float(value)

    if isinstance(value, dict):
        preferred_keys = [
            "pka",
            "pKa",
            "value",
            "predicted_pka",
            "macro_pka",
            "micro_pka",
        ]
        for key in preferred_keys:
            if key in value:
                out = _extract_first_float(value[key])
                if out is not None:
                    return out
        for nested in value.values():
            out = _extract_first_float(nested)
            if out is not None:
                return out
        return None

    if isinstance(value, (list, tuple, set)):
        for item in value:
            out = _extract_first_float(item)
            if out is not None:
                return out
        return None

    for attr in ("pka", "value"):
        if hasattr(value, attr):
            out = _extract_first_float(getattr(value, attr))
            if out is not None:
                return out

    return None


def _call_pkasolver(smiles: str) -> float | None:
    global _PKASOLVER_QUERY_MODEL

    try:
        from pkasolver.query import QueryModel, calculate_microstate_pka_values
    except Exception:
        return None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    try:
        if _PKASOLVER_QUERY_MODEL is None:
            _PKASOLVER_QUERY_MODEL = QueryModel()

        states = calculate_microstate_pka_values(
            mol,
            query_model=_PKASOLVER_QUERY_MODEL,
        )
    except Exception:
        return None

    if not states:
        return None

    # Choose a single representative pKa: closest to physiological pH.
    candidate_values = [float(state.pka) for state in states if state.pka is not None]
    if not candidate_values:
        return None

    return min(candidate_values, key=lambda value: (abs(value - 7.0), value))


def _dominant_formal_charge(smiles_variants: list[str]) -> int | None:
    charges: list[int] = []
    for variant in smiles_variants:
        mol = Chem.MolFromSmiles(variant)
        if mol is None:
            continue
        charges.append(int(Chem.GetFormalCharge(mol)))

    if not charges:
        return None

    counts = Counter(charges)
    best_count = max(counts.values())
    tied = [charge for charge, count in counts.items() if count == best_count]
    # Deterministic tie-breaker: prefer charge closest to neutral.
    return sorted(tied, key=lambda x: (abs(x), x))[0]


def _resolve_with_dimorphite(
    smiles: str,
    ph_min: float,
    ph_max: float,
    ph_step: float,
    max_variants: int,
) -> float | None:
    from dimorphite_dl import protonate_smiles

    if ph_step <= 0 or ph_max < ph_min:
        return None

    ph_points: list[float] = []
    charge_curve: list[int | None] = []

    p = ph_min
    while p <= ph_max + 1e-9:
        try:
            variants = protonate_smiles(
                smiles,
                ph_min=float(p),
                ph_max=float(p),
                precision=0.0,
                max_variants=max_variants,
                validate_output=True,
            )
        except Exception:
            return None

        dominant_charge = _dominant_formal_charge(list(variants))
        ph_points.append(float(p))
        charge_curve.append(dominant_charge)
        p += ph_step

    transitions: list[tuple[float, int, int]] = []
    for i in range(len(ph_points) - 1):
        c0 = charge_curve[i]
        c1 = charge_curve[i + 1]
        if c0 is None or c1 is None:
            continue
        if c0 != c1:
            midpoint = (ph_points[i] + ph_points[i + 1]) / 2.0
            transitions.append((midpoint, c0, c1))

    if not transitions:
        return None

    # Prefer transitions crossing/straddling neutral charge, else earliest one.
    neutral_cross = [
        item
        for item in transitions
        if (item[1] <= 0 <= item[2]) or (item[2] <= 0 <= item[1])
    ]
    chosen = neutral_cross[0] if neutral_cross else transitions[0]
    return float(chosen[0])


def resolve_pka_local(
    smiles: str,
    backend: str,
    ph_min: float,
    ph_max: float,
    ph_step: float,
    max_variants: int,
) -> tuple[float | None, str]:
    if backend in {"auto", "pkasolver"}:
        pka_value = _call_pkasolver(smiles)
        if pka_value is not None:
            return pka_value, "pkasolver"
        if backend == "pkasolver":
            return None, "pkasolver"

    if backend in {"auto", "dimorphite"}:
        try:
            pka_value = _resolve_with_dimorphite(
                smiles=smiles,
                ph_min=ph_min,
                ph_max=ph_max,
                ph_step=ph_step,
                max_variants=max_variants,
            )
        except Exception:
            return None, "dimorphite"
        return pka_value, "dimorphite"

    return None, backend


def read_dataset(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, dtype=str, keep_default_na=False, na_values=[])
    rename_map = {col: col.strip() for col in df.columns}
    df = df.rename(columns=rename_map)

    lowered = {col.lower(): col for col in df.columns}
    if "smiles" not in lowered:
        raise ValueError(f"Missing smiles column in {path}")

    if "pka" not in lowered:
        df["pka"] = ""

    if "target" not in lowered:
        df["target"] = ""

    smiles_col = lowered["smiles"]
    pka_col = lowered.get("pka", "pka")
    target_col = lowered.get("target", "target")

    out = pd.DataFrame(
        {
            "smiles": df[smiles_col],
            "pka": df[pka_col] if pka_col in df.columns else "",
            "target": df[target_col] if target_col in df.columns else "",
        }
    )
    return out


def output_path_for(input_path: Path, suffix: str) -> Path:
    return input_path.with_name(f"{input_path.stem}{suffix}{input_path.suffix}")


def process_file(
    input_path: Path,
    output_path: Path,
    cache: dict[str, CacheEntry],
    backend: str,
    max_rows: int | None,
    ph_min: float,
    ph_max: float,
    ph_step: float,
    max_variants: int,
    timeout: float,
    retries: int,
    sleep_seconds: float,
    retry_api_errors: bool,
    api_error_ttl_hours: float,
) -> dict[str, int]:
    df = read_dataset(input_path)

    stats = {
        "total_rows": len(df),
        "processed_rows": 0,
        "valid_smiles": 0,
        "resolved_pka": 0,
        "unresolved_pka": 0,
        "cache_hits": 0,
        "cache_api_error_retries": 0,
        "cache_api_error_ttl_skips": 0,
        "api_errors": 0,
        "predict_errors": 0,
        "malformed_rows": 0,
        "pkasolver_resolved": 0,
        "dimorphite_resolved": 0,
        "pubchem_resolved": 0,
    }

    limit = len(df) if max_rows is None else min(max_rows, len(df))
    now_epoch = time.time()

    for idx in range(limit):
        stats["processed_rows"] += 1
        smiles_raw = str(df.at[idx, "smiles"]).strip()

        if not smiles_raw:
            stats["malformed_rows"] += 1
            df.at[idx, "pka"] = ""
            continue

        canonical = canonicalize_smiles(smiles_raw)
        if canonical is None:
            stats["malformed_rows"] += 1
            df.at[idx, "pka"] = ""
            cache[smiles_raw] = {"status": "invalid_smiles", "error": "rdkit-parse"}
            continue

        stats["valid_smiles"] += 1

        if canonical in cache:
            stats["cache_hits"] += 1
            entry = cache[canonical]
            status = entry.get("status")

            if status == "resolved":
                cached_val = entry.get("value")
                if isinstance(cached_val, (int, float)):
                    df.at[idx, "pka"] = f"{float(cached_val):.4f}"
                    stats["resolved_pka"] += 1
                    continue

            if status == "api_error":
                should_retry = should_retry_api_error(
                    entry=entry,
                    retry_api_errors=retry_api_errors,
                    api_error_ttl_hours=api_error_ttl_hours,
                    now_epoch=now_epoch,
                )
                if should_retry:
                    stats["cache_api_error_retries"] += 1
                else:
                    df.at[idx, "pka"] = ""
                    stats["unresolved_pka"] += 1
                    stats["cache_api_error_ttl_skips"] += 1
                    continue

            if status in {"not_found", "invalid_smiles", "predict_error"}:
                df.at[idx, "pka"] = ""
                stats["unresolved_pka"] += 1
                continue

        pka_value, used_backend = resolve_pka_local(
            smiles=canonical,
            backend=backend,
            ph_min=ph_min,
            ph_max=ph_max,
            ph_step=ph_step,
            max_variants=max_variants,
        )

        if pka_value is None and backend == "auto":
            try:
                pka_value = resolve_pka_from_pubchem(
                    canonical,
                    timeout=timeout,
                    retries=retries,
                )
                used_backend = "pubchem" if pka_value is not None else used_backend
            except Exception as exc:
                cache[canonical] = {
                    "status": "api_error",
                    "error": str(exc),
                    "error_at": time.time(),
                }
                df.at[idx, "pka"] = ""
                stats["unresolved_pka"] += 1
                stats["api_errors"] += 1
                stats["predict_errors"] += 1
                if sleep_seconds > 0:
                    time.sleep(sleep_seconds)
                continue

        if pka_value is None:
            cache[canonical] = {
                "status": "not_found",
                "error": f"unresolved-{used_backend}",
            }
            df.at[idx, "pka"] = ""
            stats["unresolved_pka"] += 1
            stats["predict_errors"] += 1
            continue

        cache[canonical] = {"status": "resolved", "value": float(pka_value)}
        df.at[idx, "pka"] = f"{float(pka_value):.4f}"
        stats["resolved_pka"] += 1
        if used_backend == "pkasolver":
            stats["pkasolver_resolved"] += 1
        elif used_backend == "dimorphite":
            stats["dimorphite_resolved"] += 1
        elif used_backend == "pubchem":
            stats["pubchem_resolved"] += 1

        if used_backend == "pubchem" and sleep_seconds > 0:
            time.sleep(sleep_seconds)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)
    return stats


def print_stats(input_path: Path, output_path: Path, stats: dict[str, int]) -> None:
    print(f"\nProcessed: {input_path}")
    print(f"Output:    {output_path}")
    print(f"Rows total:           {stats['total_rows']}")
    print(f"Rows processed:       {stats['processed_rows']}")
    print(f"Valid SMILES:         {stats['valid_smiles']}")
    print(f"Resolved pKa:         {stats['resolved_pka']}")
    print(f"  via pkasolver:      {stats['pkasolver_resolved']}")
    print(f"  via dimorphite:     {stats['dimorphite_resolved']}")
    print(f"  via pubchem:        {stats['pubchem_resolved']}")
    print(f"Unresolved pKa:       {stats['unresolved_pka']}")
    print(f"Malformed rows:       {stats['malformed_rows']}")
    print(f"Cache hits:           {stats['cache_hits']}")
    print(f"Retried API errors:   {stats['cache_api_error_retries']}")
    print(f"TTL-skipped errors:   {stats['cache_api_error_ttl_skips']}")
    print(f"API errors:           {stats['api_errors']}")
    print(f"Prediction errors:    {stats['predict_errors']}")


def main() -> int:
    args = parse_args()

    inputs = args.inputs if args.inputs else discover_default_inputs(args.suffix)
    if not inputs:
        print("No input CSV files found.")
        return 1

    cache = load_cache(args.cache_file)
    print(f"Loaded cache entries: {len(cache)}")
    print(f"Backend mode: {args.backend}")

    if args.backend == "auto":
        print(
            "Note: auto backend prefers pkasolver; falls back to dimorphite, then PubChem."
        )

    for input_path in inputs:
        if not input_path.exists():
            print(f"Skipping missing input: {input_path}")
            continue

        out_path = output_path_for(input_path, args.suffix)
        stats = process_file(
            input_path=input_path,
            output_path=out_path,
            cache=cache,
            backend=args.backend,
            max_rows=args.max_rows,
            ph_min=args.ph_min,
            ph_max=args.ph_max,
            ph_step=args.ph_step,
            max_variants=args.max_variants,
            timeout=args.timeout,
            retries=args.retries,
            sleep_seconds=args.sleep,
            retry_api_errors=args.retry_api_errors,
            api_error_ttl_hours=args.api_error_ttl_hours,
        )
        print_stats(input_path, out_path, stats)

    save_cache(args.cache_file, cache)
    print(f"\nSaved cache entries: {len(cache)}")
    print(f"Cache file: {args.cache_file}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
