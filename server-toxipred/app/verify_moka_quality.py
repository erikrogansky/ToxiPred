#!/usr/bin/env python3
"""Verify the quality of moka pKa predictions against toxicity targets."""

import json
from pathlib import Path
import pandas as pd
from scipy.stats import spearmanr, pearsonr
import numpy as np


def load_cache(cache_path: Path) -> dict:
    """Load pKa cache to identify which predictions came from moka."""
    if not cache_path.exists():
        return {}
    
    try:
        return json.loads(cache_path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return {}


def analyze_dataset(
    predictions_file: Path,
    cache: dict,
    dataset_name: str
) -> dict:
    """Analyze quality of moka predictions for a single dataset."""
    
    if not predictions_file.exists():
        print(f"⚠️  {dataset_name}: File not found - {predictions_file}")
        return {}
    
    # Load predictions
    df = pd.read_csv(predictions_file, dtype={'pka': str, 'target': str})
    
    # Convert to numeric, coerce errors to NaN
    df['pka_num'] = pd.to_numeric(df['pka'], errors='coerce')
    df['target_num'] = pd.to_numeric(df['target'], errors='coerce')
    
    # Identify prediction source for each SMILES
    df['prediction_source'] = 'unknown'
    for idx, row in df.iterrows():
        canonical_smiles = row['smiles']
        if canonical_smiles in cache:
            entry = cache[canonical_smiles]
            if entry.get('status') == 'resolved':
                # We can't directly tell the source from cache alone,
                # but we can infer: if predictions were made this run, check stats
                df.at[idx, 'prediction_source'] = 'local'
    
    # Filter valid rows
    valid = df[
        (df['pka_num'].notna()) & 
        (df['target_num'].notna()) & 
        (df['pka_num'] > 0) & 
        (df['pka_num'] < 15)
    ].copy()
    
    if len(valid) == 0:
        return {}
    
    # Calculate correlations
    spearman_corr, spearman_p = spearmanr(valid['pka_num'], valid['target_num'])
    pearson_corr, pearson_p = pearsonr(valid['pka_num'], valid['target_num'])
    
    # Summary statistics
    results = {
        'dataset': dataset_name,
        'total_rows': len(df),
        'valid_predictions': len(valid),
        'pka_mean': valid['pka_num'].mean(),
        'pka_std': valid['pka_num'].std(),
        'pka_min': valid['pka_num'].min(),
        'pka_max': valid['pka_num'].max(),
        'toxicity_positive': (valid['target_num'] == 1).sum(),
        'toxicity_negative': (valid['target_num'] == 0).sum(),
        'spearman_r': spearman_corr,
        'spearman_p': spearman_p,
        'pearson_r': pearson_corr,
        'pearson_p': pearson_p,
    }
    
    # Correlation by toxicity class
    toxic = valid[valid['target_num'] == 1]['pka_num']
    non_toxic = valid[valid['target_num'] == 0]['pka_num']
    
    results['toxic_pka_mean'] = toxic.mean()
    results['toxic_pka_std'] = toxic.std()
    results['non_toxic_pka_mean'] = non_toxic.mean()
    results['non_toxic_pka_std'] = non_toxic.std()
    
    # Effect size (Cohen's d)
    if len(toxic) > 0 and len(non_toxic) > 0:
        pooled_std = np.sqrt(
            (toxic.std()**2 + non_toxic.std()**2) / 2
        )
        cohens_d = (toxic.mean() - non_toxic.mean()) / pooled_std if pooled_std > 0 else 0
        results['cohens_d'] = cohens_d
    
    return results


def print_results(all_results: list[dict]) -> None:
    """Print formatted results."""
    print("\n" + "="*80)
    print("MOKA PREDICTION QUALITY VERIFICATION")
    print("="*80)
    
    for results in all_results:
        if not results:
            continue
        
        dataset = results['dataset']
        print(f"\n📊 {dataset}")
        print("-" * 80)
        print(f"  Total rows:              {results['total_rows']}")
        print(f"  Valid predictions:       {results['valid_predictions']}")
        print(f"  Toxicity: {results['toxicity_positive']} positive, {results['toxicity_negative']} negative")
        print()
        print(f"  pKa Statistics:")
        print(f"    Mean ± SD:              {results['pka_mean']:.2f} ± {results['pka_std']:.2f}")
        print(f"    Range:                  {results['pka_min']:.2f} – {results['pka_max']:.2f}")
        print()
        print(f"  pKa by Toxicity Class:")
        print(f"    Toxic (Y=1):            {results['toxic_pka_mean']:.2f} ± {results['toxic_pka_std']:.2f}")
        print(f"    Non-toxic (Y=0):        {results['non_toxic_pka_mean']:.2f} ± {results['non_toxic_pka_std']:.2f}")
        print()
        print(f"  Correlation with Toxicity:")
        print(f"    Spearman r:             {results['spearman_r']:+.4f} (p={results['spearman_p']:.3e})")
        print(f"    Pearson r:              {results['pearson_r']:+.4f} (p={results['pearson_p']:.3e})")
        print(f"    Cohen's d (effect):     {results.get('cohens_d', 'N/A')}")
        print()
        
        # Interpretation
        r = results['spearman_r']
        p = results['spearman_p']
        
        if p < 0.05:
            sig = "✅ SIGNIFICANT"
            if abs(r) < 0.1:
                strength = "negligible"
            elif abs(r) < 0.3:
                strength = "weak"
            elif abs(r) < 0.5:
                strength = "moderate"
            elif abs(r) < 0.7:
                strength = "strong"
            else:
                strength = "very strong"
            print(f"  Interpretation: {sig} {strength} correlation")
        else:
            print(f"  Interpretation: ⚠️  NOT SIGNIFICANT (p={p:.3f})")
        print()


def main():
    cache_path = Path("new_datasets/pka_local_cache.json")
    cache = load_cache(cache_path)
    
    print(f"Loaded {len(cache)} pKa predictions from cache")
    
    # Analyze each dataset
    datasets = [
        (Path("new_datasets/in_vitro_corrosion_with_pka.csv"), "in_vitro_corrosion"),
        (Path("new_datasets/in_vitro_irritation_with_pka.csv"), "in_vitro_irritation"),
        (Path("new_datasets/in_vivo_corrosion_with_pka.csv"), "in_vivo_corrosion"),
    ]
    
    results = []
    for pred_file, name in datasets:
        result = analyze_dataset(pred_file, cache, name)
        results.append(result)
    
    print_results(results)
    
    # Summary
    print("="*80)
    print("SUMMARY & INTERPRETATION")
    print("="*80)
    print("""
Quality Assessment:
  • |r| > 0.3 with p < 0.05 → Moka captures relevant chemistry ✅
  • |r| < 0.2 with p < 0.05 → Moka has weak predictive signal ⚠️
  • p > 0.05 → No significant relationship ❌

Expected Range for pKa:
  • Most drug-like molecules: pKa 2.0–10.0
  • Toxicity-relevant pKa: Often 3–8 (carboxylic acids, amines, phenols)

Note: Even weak correlations are useful for ML models that learn non-linear
relationships. Moka's role is to provide a baseline feature, not predict toxicity
directly.
    """)


if __name__ == "__main__":
    main()
