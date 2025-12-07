import pandas as pd
import numpy as np
from pathlib import Path
import mlflow
import os

def prepare_ml_features_from_stats_engine():
    """Convert stats_engine candidates → ML training features"""
    BASE_DIR = Path.cwd()
    
    # Create ml_data directory
    (BASE_DIR / 'ml_data').mkdir(exist_ok=True)
    
    # Load stats engine candidates
    candidates = pd.read_csv(BASE_DIR / 'sar_reports' / 'candidate_signals_2025Q1_stats_engine.csv')
    
    # ML Features from stats engine (regulatory gold standard)
    features = pd.DataFrame({
        'prr': candidates['PRR'],
        'chi2': candidates['CHISQ'],
        'cases': candidates['CASES'],
        'log_prr': np.log1p(candidates['PRR'].clip(lower=0)),
        'prr_rank': candidates['PRR'].rank(pct=True),
        'chi2_rank': candidates['CHISQ'].rank(pct=True),
        'cases_log': np.log1p(candidates['CASES']),
        'candidate_signal': 1,  # All are candidates by definition
        'meets_prr': 1,
        'meets_chi2': 1,
        'meets_min_cases': 1
    })
    
    # Save for XGBoost training
    output_path = BASE_DIR / 'ml_data' / 'stats_engine_features.csv'
    features.to_csv(output_path, index=False)
    print(f'✅ ML features ready: {features.shape} rows, {features.shape[1]} features')
    print(f'💾 Saved: {output_path}')
    print('\nFeature summary:')
    print(features.describe())
    
    return features

if __name__ == '__main__':
    features = prepare_ml_features_from_stats_engine()
