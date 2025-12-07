import pandas as pd
import numpy as np
from pathlib import Path
import json
from datetime import datetime
import joblib

class PVDriftMonitor:
    def __init__(self, baseline_file='ml_data/stats_engine_features.csv'):
        self.baseline_file = Path(baseline_file)
        self.baseline_stats = self._load_baseline()
        self.drift_scores = {}
        
    def _load_baseline(self):
        df = pd.read_csv(self.baseline_file)
        return {
            'prr_mean': df['prr'].mean(),
            'prr_std': df['prr'].std(),
            'chi2_mean': df['chi2'].mean(),
            'chi2_std': df['chi2'].std(),
            'cases_mean': df['cases'].mean(),
            'feature_cols': df.columns.tolist()
        }
    
    def monitor_drift(self, new_data_file):
        """Detect statistical drift in new FAERS quarters"""
        new_df = pd.read_csv(new_data_file)
        drift_report = {}
        
        for col in ['prr', 'chi2', 'cases']:
            baseline_mean = self.baseline_stats[f'{col}_mean']
            baseline_std = self.baseline_stats[f'{col}_std']
            new_mean = new_df[col].mean()
            
            # Z-score drift detection
            z_score = abs((new_mean - baseline_mean) / baseline_std)
            drift_report[f'{col}_zscore'] = z_score
            drift_report[f'{col}_alert'] = z_score > 2.0  # 2σ threshold
        
        self.drift_scores[datetime.now().isoformat()] = drift_report
        self._save_drift_report(drift_report)
        return drift_report
    
    def _save_drift_report(self, report):
        Path('governance').mkdir(exist_ok=True)
        with open(f'governance/drift_report_{datetime.now().strftime("%Y%m%d")}.json', 'w') as f:
            json.dump(report, f, indent=2)
        
if __name__ == "__main__":
    monitor = PVDriftMonitor()
    print("✅ Baseline loaded:", len(monitor.baseline_stats['feature_cols']), "features")
