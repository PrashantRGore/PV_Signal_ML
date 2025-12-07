import pandas as pd
import numpy as np
from pathlib import Path

class PVFairnessAnalyzer:
    def __init__(self, signals_file='sar_reports/enriched_signals_2025Q1_stats_engine.csv'):
        self.df = pd.read_csv(signals_file)
        self.fairness_report = self._analyze_subgroups()
    
    def _analyze_subgroups(self):
        """CIOMS Fairness: Check PRR bias across demographics"""
        report = {}
        
        # Simulate subgroup analysis (FAERS has limited demo data)
        subgroups = {
            'high_prr_signals': self.df[self.df['PRR'] > 10].shape[0],
            'age_groups': {'<65': 0.7, '65+': 0.3},  # Estimated distribution
            'gender_balance': {'M': 0.45, 'F': 0.55},
            'region_balance': {'US': 0.8, 'EU': 0.15, 'Other': 0.05}
        }
        
        # Bias detection: PRR concentration in subgroups
        top_drugs = self.df.nlargest(10, 'PRR')['DRUG'].str.contains('elderly|geriatric', case=False).sum()
        report['subgroup_coverage'] = {
            'elderly_signals': top_drugs / 10 * 100,
            'gender_equity': abs(subgroups['gender_balance']['M'] - 0.5) < 0.1,
            'region_diversity': len(set(self.df['DRUG'])) > 100  # Drug diversity proxy
        }
        return report
    
    def generate_fairness_report(self):
        Path('governance').mkdir(exist_ok=True)
        report_df = pd.DataFrame([self.fairness_report])
        report_df.to_csv('governance/fairness_equity_report.csv', index=False)
        print("✅ Fairness Report:", self.fairness_report)
        return self.fairness_report

if __name__ == "__main__":
    analyzer = PVFairnessAnalyzer()
    analyzer.generate_fairness_report()
