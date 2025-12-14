"""
SHAP Analyzer V2 - Production Grade with Latest API
Compatible with SHAP 0.44+ (December 2025)
"""
import shap
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import io

class SHAPAnalyzerV2:
    def __init__(self):
        self.explainer = None
        self.shap_values = None
        self.base_value = None
        self.data = None
        self.feature_names = None
    
    def compute_shap_values(self, model, X_data, max_samples=1000):
        try:
            if not isinstance(X_data, pd.DataFrame):
                raise ValueError('X_data must be a pandas DataFrame')
            
            self.feature_names = X_data.columns.tolist()
            
            if len(X_data) > max_samples:
                X_sample = X_data.sample(n=max_samples, random_state=42)
            else:
                X_sample = X_data
            
            self.data = X_sample
            
            if X_sample.isnull().any().any():
                X_sample = X_sample.fillna(0)
            
            self.explainer = shap.TreeExplainer(model)
            explanation = self.explainer(X_sample)
            
            self.shap_values = explanation.values
            self.base_value = explanation.base_values
            
            return True
        except Exception as e:
            raise Exception(f'SHAP computation failed: {str(e)}')
    
    def generate_summary_plot(self):
        try:
            plt.figure(figsize=(10, 6))
            shap.summary_plot(self.shap_values, self.data, show=False, max_display=10)
            
            buf = io.BytesIO()
            plt.savefig(buf, format='png', bbox_inches='tight', dpi=150)
            buf.seek(0)
            plt.close()
            
            return buf
        except Exception as e:
            raise Exception(f'Plot generation failed: {str(e)}')
    
    def get_feature_importance(self):
        try:
            if len(self.shap_values.shape) == 2:
                importance = np.abs(self.shap_values).mean(axis=0)
            else:
                importance = np.abs(self.shap_values).mean(axis=(0, 2))
            
            importance_df = pd.DataFrame({
                'feature': self.feature_names,
                'importance': importance
            }).sort_values('importance', ascending=False)
            
            return importance_df
        except Exception as e:
            raise Exception(f'Importance calculation failed: {str(e)}')
