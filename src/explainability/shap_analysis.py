"""
SHAP Analysis - Fixed for SHAP 0.42+ API
Production-ready version
"""
import shap
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import io
from pathlib import Path

class SHAPAnalyzer:
    '''Fixed SHAP analyzer compatible with latest API'''
    
    def __init__(self):
        self.explainer = None
        self.shap_values = None
        self.base_value = None
        self.data = None
        self.feature_names = None
    
    def compute_shap_values(self, model, X, feature_names=None):
        '''
        Compute SHAP values with fixed API
        
        Parameters:
        - model: Trained XGBoost model
        - X: DataFrame or array with features
        - feature_names: Optional (auto-detected from DataFrame)
        '''
        try:
            # Convert to DataFrame if needed
            if not isinstance(X, pd.DataFrame):
                if feature_names is None:
                    feature_names = [f'feature_{i}' for i in range(X.shape[1])]
                X = pd.DataFrame(X, columns=feature_names)
            
            self.feature_names = X.columns.tolist()
            
            # Sample if too large (SHAP performance optimization)
            max_samples = 1000
            if len(X) > max_samples:
                X_sample = X.sample(n=max_samples, random_state=42)
            else:
                X_sample = X
            
            self.data = X_sample
            
            # Clean data
            if X_sample.isnull().any().any():
                X_sample = X_sample.fillna(0)
            X_sample = X_sample.replace([np.inf, -np.inf], 0)
            
            # Create explainer
            self.explainer = shap.TreeExplainer(model)
            
            # FIXED: Use new API - pass DataFrame directly, no feature_names parameter
            explanation = self.explainer(X_sample)
            
            # Extract values
            self.shap_values = explanation.values
            self.base_value = explanation.base_values
            
            return True
        
        except Exception as e:
            raise Exception(f'SHAP computation failed: {str(e)}')
    
    def generate_global_summary(self):
        '''Generate global SHAP summary'''
        if self.shap_values is None:
            raise ValueError('SHAP values not computed yet')
        
        try:
            plt.figure(figsize=(10, 6))
            
            shap.summary_plot(
                self.shap_values,
                self.data,
                show=False,
                max_display=10
            )
            
            buf = io.BytesIO()
            plt.savefig(buf, format='png', bbox_inches='tight', dpi=150)
            buf.seek(0)
            plt.close()
            
            # Calculate feature importance
            if len(self.shap_values.shape) == 2:
                importance = np.abs(self.shap_values).mean(axis=0)
            else:
                importance = np.abs(self.shap_values).mean(axis=(0, 2))
            
            top_features = pd.DataFrame({
                'feature': self.feature_names,
                'importance': importance
            }).sort_values('importance', ascending=False)
            
            return {
                'summary_plot': buf,
                'top_features': top_features.to_dict('records'),
                'computed': True,
                'n_samples_used': len(self.data)
            }
        
        except Exception as e:
            raise Exception(f'Summary generation failed: {str(e)}')
