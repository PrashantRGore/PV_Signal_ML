"""
SHAP Analyzer V2 - Production Grade with Latest API
Compatible with SHAP 0.44+ (December 2025)
"""
import shap
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import io
from pathlib import Path

class SHAPAnalyzerV2:
    '''
    Production-grade SHAP analyzer with:
    - Latest SHAP 0.44+ API
    - Automatic feature detection
    - Memory-efficient sampling
    - Multiple explanation types
    '''
    
    def __init__(self):
        self.explainer = None
        self.shap_values = None
        self.base_value = None
        self.data = None
        self.feature_names = None
    
    def compute_shap_values(self, model, X_data, max_samples=1000):
        '''
        Compute SHAP values using modern API
        
        Parameters:
        - model: Trained model (XGBoost, LightGBM, etc.)
        - X_data: DataFrame or array with features
        - max_samples: Maximum samples for efficiency
        '''
        try:
            # Ensure DataFrame
            if not isinstance(X_data, pd.DataFrame):
                raise ValueError('X_data must be a pandas DataFrame')
            
            # Store feature names from DataFrame columns
            self.feature_names = X_data.columns.tolist()
            
            # Sample if too large
            if len(X_data) > max_samples:
                X_sample = X_data.sample(n=max_samples, random_state=42)
            else:
                X_sample = X_data
            
            self.data = X_sample
            
            # Validate data
            if X_sample.isnull().any().any():
                X_sample = X_sample.fillna(0)
            
            # Create explainer - NO feature_names parameter needed!
            self.explainer = shap.TreeExplainer(model)
            
            # Compute SHAP values - feature names auto-detected from DataFrame
            explanation = self.explainer(X_sample)
            
            # Extract values and base
            self.shap_values = explanation.values
            self.base_value = explanation.base_values
            
            return True
        
        except Exception as e:
            raise Exception(f'SHAP computation failed: {str(e)}')
    
    def generate_summary_plot(self):
        '''Generate SHAP summary plot'''
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
            
            return buf
        
        except Exception as e:
            raise Exception(f'Plot generation failed: {str(e)}')
    
    def get_feature_importance(self):
        '''Calculate global feature importance'''
        try:
            # Handle both single and multi-output
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
    
    def generate_waterfall_plot(self, sample_idx=0):
        '''Generate waterfall plot for individual prediction'''
        try:
            plt.figure(figsize=(10, 6))
            
            shap.plots.waterfall(
                shap.Explanation(
                    values=self.shap_values[sample_idx],
                    base_values=self.base_value[sample_idx] if isinstance(self.base_value, np.ndarray) else self.base_value,
                    data=self.data.iloc[sample_idx],
                    feature_names=self.feature_names
                ),
                show=False
            )
            
            buf = io.BytesIO()
            plt.savefig(buf, format='png', bbox_inches='tight', dpi=150)
            buf.seek(0)
            plt.close()
            
            return buf
        
        except Exception as e:
            raise Exception(f'Waterfall plot failed: {str(e)}')
