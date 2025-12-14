"""
Fixed SHAP Analyzer with Feature Pipeline Integration
"""
import shap
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import io
from pathlib import Path

class SHAPAnalyzer:
    def __init__(self):
        self.explainer = None
        self.shap_values = None
        self.feature_names = None
        self.X_data = None
    
    def compute_shap_values(self, model, signals_df, feature_pipeline):
        '''
        Compute SHAP values using the SAME feature pipeline as training
        
        Parameters:
        - model: Trained SISA model
        - signals_df: Raw signals DataFrame
        - feature_pipeline: The FeaturePipeline instance used during training
        '''
        try:
            # CRITICAL: Use feature pipeline to transform data exactly like training
            X = feature_pipeline.transform(signals_df)
            
            # Store for plotting
            self.X_data = X.copy()
            self.feature_names = feature_pipeline.feature_names
            
            # Sample if too large
            max_samples = 1000
            if len(X) > max_samples:
                X_sample = X.sample(n=max_samples, random_state=42)
                self.X_data = X_sample
            
            # Validate model can predict
            try:
                test_pred = model.predict(X_sample.iloc[:5] if len(X) > max_samples else X.iloc[:5])
            except Exception as e:
                raise ValueError(f'Model prediction test failed: {str(e)}')
            
            # Create TreeExplainer
            self.explainer = shap.TreeExplainer(model)
            
            # Compute SHAP values
            self.shap_values = self.explainer.shap_values(self.X_data)
            
            return True
        
        except Exception as e:
            raise Exception(f'SHAP computation failed: {str(e)}')
    
    def generate_global_summary(self):
        '''Generate global SHAP summary'''
        if self.shap_values is None:
            raise ValueError('SHAP values not computed yet')
        
        try:
            # Create summary plot
            plt.figure(figsize=(10, 6))
            shap.summary_plot(
                self.shap_values,
                self.X_data,
                feature_names=self.feature_names,
                show=False,
                max_display=10
            )
            
            # Save to buffer
            buf = io.BytesIO()
            plt.savefig(buf, format='png', bbox_inches='tight', dpi=150)
            buf.seek(0)
            plt.close()
            
            # Calculate feature importance
            if isinstance(self.shap_values, list):
                feature_importance = np.abs(self.shap_values[0]).mean(axis=0)
            else:
                feature_importance = np.abs(self.shap_values).mean(axis=0)
            
            top_features = pd.DataFrame({
                'feature': self.feature_names,
                'importance': feature_importance
            }).sort_values('importance', ascending=False)
            
            return {
                'summary_plot': buf,
                'top_features': top_features.to_dict('records'),
                'computed': True,
                'n_samples_used': len(self.X_data)
            }
        
        except Exception as e:
            raise Exception(f'Summary generation failed: {str(e)}')
