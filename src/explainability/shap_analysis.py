"""
SHAP Analysis - FIXED VERSION
Bug #4: Corrected instantiation pattern (no model in __init__)
"""
import shap
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from pathlib import Path

class SHAPAnalyzer:
    """
    ✅ FIXED: Model is passed to compute_shap_values(), not __init__()
    """
    
    def __init__(self):  # ✅ No model parameter
        """Initialize SHAP analyzer."""
        self.explainer = None
        self.shap_values = None
        self.base_value = None
        self.data = None
        self.feature_names = None
    
    def compute_shap_values(self, model, X, feature_names=None):
        """
        Compute SHAP values with fixed API.
        
        ✅ FIXED: Model passed here, not in __init__
        
        Parameters:
        - model: Trained XGBoost model (SISA ensemble or single model)
        - X: DataFrame or array with features
        - feature_names: Optional (auto-detected from DataFrame)
        
        Returns:
        - bool: True if successful
        """
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
            
            # ✅ FIXED: Handle SISA ensemble (list of models) vs single model
            if isinstance(model, list):
                # SISA ensemble - use first shard for SHAP (all shards have same features)
                model = model[0]
            
            # Create explainer
            self.explainer = shap.TreeExplainer(model)
            
            # Compute SHAP values using new API
            explanation = self.explainer(X_sample)
            
            # Extract values
            self.shap_values = explanation.values
            self.base_value = explanation.base_values
            
            return True
        
        except Exception as e:
            raise Exception(f'SHAP computation failed: {str(e)}')
    
    def get_feature_importance(self):
        """
        Get feature importance from SHAP values.
        
        Returns:
        - DataFrame with features and importance scores
        """
        if self.shap_values is None:
            raise ValueError('SHAP values not computed yet')
        
        # Calculate mean absolute SHAP values
        if len(self.shap_values.shape) == 2:
            importance = np.abs(self.shap_values).mean(axis=0)
        else:
            importance = np.abs(self.shap_values).mean(axis=(0, 2))
        
        importance_df = pd.DataFrame({
            'feature': self.feature_names,
            'importance': importance
        }).sort_values('importance', ascending=False)
        
        return importance_df
    
    def plot_summary(self):
        """
        Generate Plotly summary plot (bar chart of feature importance).
        
        Returns:
        - Plotly figure
        """
        if self.shap_values is None:
            raise ValueError('SHAP values not computed yet')
        
        importance_df = self.get_feature_importance().head(10)
        
        fig = go.Figure(data=[
            go.Bar(
                x=importance_df['importance'],
                y=importance_df['feature'],
                orientation='h',
                marker=dict(
                    color=importance_df['importance'],
                    colorscale='Viridis'
                )
            )
        ])
        
        fig.update_layout(
            title='Top 10 Feature Importance (Mean |SHAP value|)',
            xaxis_title='Mean |SHAP value|',
            yaxis_title='Feature',
            height=500,
            yaxis={'categoryorder': 'total ascending'}
        )
        
        return fig