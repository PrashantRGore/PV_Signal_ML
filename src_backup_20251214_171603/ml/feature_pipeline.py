"""
Feature Pipeline Manager
Ensures consistent features across signal detection, SISA training, and SHAP
"""
import pandas as pd
import numpy as np
import pickle
from pathlib import Path

class FeaturePipeline:
    def __init__(self):
        self.feature_names = None
        self.feature_transformations = {}
        self.pipeline_path = Path('models/sisa/feature_pipeline.pkl')
        self.pipeline_path.parent.mkdir(parents=True, exist_ok=True)
    
    def fit_transform(self, signals_df):
        '''
        Transform raw signals into ML-ready features
        This is called DURING SISA training to establish the feature pipeline
        '''
        df = signals_df.copy()
        
        # Ensure required columns exist
        required_cols = ['case_count', 'prr', 'chi2']
        for col in required_cols:
            if col not in df.columns:
                # Try alternative names
                if col == 'case_count' and 'count' in df.columns:
                    df['case_count'] = df['count']
                elif col == 'chi2' and 'chi_square' in df.columns:
                    df['chi2'] = df['chi_square']
        
        # Create transformed features
        df['log_prr'] = np.log1p(df['prr'].fillna(0))
        df['log_chi2'] = np.log1p(df['chi2'].fillna(0))
        
        # Calculate ROR if not present
        if 'ror' not in df.columns:
            df['ror'] = df['prr']  # Simplified; actual ROR calculation if needed
        
        # Select final feature set
        feature_cols = ['case_count', 'prr', 'chi2', 'log_prr', 'log_chi2']
        
        # Store feature names for later use
        self.feature_names = feature_cols
        
        # Create feature matrix
        X = df[feature_cols].fillna(0).replace([np.inf, -np.inf], 0)
        
        # Save pipeline
        self.save()
        
        return X, df
    
    def transform(self, signals_df):
        '''
        Transform new data using the saved pipeline
        This is called DURING SHAP analysis to match training features
        '''
        if self.feature_names is None:
            self.load()
        
        df = signals_df.copy()
        
        # Ensure required columns
        if 'case_count' not in df.columns and 'count' in df.columns:
            df['case_count'] = df['count']
        if 'chi2' not in df.columns and 'chi_square' in df.columns:
            df['chi2'] = df['chi_square']
        
        # Apply same transformations as training
        df['log_prr'] = np.log1p(df['prr'].fillna(0))
        df['log_chi2'] = np.log1p(df['chi2'].fillna(0))
        
        if 'ror' not in df.columns:
            df['ror'] = df['prr']
        
        # Use exact same features as training
        X = df[self.feature_names].fillna(0).replace([np.inf, -np.inf], 0)
        
        return X
    
    def save(self):
        '''Save feature pipeline to disk'''
        pipeline_data = {
            'feature_names': self.feature_names,
            'transformations': self.feature_transformations
        }
        with open(self.pipeline_path, 'wb') as f:
            pickle.dump(pipeline_data, f)
    
    def load(self):
        '''Load feature pipeline from disk'''
        if self.pipeline_path.exists():
            with open(self.pipeline_path, 'rb') as f:
                pipeline_data = pickle.load(f)
                self.feature_names = pipeline_data['feature_names']
                self.feature_transformations = pipeline_data['transformations']
            return True
        return False
    
    def get_feature_info(self):
        '''Get information about the feature pipeline'''
        return {
            'feature_names': self.feature_names,
            'n_features': len(self.feature_names) if self.feature_names else 0,
            'pipeline_exists': self.pipeline_path.exists()
        }
