"""
Feature Pipeline for SISA Training
Ensures consistent features between training and inference
FIXED: Removed erroneous self.save() call from fit_transform()
"""

import pandas as pd
import numpy as np
import pickle
from pathlib import Path
from typing import List

class FeaturePipeline:
    """
    Feature pipeline for signal detection ML models
    Ensures consistent feature engineering across training and inference
    """
    
    def __init__(self):
        self.feature_names = []
        self.fitted = False
    
    def fit_transform(self, df: pd.DataFrame) -> np.ndarray:
        """
        Fit the pipeline and transform the data
        
        Args:
            df: DataFrame with signal data
        
        Returns:
            Numpy array of features
        """
        # Define features to use
        self.feature_names = ['case_count', 'prr', 'chi2']
        
        # Extract features
        X = df[self.feature_names].values
        
        # Handle missing values
        X = np.nan_to_num(X, nan=0.0, posinf=1e6, neginf=-1e6)
        
        self.fitted = True
        
        # DO NOT CALL self.save() here - pipeline is saved by sisa_trainer
        return X
    
    def transform(self, df: pd.DataFrame) -> np.ndarray:
        """
        Transform new data using fitted pipeline
        
        Args:
            df: DataFrame with signal data
        
        Returns:
            Numpy array of features
        """
        if not self.fitted:
            raise ValueError("Pipeline must be fitted before transform")
        
        # Extract same features
        X = df[self.feature_names].values
        
        # Handle missing values
        X = np.nan_to_num(X, nan=0.0, posinf=1e6, neginf=-1e6)
        
        return X
    
    def save(self, filepath: str):
        """
        Save the feature pipeline to disk
        
        Args:
            filepath: Path to save the pipeline (string or Path object)
        """
        filepath = Path(filepath)
        filepath.parent.mkdir(parents=True, exist_ok=True)
        
        with open(filepath, 'wb') as f:
            pickle.dump({
                'feature_names': self.feature_names,
                'fitted': self.fitted
            }, f)
    
    def load(self, filepath: str):
        """
        Load the feature pipeline from disk
        
        Args:
            filepath: Path to load the pipeline from
        """
        with open(filepath, 'rb') as f:
            data = pickle.load(f)
            self.feature_names = data['feature_names']
            self.fitted = data['fitted']
    
    @staticmethod
    def from_file(filepath: str) -> 'FeaturePipeline':
        """
        Load pipeline from file
        
        Args:
            filepath: Path to pipeline file
        
        Returns:
            Loaded FeaturePipeline instance
        """
        pipeline = FeaturePipeline()
        pipeline.load(filepath)
        return pipeline
