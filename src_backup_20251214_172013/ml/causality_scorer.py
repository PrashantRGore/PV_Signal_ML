"""
Causality Assessment: Apply ML model to rank signals
"""

import pandas as pd
import pickle
from pathlib import Path
import streamlit as st
from src.utils.logger import setup_logger

logger = setup_logger(__name__)

class CausalityScorer:
    """
    Applies trained ML model to score signal causality
    Integrates with SISA sharding for GDPR compliance
    """
    
    def __init__(self, model_path='models/sisa_model.pkl'):
        self.model_path = Path(model_path)
        self.model = None
    
    def load_model(self):
        """Load trained SISA model"""
        if self.model_path.exists():
            with open(self.model_path, 'rb') as f:
                self.model = pickle.load(f)
            logger.info(f"Loaded causality model from {self.model_path}")
            return True
        else:
            logger.warning(f"Model not found: {self.model_path}")
            return False
    
    def prepare_features(self, signals_df):
        """Prepare features for ML model"""
        features = pd.DataFrame({
            'prr': signals_df['prr'],
            'chi2': signals_df['chi2'],
            'case_count': signals_df['case_count'],
            'ror': signals_df.get('ror', 0)
        })
        return features
    
    def score_signals(self, signals_df):
        """Apply ML model to score causality probability"""
        if self.model is None:
            if not self.load_model():
                st.warning("⚠️ ML model not available. Using statistical signals only.")
                signals_df['causality_score'] = 0.5
                return signals_df
        
        try:
            X = self.prepare_features(signals_df)
            
            if hasattr(self.model, 'predict_proba'):
                probs = self.model.predict_proba(X)
                causality_scores = probs[:, 1] if probs.shape[1] > 1 else probs[:, 0]
            else:
                causality_scores = self.model.predict(X)
            
            signals_df['causality_score'] = causality_scores
            
            signals_df['causality_level'] = pd.cut(
                causality_scores,
                bins=[0, 0.3, 0.7, 1.0],
                labels=['Low', 'Moderate', 'High']
            )
            
            logger.info(f"Scored {len(signals_df)} signals for causality")
            return signals_df
            
        except Exception as e:
            logger.error(f"Error scoring causality: {e}")
            st.error(f"Causality scoring failed: {e}")
            signals_df['causality_score'] = 0.5
            return signals_df
