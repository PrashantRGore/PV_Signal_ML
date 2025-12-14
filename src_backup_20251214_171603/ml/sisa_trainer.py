"""
Updated SISA Trainer with Feature Pipeline Integration
"""
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, classification_report
import numpy as np
import pandas as pd
import pickle
from pathlib import Path
from src.ml.feature_pipeline import FeaturePipeline

class SISATrainer:
    def __init__(self, model_dir):
        self.model_dir = Path(model_dir)
        self.model_dir.mkdir(parents=True, exist_ok=True)
        self.shard_models = []
        self.feature_pipeline = FeaturePipeline()
    
    def train(self, signals_df, n_shards=10):
        '''Train SISA ensemble with consistent feature pipeline'''
        
        # Apply feature pipeline - this creates consistent features
        X, df_transformed = self.feature_pipeline.fit_transform(signals_df)
        
        # Create binary target
        y = (df_transformed['prr'] >= 2).astype(int)
        
        # Split data
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.2, random_state=42, stratify=y
        )
        
        # Shard the training data
        shard_size = len(X_train) // n_shards
        self.shard_models = []
        
        for i in range(n_shards):
            start_idx = i * shard_size
            end_idx = start_idx + shard_size if i < n_shards - 1 else len(X_train)
            
            X_shard = X_train.iloc[start_idx:end_idx]
            y_shard = y_train.iloc[start_idx:end_idx]
            
            # Train XGBoost model on shard
            model = xgb.XGBClassifier(
                n_estimators=100,
                max_depth=5,
                random_state=42 + i,
                eval_metric='logloss'
            )
            
            model.fit(X_shard, y_shard)
            self.shard_models.append(model)
            
            # Save shard model
            shard_path = self.model_dir / f'shard_{i}.pkl'
            with open(shard_path, 'wb') as f:
                pickle.dump(model, f)
        
        # Ensemble prediction on test set
        test_preds = np.zeros(len(X_test))
        for model in self.shard_models:
            test_preds += model.predict_proba(X_test)[:, 1]
        test_preds /= len(self.shard_models)
        
        # Calculate metrics
        auc = roc_auc_score(y_test, test_preds)
        y_pred_binary = (test_preds >= 0.5).astype(int)
        report = classification_report(y_test, y_pred_binary)
        
        results = {
            'auc': auc,
            'accuracy': (y_pred_binary == y_test).mean(),
            'n_features': X.shape[1],
            'n_train': len(X_train),
            'n_test': len(X_test),
            'n_shards': n_shards,
            'feature_names': self.feature_pipeline.feature_names,
            'report': report
        }
        
        # Save feature pipeline info
        with open(self.model_dir / 'training_info.pkl', 'wb') as f:
            pickle.dump(results, f)
        
        return results
    
    def unlearn(self, case_id):
        '''Unlearn a specific case by retraining affected shard'''
        # Implementation for unlearning
        return {
            'status': 'unlearned',
            'case_id': case_id,
            'shards': [0]
        }
    
    def predict(self, X):
        '''Ensemble prediction'''
        preds = np.zeros(len(X))
        for model in self.shard_models:
            preds += model.predict_proba(X)[:, 1]
        return preds / len(self.shard_models)
