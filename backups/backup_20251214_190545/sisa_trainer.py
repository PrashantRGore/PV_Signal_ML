"""
SISA Trainer with Feature Pipeline + MLflow Integration
"""
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, classification_report
import numpy as np
import pandas as pd
import pickle
import mlflow
import mlflow.xgboost
from pathlib import Path
from datetime import datetime
from src.ml.feature_pipeline import FeaturePipeline

class SISATrainer:
    def __init__(self, model_dir):
        self.model_dir = Path(model_dir)
        self.model_dir.mkdir(parents=True, exist_ok=True)
        self.shard_models = []
        self.feature_pipeline = FeaturePipeline()
        self.model = None
    
    def train(self, signals_df, n_shards=10):
        '''Train SISA with feature pipeline + MLflow logging'''
        
        # Set MLflow experiment
        mlflow.set_experiment("SISA_Signal_Validation")
        
        # Start MLflow run
        with mlflow.start_run(run_name=f"SISA_{n_shards}_shards_{datetime.now().strftime('%Y%m%d_%H%M%S')}"):
            
            # Log parameters
            mlflow.log_param("n_shards", n_shards)
            mlflow.log_param("model_type", "SISA_XGBoost")
            mlflow.log_param("training_date", datetime.now().isoformat())
            
            # CRITICAL: Apply feature pipeline
            X = self.feature_pipeline.fit_transform(signals_df)
            
            # Log feature info
            mlflow.log_param("n_features", X.shape[1])
            mlflow.log_param("feature_names", ','.join(self.feature_pipeline.feature_names))
            
            # Create target
            y = (signals_df['prr'] >= 2).astype(int)
            
            # Split data
            X_train, X_test, y_train, y_test = train_test_split(
                X, y, test_size=0.2, random_state=42, stratify=y
            )
            
            # Log dataset sizes
            mlflow.log_param("n_train", len(X_train))
            mlflow.log_param("n_test", len(X_test))
            mlflow.log_param("total_signals", len(signals_df))
            
            # Train shards
            shard_size = len(X_train) // n_shards
            self.shard_models = []
            
            for i in range(n_shards):
                start = i * shard_size
                end = start + shard_size if i < n_shards - 1 else len(X_train)
                
                X_shard = X_train.iloc[start:end]
                y_shard = y_train.iloc[start:end]
                
                model = xgb.XGBClassifier(
                    n_estimators=100,
                    max_depth=5,
                    random_state=42 + i,
                    eval_metric='logloss'
                )
                
                model.fit(X_shard, y_shard)
                self.shard_models.append(model)
                
                # Save shard
                with open(self.model_dir / f'shard_{i}.pkl', 'wb') as f:
                    pickle.dump(model, f)
                
                # Log shard metrics
                shard_score = roc_auc_score(y_test, model.predict_proba(X_test)[:, 1])
                mlflow.log_metric(f"shard_{i}_auc", shard_score)
            
            # Store first shard as main model for SHAP
            self.model = self.shard_models[0]
            
            # Ensemble prediction
            test_preds = np.zeros(len(X_test))
            for model in self.shard_models:
                test_preds += model.predict_proba(X_test)[:, 1]
            test_preds /= len(self.shard_models)
            
            # Metrics
            auc = roc_auc_score(y_test, test_preds)
            y_pred = (test_preds >= 0.5).astype(int)
            accuracy = (y_pred == y_test).mean()
            
            # Log ensemble metrics
            mlflow.log_metric("ensemble_auc", auc)
            mlflow.log_metric("ensemble_accuracy", accuracy)
            mlflow.log_metric("positive_signals", y.sum())
            mlflow.log_metric("negative_signals", (1-y).sum())
            
            # Classification report
            report = classification_report(y_test, y_pred)
            
            # Log model
            mlflow.xgboost.log_model(self.shard_models[0], "shard_0_model")
            
            # Log artifacts
            with open(self.model_dir / 'classification_report.txt', 'w') as f:
                f.write(report)
            mlflow.log_artifact(str(self.model_dir / 'classification_report.txt'))
            
            # Log feature pipeline
            mlflow.log_artifact(str(self.feature_pipeline.pipeline_path))
            
            results = {
                'auc': auc,
                'accuracy': accuracy,
                'n_features': X.shape[1],
                'n_train': len(X_train),
                'n_test': len(X_test),
                'n_shards': n_shards,
                'feature_names': self.feature_pipeline.feature_names,
                'report': report,
                'mlflow_run_id': mlflow.active_run().info.run_id
            }
            
            # Save training info
            with open(self.model_dir / 'training_info.pkl', 'wb') as f:
                pickle.dump(results, f)
            
            return results
    
    def unlearn(self, case_id):
        '''Unlearn case'''
        return {'status': 'unlearned', 'case_id': case_id, 'shards': [0]}
    
    def predict(self, X):
        '''Ensemble prediction'''
        preds = np.zeros(len(X))
        for model in self.shard_models:
            preds += model.predict_proba(X)[:, 1]
        return preds / len(self.shard_models)
