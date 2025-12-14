"""
Enhanced SISA Trainer with Progress Tracking
Shows step-by-step unlearning progress
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
import streamlit as st
import time

class SISATrainer:
    def __init__(self, model_dir):
        self.model_dir = Path(model_dir)
        self.model_dir.mkdir(parents=True, exist_ok=True)
        self.shard_models = []
        self.feature_pipeline = FeaturePipeline()
        self.model = None
        self.training_data = None  # Store for unlearning
        self.shard_assignments = {}  # Track which cases are in which shards
    
    def train(self, signals_df, n_shards=10):
        '''Train SISA with feature pipeline + MLflow logging'''
        
        # Set MLflow experiment
        mlflow.set_experiment("2_ML_Model_Training")
        
        # Start MLflow run
        with mlflow.start_run(run_name=f"SISA_{n_shards}_shards_{datetime.now().strftime('%Y%m%d_%H%M%S')}"):
            
            # Log parameters
            mlflow.log_param("n_shards", n_shards)
            mlflow.log_param("model_type", "SISA_XGBoost")
            mlflow.log_param("training_date", datetime.now().isoformat())
            
            # Apply feature pipeline
            X = self.feature_pipeline.fit_transform(signals_df)
            
            # Log feature info
            mlflow.log_param("n_features", X.shape[1])
            mlflow.log_param("feature_names", ','.join(self.feature_pipeline.feature_names))
            
            # Create target
            y = (signals_df['prr'] >= 2).astype(int)
            
            # Create case IDs for tracking
            case_ids = signals_df['drug_name'] + '_' + signals_df['event_term']
            
            # Split data
            X_train, X_test, y_train, y_test, ids_train, ids_test = train_test_split(
                X, y, case_ids, test_size=0.2, random_state=42, stratify=y
            )
            
            # Store training data for unlearning
            self.training_data = pd.DataFrame(X_train, columns=self.feature_pipeline.feature_names)
            self.training_data['target'] = y_train.values
            self.training_data['case_id'] = ids_train.values
            
            # Log dataset sizes
            mlflow.log_param("n_train", len(X_train))
            mlflow.log_param("n_test", len(X_test))
            mlflow.log_param("total_signals", len(signals_df))
            
            # Train shards
            shard_size = len(X_train) // n_shards
            self.shard_models = []
            self.shard_assignments = {}
            
            for i in range(n_shards):
                start = i * shard_size
                end = start + shard_size if i < n_shards - 1 else len(X_train)
                
                X_shard = X_train.iloc[start:end]
                y_shard = y_train.iloc[start:end]
                shard_ids = ids_train.iloc[start:end]
                
                # Track case assignments
                for case_id in shard_ids:
                    self.shard_assignments[case_id] = i
                
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
            
            # Save metadata
            metadata = {
                'shard_assignments': self.shard_assignments,
                'n_shards': n_shards,
                'training_date': datetime.now().isoformat()
            }
            metadata_file = self.model_dir / 'metadata.pkl'
            with open(metadata_file, 'wb') as f:
                pickle.dump(metadata, f)
            mlflow.log_artifact(str(metadata_file))
            
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
            
            return results
    
    def unlearn(self, case_id, progress_callback=None):
        '''
        Unlearn case with PROGRESS TRACKING
        
        Architecture Flow:
        1. Identify case in training data
        2. Log unlearning event to MLflow
        3. Identify affected shard
        4. Remove case from shard data
        5. Retrain affected shard
        6. Update ensemble
        7. Log updated model
        8. Complete unlearning
        '''
        
        try:
            # Step 1: Identify case
            if progress_callback:
                progress_callback(0.0, "🔍 Identifying case in training data...")
            time.sleep(0.5)
            
            if case_id not in self.shard_assignments:
                if progress_callback:
                    progress_callback(1.0, f"⚠️ Case {case_id} not found in training data")
                return {'status': 'not_found', 'case_id': case_id}
            
            affected_shard = self.shard_assignments[case_id]
            
            if progress_callback:
                progress_callback(0.15, f"✅ Case found in Shard {affected_shard}")
            time.sleep(0.5)
            
            # Step 2: Log to MLflow
            if progress_callback:
                progress_callback(0.25, "📊 Logging unlearning event to MLflow...")
            time.sleep(0.5)
            
            mlflow.set_experiment("5_Data_Unlearning")
            
            with mlflow.start_run(run_name=f"Unlearning_{case_id}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"):
                
                mlflow.set_tag("process.type", "data_unlearning")
                mlflow.set_tag("regulatory.compliance", "GDPR_Article17_HIPAA")
                mlflow.log_param("case_id", case_id)
                mlflow.log_param("affected_shard", affected_shard)
                mlflow.log_param("unlearning_date", datetime.now().isoformat())
                
                if progress_callback:
                    progress_callback(0.35, "✅ MLflow logging initiated")
                time.sleep(0.5)
                
                # Step 3: Remove from shard data
                if progress_callback:
                    progress_callback(0.45, f"🗑️ Removing case from Shard {affected_shard}...")
                time.sleep(0.5)
                
                # Filter training data
                shard_data = self.training_data[
                    self.training_data['case_id'].isin(
                        [k for k, v in self.shard_assignments.items() if v == affected_shard]
                    )
                ]
                shard_data_filtered = shard_data[shard_data['case_id'] != case_id]
                
                original_size = len(shard_data)
                new_size = len(shard_data_filtered)
                
                mlflow.log_metric("original_shard_size", original_size)
                mlflow.log_metric("new_shard_size", new_size)
                
                if progress_callback:
                    progress_callback(0.55, f"✅ Case removed ({original_size} → {new_size} samples)")
                time.sleep(0.5)
                
                # Step 4: Retrain affected shard
                if progress_callback:
                    progress_callback(0.65, f"🔄 Re-training Shard {affected_shard}...")
                time.sleep(0.5)
                
                X_shard = shard_data_filtered.drop(['target', 'case_id'], axis=1)
                y_shard = shard_data_filtered['target']
                
                new_model = xgb.XGBClassifier(
                    n_estimators=100,
                    max_depth=5,
                    random_state=42 + affected_shard,
                    eval_metric='logloss'
                )
                
                new_model.fit(X_shard, y_shard)
                
                if progress_callback:
                    progress_callback(0.80, f"✅ Shard {affected_shard} retrained successfully")
                time.sleep(0.5)
                
                # Step 5: Update ensemble
                if progress_callback:
                    progress_callback(0.90, "🔄 Updating model ensemble...")
                time.sleep(0.5)
                
                self.shard_models[affected_shard] = new_model
                
                # Save updated shard
                with open(self.model_dir / f'shard_{affected_shard}.pkl', 'wb') as f:
                    pickle.dump(new_model, f)
                
                # Update assignments
                del self.shard_assignments[case_id]
                
                # Step 6: Log updated model
                if progress_callback:
                    progress_callback(0.95, "📋 Logging updated model to MLflow...")
                time.sleep(0.5)
                
                mlflow.xgboost.log_model(new_model, f"updated_shard_{affected_shard}")
                mlflow.set_tag("unlearning.status", "complete")
                
                run_id = mlflow.active_run().info.run_id
                mlflow.log_param("mlflow_run_id", run_id)
                
                if progress_callback:
                    progress_callback(1.0, "✅ Unlearning complete!")
                time.sleep(0.5)
                
                return {
                    'status': 'unlearned',
                    'case_id': case_id,
                    'affected_shard': affected_shard,
                    'original_shard_size': original_size,
                    'new_shard_size': new_size,
                    'mlflow_run_id': run_id
                }
        
        except Exception as e:
            if progress_callback:
                progress_callback(1.0, f"❌ Error: {str(e)}")
            return {'status': 'error', 'error': str(e)}
    
    def predict(self, X):
        '''Ensemble prediction'''
        preds = np.zeros(len(X))
        for model in self.shard_models:
            preds += model.predict_proba(X)[:, 1]
        return preds / len(self.shard_models)
