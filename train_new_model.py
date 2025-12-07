#!/usr/bin/env python3
"""
Train a new ML model and log to MLflow

This script:
1. Loads sample FAERS data
2. Trains an XGBoost model
3. Logs to MLflow with today's date
4. Saves metrics and artifacts

Usage:
    python train_new_model.py
"""

import os
import sys
import json
import sqlite3
import numpy as np
import pandas as pd
from datetime import datetime
from pathlib import Path

# MLflow
import mlflow
import mlflow.xgboost
from mlflow.models import infer_signature

# ML Libraries
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score
from sklearn.preprocessing import LabelEncoder

def setup_mlflow():
    """Configure MLflow tracking."""
    mlflow_uri = "file:///./mlruns"
    mlflow.set_tracking_uri(mlflow_uri)
    mlflow.set_experiment("pv-signal-detection")
    print(f"✅ MLflow configured: {mlflow_uri}")

def load_sample_data():
    """Load sample FAERS data from database or create synthetic data."""
    print("\n📊 Loading Sample Data...")
    
    db_path = Path("pv_signal.db")
    
    if db_path.exists():
        try:
            conn = sqlite3.connect(db_path)
            query = "SELECT * FROM icsr LIMIT 100"
            df = pd.read_sql_query(query, conn)
            conn.close()
            
            if len(df) > 0:
                print(f"✅ Loaded {len(df)} records from database")
                return df
        except Exception as e:
            print(f"⚠️  Could not load from database: {e}")
    
    # Create synthetic data if database is empty
    print("📝 Creating synthetic sample data...")
    
    np.random.seed(42)
    n_samples = 100
    
    drugs = ["Aspirin", "Ibuprofen", "Acetaminophen", "Naproxen", "Diclofenac"]
    reactions = ["Headache", "Nausea", "Dizziness", "Rash", "Fever", "Fatigue"]
    outcomes = ["Recovered", "Recovering", "Not recovered", "Fatal", "Unknown"]
    
    data = {
        "icsr_id": [f"ICSR_2025_{i:04d}" for i in range(n_samples)],
        "drug_name": np.random.choice(drugs, n_samples),
        "reaction": np.random.choice(reactions, n_samples),
        "outcome": np.random.choice(outcomes, n_samples),
        "report_date": pd.date_range("2025-12-01", periods=n_samples, freq="H"),
        "severity": np.random.choice([1, 2, 3, 4, 5], n_samples),  # 1-5 scale
        "age": np.random.randint(18, 85, n_samples),
        "gender": np.random.choice(["M", "F"], n_samples),
    }
    
    df = pd.DataFrame(data)
    print(f"✅ Created {len(df)} synthetic records")
    return df

def prepare_features(df):
    """Prepare features for ML model."""
    print("\n🔧 Preparing Features...")
    
    # Create a target variable (signal detection: 1 if serious, 0 otherwise)
    serious_outcomes = ["Fatal", "Not recovered"]
    df["is_signal"] = df["outcome"].isin(serious_outcomes).astype(int)
    
    # Encode categorical variables
    le_drug = LabelEncoder()
    le_reaction = LabelEncoder()
    le_gender = LabelEncoder()
    
    df["drug_encoded"] = le_drug.fit_transform(df["drug_name"])
    df["reaction_encoded"] = le_reaction.fit_transform(df["reaction"])
    df["gender_encoded"] = le_gender.fit_transform(df["gender"])
    
    # Select features for model
    feature_cols = ["drug_encoded", "reaction_encoded", "severity", "age", "gender_encoded"]
    X = df[feature_cols]
    y = df["is_signal"]
    
    print(f"✅ Features prepared: {len(feature_cols)} features, {len(X)} samples")
    print(f"   Signal distribution: {y.value_counts().to_dict()}")
    
    return X, y, feature_cols

def train_model(X, y):
    """Train XGBoost model."""
    print("\n🤖 Training XGBoost Model...")
    
    # Split data
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42
    )
    
    print(f"   Training set: {len(X_train)} samples")
    print(f"   Test set: {len(X_test)} samples")
    
    # Train model
    model = xgb.XGBClassifier(
        n_estimators=100,
        max_depth=5,
        learning_rate=0.1,
        random_state=42,
        verbosity=0
    )
    
    model.fit(X_train, y_train)
    print("✅ Model trained")
    
    # Evaluate
    y_pred = model.predict(X_test)
    y_pred_proba = model.predict_proba(X_test)[:, 1]
    
    metrics = {
        "accuracy": accuracy_score(y_test, y_pred),
        "precision": precision_score(y_test, y_pred, zero_division=0),
        "recall": recall_score(y_test, y_pred, zero_division=0),
        "f1": f1_score(y_test, y_pred, zero_division=0),
        "roc_auc": roc_auc_score(y_test, y_pred_proba),
    }
    
    print("✅ Model evaluated:")
    for metric, value in metrics.items():
        print(f"   {metric}: {value:.4f}")
    
    return model, metrics, X_test, y_test, y_pred

def log_to_mlflow(model, metrics, X_test, y_test, y_pred, feature_cols):
    """Log model and metrics to MLflow."""
    print("\n📝 Logging to MLflow...")
    
    with mlflow.start_run(run_name=f"xgboost_run_{datetime.now().strftime('%Y%m%d_%H%M%S')}"):
        # Log parameters
        mlflow.log_params({
            "model_type": "XGBoost",
            "n_estimators": 100,
            "max_depth": 5,
            "learning_rate": 0.1,
            "test_size": 0.2,
        })
        
        # Log metrics
        for metric_name, metric_value in metrics.items():
            mlflow.log_metric(metric_name, metric_value)
        
        # Log model
        signature = infer_signature(X_test, y_pred)
        mlflow.xgboost.log_model(
            model,
            "model",
            signature=signature,
            input_example=X_test.iloc[:5]
        )
        
        # Log feature importance
        feature_importance = pd.DataFrame({
            "feature": feature_cols,
            "importance": model.feature_importances_
        }).sort_values("importance", ascending=False)
        
        mlflow.log_table(
            feature_importance,
            artifact_file="feature_importance.json"
        )
        
        # Log additional artifacts
        artifacts = {
            "model_info": {
                "model_type": "XGBoost",
                "n_features": len(feature_cols),
                "features": feature_cols,
                "training_date": datetime.now().isoformat(),
            },
            "metrics": metrics,
            "feature_importance": feature_importance.to_dict(orient="records"),
        }
        
        artifacts_json = json.dumps(artifacts, indent=2, default=str)
        with open("mlflow_artifacts.json", "w") as f:
            f.write(artifacts_json)
        
        mlflow.log_artifact("mlflow_artifacts.json")
        
        print("✅ Model logged to MLflow")
        print(f"   Run ID: {mlflow.active_run().info.run_id}")
        print(f"   Experiment: pv-signal-detection")

def main():
    """Main training pipeline."""
    print("\n" + "="*70)
    print("🚀 Train New ML Model and Log to MLflow")
    print("="*70)
    
    try:
        # Setup
        setup_mlflow()
        
        # Load data
        df = load_sample_data()
        
        # Prepare features
        X, y, feature_cols = prepare_features(df)
        
        # Train model
        model, metrics, X_test, y_test, y_pred = train_model(X, y)
        
        # Log to MLflow
        log_to_mlflow(model, metrics, X_test, y_test, y_pred, feature_cols)
        
        print("\n" + "="*70)
        print("✅ Training Complete!")
        print("="*70)
        print("\n📊 View results at: http://localhost:5000")
        print("   Look for experiment: pv-signal-detection")
        print("   Latest run: Today's date\n")
        
        return 0
        
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())
