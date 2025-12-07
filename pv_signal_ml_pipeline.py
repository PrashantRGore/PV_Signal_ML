"""
PV-Signal-ML end-to-end baseline pipeline (audit-ready)
1. Load features_2020_prr_labeled.parquet
2. Train/test split (March train, June test)
3. XGBoost + MLflow
4. Save ranked June signals with evidence columns
5. Log run metadata for audits
"""

from pathlib import Path
import pandas as pd
import numpy as np
import xgboost as xgb
from sklearn.metrics import average_precision_score, roc_auc_score
import mlflow

from run_metadata import log_run_metadata
from audit_logging import audit_logger

PIPELINE_VERSION = "1.0.0"

# Configuration
BASE_DIR = Path(r"C:\Users\koreo\Downloads\pv-signal-ml")
DATA_DIR = BASE_DIR / "data"
PROCESSED_DIR = DATA_DIR / "processed"
RESULTS_DIR = BASE_DIR / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

FEATURES_PATH = PROCESSED_DIR / "features_2020_prr_labeled.parquet"

MLFLOW_URI = "file:///C:/Users/koreo/mlruns"
MLFLOW_EXPERIMENT = "pv-signal-ml-prr"


print("[INFO] Loading features...")
features = pd.read_parquet(FEATURES_PATH)
features["reference_date"] = pd.to_datetime(features["reference_date"])
print("[INFO] Features shape:", features.shape)

# 2. Train/test split
train_cut = pd.Timestamp("2020-03-01")
test_month = pd.Timestamp("2020-06-01")

train_mask = features["reference_date"] <= train_cut
test_mask = features["reference_date"] == test_month

exclude_cols = {
    "DRUGNAME_NORM",
    "pt",
    "reference_date",
    "future_A_6m",
    "future_N_6m",
    "future_prr_6m",
    "label_prr",
}
feature_cols = [c for c in features.columns if c not in exclude_cols]

X_train = features.loc[train_mask, feature_cols]
y_train = features.loc[train_mask, "label_prr"]
X_test = features.loc[test_mask, feature_cols]
y_test = features.loc[test_mask, "label_prr"]

print("[INFO] Train size:", X_train.shape, "Test size:", X_test.shape)
print("[INFO] Train labels:\n", y_train.value_counts())
print("[INFO] Test labels:\n", y_test.value_counts())

# 3. Train XGBoost + MLflow
mlflow.end_run()
mlflow.set_tracking_uri(MLFLOW_URI)
mlflow.set_experiment(MLFLOW_EXPERIMENT)

with mlflow.start_run():
    params = {
        "objective": "binary:logistic",
        "eval_metric": "aucpr",
        "max_depth": 6,
        "eta": 0.05,
        "subsample": 0.8,
        "colsample_bytree": 0.8,
        "lambda": 1.0,
        "alpha": 0.0,
        "n_estimators": 500,
        "tree_method": "hist",
        "base_score": 0.5,
    }

    print("[INFO] Training XGBoost...")
    model = xgb.XGBClassifier(**params)
    model.fit(X_train, y_train)

    y_pred_proba = model.predict_proba(X_test)[:, 1]

    if len(y_test.unique()) < 2:
        print("[WARN] Test has single class; AP/AUC=NaN (expected for rare events)")
        ap = float("nan")
        auc = float("nan")
    else:
        ap = average_precision_score(y_test, y_pred_proba)
        auc = roc_auc_score(y_test, y_pred_proba)

    mlflow.log_params(params)
    mlflow.log_param("pipeline_version", PIPELINE_VERSION)
    mlflow.log_metric("AP_test", float(ap))
    mlflow.log_metric("AUC_test", float(auc))
    mlflow.xgboost.log_model(model, "model")

    print("[INFO] AP_test:", ap)
    print("[INFO] AUC_test:", auc)

    # 4. Save ranked signals with evidence columns
    print("[INFO] Saving ranked June signals...")
    cols = [
        "DRUGNAME_NORM",
        "pt",
        "label_prr",
        "A_6m",
        "N_6m",
        "prr_6m",
        "A_12m",
        "N_12m",
    ]
    test_df = features.loc[test_mask, cols].copy()
    test_df["score"] = y_pred_proba
    test_df["rank"] = test_df["score"].rank(ascending=False)
    test_df["signal_period"] = "2020-06"
    test_df["model_version"] = "xgboost_prr_v1"

    ranked = test_df.sort_values("score", ascending=False).reset_index(drop=True)
    ranked_path = RESULTS_DIR / "ranked_signals_2020Q2.csv"
    ranked.to_csv(ranked_path, index=False)
    print("[INFO] Saved:", ranked_path)
    print("\n[INFO] Top 10 predictions:")
    print(ranked.head(10)[["DRUGNAME_NORM", "pt", "score", "label_prr"]])

    # 5. Log run metadata for audits
    run_id = mlflow.active_run().info.run_id
    label_def = "label_prr = 1 if prr_6m >= 0.1 and A_6m >= 3"
    data_period = "FAERS 2020, signals in 2020-06"

    log_run_metadata(
        model_version="xgboost_prr_v1",
        data_period=data_period,
        label_definition=label_def,
        feature_list=feature_cols,
        train_size=len(X_train),
        test_size=len(X_test),
        train_pos=int(y_train.sum()),
        test_pos=int(y_test.sum()),
        ap_test=ap,
        auc_test=auc,
        mlflow_run_id=run_id,
    )
    
    # 6. Log MLflow run to audit trail (FDA 21 CFR Part 11 compliance)
    audit_logger.log_mlflow_run(
        run_id=run_id,
        experiment_name=MLFLOW_EXPERIMENT,
        model_version="xgboost_prr_v1",
        metrics={
            "AP_test": float(ap) if not (ap != ap) else None,
            "AUC_test": float(auc) if not (auc != auc) else None,
            "train_size": len(X_train),
            "test_size": len(X_test),
        },
        parameters=params,
        user_id="system",
        status="completed"
    )
    print(f"✅ MLflow run {run_id} logged to audit trail")

print("\n🎉 PIPELINE COMPLETE!")
print("✅ MLflow: mlflow ui")
print("✅ Ranked signals:", ranked_path)
print("✅ Run metadata: results/run_summaries.jsonl")