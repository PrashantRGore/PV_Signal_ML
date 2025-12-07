"""
SHAP Analysis - RETRAIN + SHAP (100% guaranteed to work)
No MLflow model loading issues
"""

import pandas as pd
import numpy as np
import xgboost as xgb
from sklearn.metrics import average_precision_score, roc_auc_score
import shap
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

# Config
BASE_DIR = Path(r"C:\Users\koreo\Downloads\pv-signal-ml")
DATA_DIR = BASE_DIR / "data" / "processed"
RESULTS_DIR = BASE_DIR / "results"
SHAP_DIR = RESULTS_DIR / "shap_analysis"
SHAP_DIR.mkdir(parents=True, exist_ok=True)

FEATURES_PATH = DATA_DIR / "features_2020_prr_labeled.parquet"

print("🔍 SHAP Analysis (SIMPLE VERSION)")
print("="*50)

# 1. Load + split data (same as pipeline)
print("[INFO] Loading features...")
features = pd.read_parquet(FEATURES_PATH)
features["reference_date"] = pd.to_datetime(features["reference_date"])

train_cut = pd.Timestamp("2020-03-01")
test_month = pd.Timestamp("2020-06-01")

train_mask = features["reference_date"] <= train_cut
test_mask = features["reference_date"] == test_month

exclude_cols = {"DRUGNAME_NORM", "pt", "reference_date", "future_A_6m", 
                "future_N_6m", "future_prr_6m", "label_prr"}
feature_cols = [c for c in features.columns if c not in exclude_cols]

X_train = features.loc[train_mask, feature_cols]
y_train = features.loc[train_mask, "label_prr"]
X_test = features.loc[test_mask, feature_cols]
y_test = features.loc[test_mask, "label_prr"]

print(f"[INFO] Train: {X_train.shape}, Test: {X_test.shape}")
print(f"[INFO] Train positives: {y_train.sum()}, Test positives: {y_test.sum()}")

# 2. SAME MODEL PARAMS as pipeline
params = {
    "objective": "binary:logistic", "eval_metric": "aucpr",
    "max_depth": 6, "eta": 0.05, "subsample": 0.8,
    "colsample_bytree": 0.8, "lambda": 1.0, "alpha": 0.0,
    "n_estimators": 500, "tree_method": "hist", "base_score": 0.5,
}

print("[INFO] Training model for SHAP...")
model = xgb.XGBClassifier(**params)
model.fit(X_train, y_train)

print("[INFO] Model trained ✓")

# 3. Top 100 test predictions
y_pred_proba = model.predict_proba(X_test)[:, 1]
test_df = features.loc[test_mask, ["DRUGNAME_NORM", "pt", "label_prr"]].copy()
test_df["score"] = y_pred_proba
test_df["rank"] = test_df["score"].rank(ascending=False)

top100_idx = test_df.nsmallest(100, "rank").index.get_level_values(0)
X_top100 = X_test.loc[top100_idx]
print(f"[INFO] Top 100 shape: {X_top100.shape}")

# 4. SHAP
print("[INFO] Computing SHAP values...")
explainer = shap.TreeExplainer(model)
shap_values = explainer.shap_values(X_top100)

# 5. PLOTS
plt.style.use('default')

# Global summary
plt.figure(figsize=(12, 8))
shap.summary_plot(shap_values, X_top100, show=False)
plt.tight_layout()
plt.savefig(SHAP_DIR / "shap_summary_global.png", dpi=300, bbox_inches='tight')
plt.close()
print("✅ shap_summary_global.png")

# Feature importance bar
plt.figure(figsize=(10, 6))
shap.summary_plot(shap_values, X_top100, plot_type="bar", show=False)
plt.tight_layout()
plt.savefig(SHAP_DIR / "shap_feature_importance.png", dpi=300, bbox_inches='tight')
plt.close()
print("✅ shap_feature_importance.png")

# 6. TABLES
shap_importance = np.abs(shap_values).mean(axis=0)
importance_df = pd.DataFrame({
    "feature": feature_cols,
    "shap_importance": shap_importance
}).sort_values("shap_importance", ascending=False)

importance_df.to_csv(SHAP_DIR / "shap_feature_importance.csv", index=False)
top10 = test_df.nsmallest(10, "rank")[["DRUGNAME_NORM", "pt", "score", "label_prr"]].reset_index(drop=True)
top10["rank"] = range(1, 11)
top10.to_csv(SHAP_DIR / "top10_predictions.csv", index=False)

# 7. RESULTS
print("\n" + "="*60)
print("🏆 TOP 10 PREDICTIONS")
print("="*60)
print(top10.round(4).to_string(index=False))

print("\n" + "="*60)
print("📊 TOP 10 SHAP FEATURES")
print("="*60)
print(importance_df.head(10).round(4).to_string(index=False))

print("\n🎉 SHAP ANALYSIS COMPLETE!")
print(f"📁 Results: {SHAP_DIR}")
