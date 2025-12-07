from pathlib import Path
import json

import mlflow
import pandas as pd
from audit_logging import audit_logger
# import xgboost as xgb  # uncomment and use when you plug in real training

# ---------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------
BASE_DIR = Path(r"C:\Users\koreo\Downloads\pv-signal-ml")
SAR_DIR = BASE_DIR / "sar_reports"

# -----------------------------------------------------------------------
# MLflow: unified file-based tracking (portable, no external dependencies)
# -----------------------------------------------------------------------
MLFLOW_URI = "file:///C:/Users/koreo/mlruns"
mlflow.set_tracking_uri(MLFLOW_URI)
mlflow.set_experiment("pv-signal-ml")


def train_model():
    """
    Minimal training wrapper that creates an MLflow run and writes
    model_run_metadata.json used by the signal reports.

    Replace the commented section with your real XGBoost training.
    """
    with mlflow.start_run() as run:
        run_id = run.info.run_id
        model_version = "1.0"   # bump when you change the model
        data_period = "FAERS-ALL"

        # ---- log high-level metadata ----
        mlflow.log_param("model_type", "xgboost")
        mlflow.log_param("data_period", data_period)
        mlflow.log_param("dp_enabled", False)
        mlflow.log_param("epsilon", None)

        # ---- your real training code goes here (optional placeholder) ----
        # example:
        # data = pd.read_parquet(SAR_DIR / "training_data.parquet")
        # X = data.drop(columns=["label"])
        # y = data["label"]
        # model = xgb.XGBClassifier(...)
        # model.fit(X, y)
        # mlflow.xgboost.log_model(model, "model")

        # ---- write lightweight metadata for report builder ----
        meta = {
            "run_id": run_id,
            "model_version": model_version,
            "tracking_uri": MLFLOW_URI,
            "experiment_name": "pv-signal-ml",
            "data_period": data_period,
            "dp_enabled": False,
            "epsilon": None,
        }
        out = SAR_DIR / "model_run_metadata.json"
        out.parent.mkdir(exist_ok=True, parents=True)
        with open(out, "w", encoding="utf-8") as f:
            json.dump(meta, f, indent=2)

        # simple methods artifact for audit trails
        methods = {
            "input_signals_file": "full_signals_1M.csv",
            "faers_processing": "faers_build_signals.py exact PRR",
            "prr_threshold": 2.0,
            "chisq_threshold": 4.0,
        }
        methods_path = SAR_DIR / f"methods_{run_id}.json"
        with open(methods_path, "w", encoding="utf-8") as f:
            json.dump(methods, f, indent=2)

        # ---- Log MLflow run to audit trail (FDA 21 CFR Part 11 compliance) ----
        audit_logger.log_mlflow_run(
            run_id=run_id,
            experiment_name="pv-signal-ml",
            model_version=model_version,
            metrics={},
            parameters={
                "model_type": "xgboost",
                "data_period": data_period,
                "dp_enabled": False,
            },
            user_id="system",
            status="completed"
        )

        print("✅ MLflow run started and metadata saved:", run_id)
        print(f"✅ MLflow run {run_id} logged to audit trail")


if __name__ == "__main__":
    train_model()
