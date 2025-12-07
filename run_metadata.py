from pathlib import Path
import json
from datetime import datetime

BASE_DIR = Path(r"C:\Users\koreo\Downloads\pv-signal-ml")
RESULTS_DIR = BASE_DIR / "results"
RESULTS_DIR.mkdir(exist_ok=True, parents=True)
RUN_LOG = RESULTS_DIR / "run_summaries.jsonl"


def log_run_metadata(
    model_version: str,
    data_period: str,
    label_definition: str,
    feature_list,
    train_size: int,
    test_size: int,
    train_pos: int,
    test_pos: int,
    ap_test,
    auc_test,
    mlflow_run_id: str,
):
    record = {
        "timestamp": datetime.utcnow().isoformat() + "Z",
        "model_version": model_version,
        "data_period": data_period,
        "label_definition": label_definition,
        "features": list(feature_list),
        "train_size": int(train_size),
        "test_size": int(test_size),
        "train_positives": int(train_pos),
        "test_positives": int(test_pos),
        "AP_test": None if ap_test != ap_test else float(ap_test),
        "AUC_test": None if auc_test != auc_test else float(auc_test),
        "mlflow_run_id": mlflow_run_id,
    }
    with RUN_LOG.open("a", encoding="utf-8") as f:
        f.write(json.dumps(record) + "\n")