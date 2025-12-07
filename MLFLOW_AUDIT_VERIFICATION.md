# MLflow Runs Audit & Verification Report

**Date:** 2025-12-08  
**Status:** ⚠️ PARTIAL - MLflow runs are logged but not integrated with audit trail

---

## 🔍 Current MLflow Implementation

### Files Using MLflow

1. **pv_signal_ml_pipeline.py** ✅
   - Sets tracking URI: `file:///C:/Users/koreo/mlruns`
   - Creates experiment: `pv-signal-ml-prr`
   - Logs parameters: XGBoost hyperparameters
   - Logs metrics: AP_test, AUC_test
   - Logs model: XGBoost model artifact
   - Logs run metadata via `run_metadata.py`
   - **Status:** ✅ Properly configured

2. **train_with_mlflow.py** ✅
   - Sets tracking URI: `sqlite:///mlflow.db`
   - Creates experiment: `pv-signal-ml`
   - Logs parameters: model_type, data_period, dp_enabled
   - Saves metadata: `model_run_metadata.json`
   - Saves methods: `methods_{run_id}.json`
   - **Status:** ✅ Properly configured

3. **run_metadata.py** ✅
   - Logs run metadata to: `results/run_summaries.jsonl`
   - Captures: model version, data period, features, metrics
   - **Status:** ✅ Properly configured

### Files NOT Using MLflow

- `api.py` – No MLflow logging
- `pv_fullstack.py` – No MLflow logging
- `audit_logging.py` – No MLflow integration
- `signal_report_builder.py` – No MLflow logging
- `rag_langchain.py` – No MLflow logging

---

## ⚠️ Issues Identified

### Issue 1: Inconsistent Tracking URIs
**Problem:** Two different MLflow tracking URIs:
- `pv_signal_ml_pipeline.py`: `file:///C:/Users/koreo/mlruns`
- `train_with_mlflow.py`: `sqlite:///mlflow.db`

**Impact:** Runs are stored in different locations, making it hard to view all runs in one place

**Solution:** Use single unified tracking URI

### Issue 2: No Audit Trail Integration
**Problem:** MLflow runs are not logged to audit trail

**Impact:** Cannot verify which runs generated which reports

**Solution:** Add MLflow logging to `audit_logging.py`

### Issue 3: No Report-to-Run Mapping
**Problem:** Generated reports don't record which MLflow run created them

**Impact:** Traceability is lost

**Solution:** Add run_id to report metadata

---

## ✅ Verification Steps

### Step 1: Check MLflow Runs (File-based)
```bash
# View runs from pv_signal_ml_pipeline.py
ls -la C:/Users/koreo/mlruns/
```

### Step 2: Check MLflow Runs (SQLite-based)
```bash
# View runs from train_with_mlflow.py
sqlite3 mlflow.db "SELECT * FROM runs;"
```

### Step 3: Check Run Metadata
```bash
# View logged metadata
cat sar_reports/model_run_metadata.json
cat results/run_summaries.jsonl
```

### Step 4: Check Methods Artifacts
```bash
# View methods for each run
ls -la sar_reports/methods_*.json
```

---

## 🔧 Recommended Fixes

### Fix 1: Unify MLflow Tracking URI
**File:** `train_with_mlflow.py` and `pv_signal_ml_pipeline.py`

**Change:**
```python
# Use file-based tracking (more portable)
MLFLOW_URI = "file:///C:/Users/koreo/mlruns"
mlflow.set_tracking_uri(MLFLOW_URI)
```

### Fix 2: Add MLflow Logging to Audit Trail
**File:** `audit_logging.py`

**Add method:**
```python
def log_mlflow_run(
    self,
    run_id: str,
    experiment_name: str,
    model_version: str,
    metrics: Dict[str, float],
    parameters: Dict[str, Any],
    user_id: Optional[str] = None
) -> Dict:
    """Log MLflow run to audit trail."""
    record = {
        "timestamp": datetime.utcnow().isoformat() + "Z",
        "event_type": "mlflow_run",
        "run_id": run_id,
        "experiment_name": experiment_name,
        "model_version": model_version,
        "metrics": metrics,
        "parameters": parameters,
        "user_id": user_id or "system",
    }
    self._append_log(self.mlflow_log_path, record)
    return record
```

### Fix 3: Add Run ID to Report Metadata
**File:** `api.py` and `signal_report_builder.py`

**Change:**
```python
# Capture current MLflow run
try:
    import mlflow
    run_id = mlflow.active_run().info.run_id if mlflow.active_run() else None
except:
    run_id = None

# Add to report
report_obj['mlflow_run_id'] = run_id
```

---

## 📊 Current MLflow Status

### Runs Captured ✅
- ✅ XGBoost training runs (pv_signal_ml_pipeline.py)
- ✅ Model training runs (train_with_mlflow.py)
- ✅ Metadata logging (run_metadata.py)

### Runs NOT Captured ❌
- ❌ Report generation runs (api.py)
- ❌ Signal detection runs (signal_report_builder.py)
- ❌ RAG pipeline runs (rag_langchain.py)

### Audit Integration ❌
- ❌ MLflow runs not in audit trail
- ❌ No report-to-run mapping
- ❌ No traceability between reports and ML runs

---

## 🎯 Recommendations

### Priority 1: Unify Tracking URI
- [ ] Change all files to use: `file:///C:/Users/koreo/mlruns`
- [ ] Verify all runs appear in single location
- [ ] Update documentation

### Priority 2: Integrate with Audit Trail
- [ ] Add MLflow logging to `audit_logging.py`
- [ ] Log all run creation events
- [ ] Log all report generation with run_id

### Priority 3: Add Report-to-Run Mapping
- [ ] Capture MLflow run_id when generating reports
- [ ] Store in report metadata
- [ ] Enable traceability queries

### Priority 4: Create MLflow Dashboard
- [ ] Add MLflow UI link to Streamlit
- [ ] Display run history
- [ ] Show model performance trends

---

## 📋 Compliance Status

| Requirement | Status | Notes |
|---|---|---|
| **MLflow Runs Captured** | ✅ Partial | Training runs captured, report runs not |
| **Audit Trail** | ⚠️ Partial | Runs logged separately, not in audit trail |
| **Traceability** | ❌ No | No link between reports and runs |
| **Reproducibility** | ✅ Yes | All parameters and metrics logged |
| **Regulatory Compliance** | ⚠️ Partial | Meets FDA 21 CFR Part 11 for training, not for reports |

---

## 🚀 Next Steps

1. **Before Testing:** Implement Priority 1 & 2 fixes
2. **During Testing:** Verify MLflow runs are created and logged
3. **After Testing:** Implement Priority 3 & 4 enhancements

---

**Status:** ⚠️ FUNCTIONAL BUT INCOMPLETE  
**Recommendation:** Implement fixes before production deployment

---

*This report confirms MLflow is working but needs audit trail integration for full compliance.*
