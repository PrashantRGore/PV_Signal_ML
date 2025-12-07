# Critical Fixes Implementation Guide

**Date:** 2025-12-08  
**Status:** IN PROGRESS - Phase 1 (MLflow Audit Integration)

---

## ✅ COMPLETED FIXES

### Fix 1: pv_signal_ml_pipeline.py ✅
**What:** Added MLflow audit logging to XGBoost training pipeline
**Changes:**
- Added `from audit_logging import audit_logger`
- Added `audit_logger.log_mlflow_run()` call after MLflow logging
- Logs: run_id, experiment, metrics (AP, AUC), parameters
- **Status:** ✅ COMPLETE

**Code Added:**
```python
# 6. Log MLflow run to audit trail (FDA 21 CFR Part 11 compliance)
audit_logger.log_mlflow_run(
    run_id=run_id,
    experiment_name=MLFLOW_EXPERIMENT,
    model_version="xgboost_prr_v1",
    metrics={...},
    parameters=params,
    user_id="system",
    status="completed"
)
```

---

### Fix 2: train_with_mlflow.py ✅
**What:** Added MLflow audit logging to training wrapper
**Changes:**
- Added `from audit_logging import audit_logger`
- Changed tracking URI to unified file-based: `file:///C:/Users/koreo/mlruns`
- Added `audit_logger.log_mlflow_run()` call
- **Status:** ✅ COMPLETE

**Code Added:**
```python
# Unified MLflow tracking URI
MLFLOW_URI = "file:///C:/Users/koreo/mlruns"
mlflow.set_tracking_uri(MLFLOW_URI)

# Log to audit trail
audit_logger.log_mlflow_run(
    run_id=run_id,
    experiment_name="pv-signal-ml",
    model_version=model_version,
    metrics={},
    parameters={...},
    user_id="system",
    status="completed"
)
```

---

## 🔄 IN PROGRESS FIXES

### Fix 3: api.py (NEXT)
**What:** Add MLflow and audit logging to `/signal-report` endpoint
**Location:** Lines 44-127
**Changes Needed:**
1. Import MLflow and audit_logger
2. Log report generation to MLflow
3. Log to audit trail
4. Track request parameters and response metrics

**Code to Add:**
```python
import mlflow
from audit_logging import audit_logger

@app.post("/signal-report")
def generate_signal_report(request: SignalRequest):
    try:
        # Start MLflow run for this report generation
        mlflow.set_tracking_uri("file:///C:/Users/koreo/mlruns")
        mlflow.set_experiment("pv-signal-ml-reports")
        
        with mlflow.start_run():
            # ... existing report generation code ...
            
            # Log to MLflow
            mlflow.log_param("drug", request.drug)
            mlflow.log_param("event", request.event)
            mlflow.log_param("period", request.period)
            mlflow.log_metric("report_size", len(json.dumps(sar)))
            
            # Log to audit trail
            audit_logger.log_report_generation(
                report_type="SAR",
                drug=request.drug,
                event=request.event,
                period=request.period,
                user_id="api",
                output_file=str(report_path),
                success=True
            )
            
            return {...}
    except Exception as e:
        # Log failure
        audit_logger.log_report_generation(
            report_type="SAR",
            drug=request.drug,
            event=request.event,
            period=request.period,
            user_id="api",
            success=False,
            error_message=str(e)
        )
        return {...}
```

---

### Fix 4: rag_langchain.py (NEXT)
**What:** Add MLflow logging to RAG pipeline
**Location:** Lines 35-69
**Changes Needed:**
1. Import MLflow
2. Log RAG operations
3. Track embedding model, LLM parameters
4. Track similarity search results

**Code to Add:**
```python
import mlflow

class PVSignalRAGLangChain:
    def __init__(self):
        # ... existing code ...
        mlflow.set_tracking_uri("file:///C:/Users/koreo/mlruns")
        mlflow.set_experiment("pv-signal-ml-rag")
    
    def generate_sar_report(self, signal_row):
        with mlflow.start_run():
            # Log RAG parameters
            mlflow.log_param("embedding_model", "all-MiniLM-L6-v2")
            mlflow.log_param("llm_model", "llama3.2:3b")
            mlflow.log_param("temperature", 0.1)
            mlflow.log_param("drug", signal_row['DRUG'])
            mlflow.log_param("event", signal_row['EVENT'])
            
            # ... existing RAG code ...
            
            # Log metrics
            mlflow.log_metric("prr", signal_row['PRR'])
            mlflow.log_metric("chi2", signal_row['CHISQ'])
            mlflow.log_metric("cases", signal_row['CASES'])
            mlflow.log_metric("response_length", len(response))
            
            return {...}
```

---

### Fix 5: signal_report_builder.py (NEXT)
**What:** Add MLflow logging to SAR builder
**Location:** Lines 66-205
**Changes Needed:**
1. Import MLflow
2. Log report generation
3. Track input parameters and output metrics

**Code to Add:**
```python
import mlflow

def build_signal_report(drug: str, event: str, period: str = "FAERS-ALL"):
    mlflow.set_tracking_uri("file:///C:/Users/koreo/mlruns")
    mlflow.set_experiment("pv-signal-ml-reports")
    
    with mlflow.start_run():
        # Log input parameters
        mlflow.log_param("drug", drug)
        mlflow.log_param("event", event)
        mlflow.log_param("period", period)
        
        # ... existing report generation code ...
        
        # Log output metrics
        mlflow.log_metric("cases", statistics['cases'])
        mlflow.log_metric("prr", statistics['prr'])
        mlflow.log_metric("chisq", statistics['chisq'])
        mlflow.log_metric("evidence_items", len(evidence_items))
        mlflow.log_artifact(str(json_path))
        mlflow.log_artifact(str(md_path))
        
        return json_path, md_path
```

---

## 🎯 ENHANCEMENT OPPORTUNITIES

### Enhancement 1: Real-Time Monitoring Dashboard
**File:** `signal_monitoring_dashboard.py` (NEW)
**Features:**
- Live signal trends
- New signal alerts
- Signal status changes
- Performance metrics

### Enhancement 2: Signal Causality Scoring
**File:** `causality_scorer.py` (NEW)
**Features:**
- WHO-UMC causality assessment
- Automated scoring
- Evidence-based recommendations

### Enhancement 3: Benefit-Risk Visualization
**File:** `benefit_risk_analyzer.py` (NEW)
**Features:**
- Interactive charts
- Risk-benefit tradeoffs
- Regulatory decision support

### Enhancement 4: Regulatory Submission Generator
**File:** `regulatory_submission_generator.py` (NEW)
**Features:**
- Auto-generate EMA/FDA documents
- Template-based generation
- Compliance checking

### Enhancement 5: Signal Comparison
**File:** `signal_comparison.py` (NEW)
**Features:**
- Time-series analysis
- Trend detection
- Period comparison

### Enhancement 6: Automated Alerts
**File:** `alert_system.py` (NEW)
**Features:**
- Email notifications
- Slack integration
- Configurable thresholds

### Enhancement 7: Advanced Search
**File:** `advanced_search.py` (NEW)
**Features:**
- Full-text search
- Faceted filtering
- Saved searches

---

## 📊 IMPLEMENTATION TIMELINE

### Phase 1: Critical Fixes (2-3 hours)
- [x] Fix pv_signal_ml_pipeline.py
- [x] Fix train_with_mlflow.py
- [ ] Fix api.py
- [ ] Fix rag_langchain.py
- [ ] Fix signal_report_builder.py

### Phase 2: RAG Enhancements (2-3 hours)
- [ ] Proper LangChain chains
- [ ] Prompt templates
- [ ] Memory management
- [ ] Better RAG quality

### Phase 3: Feature Enhancements (4-6 hours)
- [ ] Real-time dashboard
- [ ] Causality scoring
- [ ] Benefit-risk analysis
- [ ] Regulatory submission generator
- [ ] Signal comparison
- [ ] Automated alerts
- [ ] Advanced search

### Phase 4: Testing & Deployment (2-4 hours)
- [ ] Comprehensive testing
- [ ] GitHub deployment
- [ ] Streamlit Cloud deployment
- [ ] Documentation

---

## 🔍 VERIFICATION CHECKLIST

### MLflow Audit Integration
- [ ] All MLflow runs logged to audit trail
- [ ] Audit logs created in `audit_logs/mlflow_runs.jsonl`
- [ ] All metrics captured
- [ ] All parameters captured
- [ ] Timestamps correct
- [ ] User IDs tracked

### RAG Pipeline
- [ ] RAG operations logged to MLflow
- [ ] Embedding model tracked
- [ ] LLM parameters tracked
- [ ] Response quality metrics logged
- [ ] Similarity search results tracked

### Report Generation
- [ ] All reports logged to MLflow
- [ ] All reports logged to audit trail
- [ ] Input parameters captured
- [ ] Output metrics captured
- [ ] File paths tracked
- [ ] Success/failure status logged

### Compliance
- [ ] FDA 21 CFR Part 11 audit trail complete
- [ ] EMA GVP Module IX requirements met
- [ ] CIOMS XIV requirements met
- [ ] GDPR compliance maintained
- [ ] ICH E2A requirements met

---

## 🚀 DEPLOYMENT READINESS

**Current Status:** 🟡 PHASE 1 IN PROGRESS

**Blockers:** None (can proceed with remaining fixes)

**Next Actions:**
1. Complete api.py fix
2. Complete rag_langchain.py fix
3. Complete signal_report_builder.py fix
4. Test all MLflow logging
5. Proceed to Phase 2 enhancements

---

**Estimated Time to Production:** 8-12 hours  
**Risk Level:** LOW (core fixes are straightforward)

---

*This guide ensures systematic implementation of critical fixes and enhancements for production deployment.*
