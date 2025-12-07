# Setup and Data Guide - PV-Signal-ML

**Date:** 2025-12-08  
**Status:** 🔬 Research Prototype

---

## 📋 Overview

This guide explains how to set up PV-Signal-ML with all necessary data, models, and configurations. The application requires:

1. **FAERS Data** - FDA adverse event reports (downloaded on first run)
2. **MLflow Artifacts** - Trained ML models and experiment tracking
3. **Generated Reports** - Signal Assessment Reports (SARs) and PSMFs
4. **Database** - SQLite for data storage and audit logs

---

## 🚀 Quick Start (5 minutes)

### Step 1: Clone the Repository

```bash
git clone https://github.com/PrashantRGore/PV_Signal_ML.git
cd PV_Signal_ML
```

### Step 2: Create Virtual Environment

```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

### Step 3: Install Dependencies

```bash
pip install -r requirements.txt
```

### Step 4: Run Setup Script

```bash
python setup_data.py
```

This will:
- ✅ Create all necessary directories
- ✅ Initialize SQLite database
- ✅ Create sample FAERS data
- ✅ Set up MLflow tracking
- ✅ Configure Streamlit

### Step 5: Start the Application

```bash
streamlit run pv_fullstack.py
```

The app will open at `http://localhost:8501`

---

## 📊 Data Management

### FAERS Data

FAERS (FDA Adverse Event Reporting System) data is automatically downloaded when needed.

**Location:** `data/raw/faers/`

**Structure:**
```
data/raw/faers/
├── 2020_Q1/
│   ├── ASCII/
│   │   ├── DEMO20Q1.txt
│   │   ├── DRUG20Q1.txt
│   │   ├── REAC20Q1.txt
│   │   └── ...
│   └── DELETED/
└── 2020_Q2/
    └── ...
```

**Data Download:**
- Automatic: The app downloads FAERS data on first use
- Manual: Download from https://fis.fda.gov/sense/app/9524532e-2eb4-490e-b914-0a5f8e970e2d/sheet/7a5acf3b-72d4-4b5d-99a7-4ca3e9ca74d4/state/analysis

**Size:** ~500 MB for full dataset (2020-2025)

### Sample Data

For quick testing without downloading full FAERS:

```bash
python setup_data.py
```

This creates sample data in `data/raw/faers/` with 3 test cases.

---

## 🤖 MLflow Artifacts

MLflow tracks all model training and predictions.

**Location:** `mlruns/`

**Structure:**
```
mlruns/
├── 0/
│   ├── meta.yaml
│   └── artifacts/
│       ├── model.pkl
│       ├── metrics.json
│       └── params.json
└── 1/
    └── ...
```

**Auto-Generated:**
- Models are trained and logged automatically
- Metrics are tracked per experiment
- Artifacts are stored locally

**Access MLflow UI:**
```bash
mlflow ui
```

Then open `http://localhost:5000`

---

## 📄 Generated Reports

Reports are generated dynamically when you run signal detection.

**Location:** `sar_reports/reports/`

**Types:**
1. **Signal Assessment Reports (SARs)** - Per drug-reaction pair
2. **PSMFs** - Periodic Safety Update Format (EMA compliant)
3. **Audit Logs** - Complete audit trail

**Example:**
```
sar_reports/reports/
├── Aspirin__Headache__2024-01-01_2024-03-31.json
├── Ibuprofen__Nausea__2024-01-01_2024-03-31.json
└── ...
```

**Auto-Generated:**
- Reports are created when you run signal detection
- Cached for performance
- Regenerated if data changes

---

## 🗄️ Database

SQLite database stores all data and audit logs.

**Location:** `pv_signal.db`

**Tables:**
```sql
-- Individual Case Safety Reports
CREATE TABLE icsr (
    icsr_id TEXT PRIMARY KEY,
    report_id TEXT,
    case_id TEXT,
    drug_name TEXT,
    reaction TEXT,
    outcome TEXT,
    report_date TEXT,
    created_at TIMESTAMP
);

-- Detected Signals
CREATE TABLE signals (
    signal_id TEXT PRIMARY KEY,
    drug_name TEXT,
    reaction TEXT,
    prr REAL,
    chi_square REAL,
    expected_cases REAL,
    observed_cases REAL,
    confidence_interval TEXT,
    signal_strength TEXT,
    created_at TIMESTAMP
);

-- Audit Logs
CREATE TABLE audit_logs (
    log_id TEXT PRIMARY KEY,
    action TEXT,
    user_id TEXT,
    timestamp TIMESTAMP,
    details TEXT
);
```

**Auto-Created:**
- Database is initialized on first run
- Tables are created automatically
- Data is populated from FAERS

---

## 🔧 Configuration Files

### Streamlit Config (`.streamlit/config.toml`)

```toml
[theme]
primaryColor = "#1f77b4"
backgroundColor = "#ffffff"

[client]
showErrorDetails = true
```

### MLflow Config (`.mlflow_config.json`)

```json
{
  "tracking_uri": "file:///./mlruns",
  "experiment_name": "pv-signal-detection",
  "artifact_location": "./mlruns/artifacts"
}
```

### Environment Variables

Create `.env` file (optional):

```bash
# Data paths
FAERS_DATA_PATH=data/raw/faers
MLFLOW_TRACKING_URI=file:///./mlruns

# API settings
API_HOST=127.0.0.1
API_PORT=8000

# Streamlit settings
STREAMLIT_SERVER_PORT=8501
```

---

## 📦 Directory Structure

```
pv-signal-ml/
├── .streamlit/
│   └── config.toml
├── data/
│   ├── raw/
│   │   └── faers/          # FAERS data (auto-downloaded)
│   └── processed/          # Processed data
├── mlruns/                 # MLflow experiments & models
├── sar_reports/
│   └── reports/            # Generated SARs & PSMFs
├── audit_logs/             # Audit trail
├── gdpr_registry/          # GDPR deletion registry
├── lineage/                # Data lineage tracking
├── src/                    # Source code
├── templates/              # Report templates
├── pv_signal.db            # SQLite database
├── pv_fullstack.py         # Streamlit app
├── api.py                  # FastAPI backend
├── setup_data.py           # Setup script
├── requirements.txt        # Dependencies
└── README.md               # Documentation
```

---

## 🔄 Data Flow

```
FAERS Data (FDA)
    ↓
setup_data.py (Download & Parse)
    ↓
SQLite Database (pv_signal.db)
    ↓
Signal Detection (PRR, Chi-square)
    ↓
ML Triage (XGBoost)
    ↓
MLflow Tracking (mlruns/)
    ↓
Report Generation (SARs, PSMFs)
    ↓
Streamlit UI / FastAPI
```

---

## 🚨 Troubleshooting

### Issue: "FAERS data not found"

**Solution:**
```bash
python setup_data.py
```

This downloads sample data. For full data, download from FDA website.

### Issue: "MLflow experiments not showing"

**Solution:**
```bash
mlflow ui
```

Then open `http://localhost:5000` to view experiments.

### Issue: "Database locked"

**Solution:**
```bash
rm pv_signal.db
python setup_data.py
```

This recreates the database.

### Issue: "Reports not generating"

**Solution:**
1. Ensure `sar_reports/reports/` directory exists
2. Check that Ollama is running (for RAG)
3. Check logs in `audit_logs/`

---

## 📥 Downloading Full FAERS Data

For production use with real data:

1. Go to https://fis.fda.gov/sense/app/9524532e-2eb4-490e-b914-0a5f8e970e2d/sheet/7a5acf3b-72d4-4b5d-99a7-4ca3e9ca74d4/state/analysis
2. Download quarterly FAERS files (2020-2025)
3. Extract to `data/raw/faers/`
4. Run `python setup_data.py` to process

**Note:** Full dataset is ~500 MB. Sample data (3 cases) is included for quick testing.

---

## 🔐 Security & Privacy

### Data Protection
- ✅ Aggregated data only (no PII)
- ✅ GDPR deletion registry implemented
- ✅ Audit logs for all operations
- ✅ Data lineage tracking

### Credentials
- ❌ No API keys in code
- ❌ No passwords in config files
- ✅ Environment variables for sensitive data

---

## 📚 Additional Resources

- **README.md** - Project overview
- **VALIDATION_STATUS.md** - Validation & compliance status
- **SECURITY_REVIEW.md** - Security assessment
- **governance_dpia.md** - Data protection impact assessment

---

## ✅ Verification Checklist

After setup, verify everything works:

```bash
# 1. Check directories created
ls -la data/ mlruns/ sar_reports/ audit_logs/

# 2. Check database initialized
sqlite3 pv_signal.db ".tables"

# 3. Check sample data exists
ls data/raw/faers/

# 4. Check dependencies installed
pip list | grep -E "streamlit|fastapi|xgboost|mlflow"

# 5. Start the app
streamlit run pv_fullstack.py

# 6. Open browser to http://localhost:8501
```

---

## 🎯 Next Steps

1. **Run Setup:** `python setup_data.py`
2. **Start App:** `streamlit run pv_fullstack.py`
3. **Explore Data:** Use the Streamlit UI to view signals
4. **Check Reports:** View generated SARs in `sar_reports/reports/`
5. **Review Logs:** Check audit logs in `audit_logs/`

---

## 📞 Support

For issues or questions:

1. Check this guide
2. Review README.md
3. Check VALIDATION_STATUS.md
4. Open an issue on GitHub

---

**Last Updated:** 2025-12-08  
**Status:** 🔬 Research Prototype

---

*This guide ensures all necessary data, models, and configurations are properly set up for the PV-Signal-ML application.*
