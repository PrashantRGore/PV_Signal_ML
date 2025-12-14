"""
Central Configuration for PV_Signal_ML Enhanced System
Implements SISA, DP-SGD, and SHAP configurations
"""
import os
from pathlib import Path

# ===== PROJECT PATHS =====
PROJECT_ROOT = Path(__file__).parent
DATA_DIR = PROJECT_ROOT / "data"
MODELS_DIR = PROJECT_ROOT / "models"
RESULTS_DIR = PROJECT_ROOT / "results"
LOGS_DIR = PROJECT_ROOT / "logs"

# Create directories if they don't exist
for dir_path in [DATA_DIR, MODELS_DIR, RESULTS_DIR, LOGS_DIR]:
    dir_path.mkdir(parents=True, exist_ok=True)

# ===== DATA SOURCE CONFIGURATION =====
DATA_SOURCE_OPTIONS = ["local", "demo_hf", "faers_live"]
DEFAULT_DATA_SOURCE = "demo_hf"

# HuggingFace Demo Dataset
HF_REPO_ID = "PrashantRGore/pv-signal-ml-data"
HF_CACHE_DIR = DATA_DIR / "cache" / "huggingface"
HF_FILES = ["demo.parquet", "drug.parquet", "reac.parquet", "outc.parquet"]

# FAERS Configuration
FAERS_BASE_URL = "https://fis.fda.gov/content/Exports/"
FAERS_RAW_DIR = DATA_DIR / "faers" / "raw"
FAERS_PROCESSED_DIR = DATA_DIR / "faers" / "processed"
FAERS_QUARTERS = {
    1: (1, 3),   # Q1: Jan-Mar
    2: (4, 6),   # Q2: Apr-Jun
    3: (7, 9),   # Q3: Jul-Sep
    4: (10, 12)  # Q4: Oct-Dec
}

# ===== SISA ARCHITECTURE CONFIGURATION =====
ENABLE_SISA = True  # Toggle SISA machine unlearning
NUM_SHARDS = 10  # Number of data shards for SISA
RECORDS_PER_SLICE = 10000  # Records per checkpoint slice
SISA_CHECKPOINT_DIR = MODELS_DIR / "sisa_checkpoints"
SISA_MAPPING_FILE = MODELS_DIR / "sisa_shard_mapping.json"
SISA_ENSEMBLE_METHOD = "average"  # Options: "average", "vote", "weighted"

# ===== DP-SGD CONFIGURATION =====
ENABLE_DP_SGD = False  # Toggle differential privacy (EXPERIMENTAL)
DP_EPSILON = 3.0  # Privacy budget (lower = more private, less accurate)
DP_DELTA = 1e-5  # Probability of privacy breach
DP_CLIP_NORM = 1.0  # Gradient clipping threshold
DP_NOISE_MULTIPLIER = 1.1  # Noise scale: σ = noise_multiplier * clip_norm

# ===== ML ENGINE CONFIGURATION =====
ML_MODEL_TYPE = "xgboost"
XGBOOST_PARAMS = {
    "max_depth": 6,
    "learning_rate": 0.1,
    "n_estimators": 100,
    "objective": "binary:logistic",
    "eval_metric": "auc",
    "random_state": 42,
    "n_jobs": -1
}

# MLflow Configuration
MLFLOW_TRACKING_URI = "file:///" + str(LOGS_DIR / "mlruns")
MLFLOW_EXPERIMENT_NAME = "PV_Signal_Detection"

# ===== SHAP EXPLAINABILITY CONFIGURATION =====
ENABLE_SHAP = True
SHAP_TOP_N_SIGNALS = 50  # Generate force plots for top N signals
SHAP_OUTPUT_DIR = RESULTS_DIR / "shap"
SHAP_SAMPLE_SIZE = 500  # Max samples for SHAP computation (performance)

# ===== STATISTICS ENGINE CONFIGURATION =====
STATS_MIN_CASE_COUNT = 3  # Minimum cases for signal detection
STATS_PRR_THRESHOLD = 2.0  # Proportional Reporting Ratio threshold
STATS_CHI2_THRESHOLD = 3.841  # Chi-square threshold (p < 0.05, df=1)
STATS_CONFIDENCE_LEVEL = 0.95

# ===== RAG CONFIGURATION =====
RAG_ENABLE = True
RAG_MODEL = "llama3.2"  # Ollama model
RAG_TEMPERATURE = 0.1  # Low temperature for deterministic outputs
RAG_MAX_TOKENS = 1000
RAG_CHROMA_PERSIST_DIR = DATA_DIR / "chroma_db"

# ===== REGULATORY COMPLIANCE =====
REGULATORY_FRAMEWORKS = ["EMA_GVP_IX", "FDA_21CFR11", "GDPR", "HIPAA"]
ENABLE_AUDIT_LOGGING = True
AUDIT_LOG_FILE = LOGS_DIR / "audit.log"

# ===== DISCLAIMERS =====
DISCLAIMER_TEXT = """
⚠️ RESEARCH/EDUCATIONAL PROTOTYPE
- Not for direct clinical or regulatory decisions
- Data must be de-identified (no PII/PHI)
- FAERS and HuggingFace data subject to respective licenses
- SISA and DP-SGD are experimental privacy-preserving techniques
- For production use, validation by QPPV and DPO required
"""

# ===== PERFORMANCE THRESHOLDS =====
MAX_MEMORY_GB = 8  # System RAM constraint
DASHBOARD_REFRESH_MAX_SECONDS = 10
SIGNAL_DETECTION_SENSITIVITY_THRESHOLD = 0.85  # Minimum required

# ===== VALIDATION STATUS =====
VALIDATION_STATUS = {
    "code_version": "v2.0-enhanced",
    "last_validation_date": None,
    "validated_by": None,
    "gxp_compliant": False,
    "production_ready": False
}

# Export all configurations
__all__ = [
    "PROJECT_ROOT", "DATA_DIR", "MODELS_DIR", "RESULTS_DIR", "LOGS_DIR",
    "DATA_SOURCE_OPTIONS", "HF_REPO_ID", "FAERS_BASE_URL",
    "ENABLE_SISA", "NUM_SHARDS", "ENABLE_DP_SGD", "DP_EPSILON",
    "ENABLE_SHAP", "MLFLOW_TRACKING_URI", "DISCLAIMER_TEXT",
    "REGULATORY_FRAMEWORKS", "VALIDATION_STATUS"
]
