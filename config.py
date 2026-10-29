"""
Configuration settings for PV Signal ML System
Research & Educational Prototype for Pharmacovigilance Signal Detection
Version: 3.0.0-prototype
"""

from pathlib import Path

# ===== PROJECT PATHS =====
BASE_DIR = Path(__file__).parent
DATA_DIR = BASE_DIR / "data"
MODELS_DIR = BASE_DIR / "models"
LOGS_DIR = BASE_DIR / "logs"
AUDIT_DIR = LOGS_DIR / "audit"
MLRUNS_DIR = BASE_DIR / "mlruns"

# Create directories if they don't exist
for directory in [DATA_DIR, MODELS_DIR, LOGS_DIR, AUDIT_DIR, MLRUNS_DIR]:
    directory.mkdir(exist_ok=True, parents=True)

# ===== AUDIT & PROTOTYPE LOGGING =====
ENABLE_AUDIT_LOGGING = True  # Audit trail simulation (demonstrates 21 CFR Part 11 concepts)
AUDIT_LOG_FILE = AUDIT_DIR / "pv_signal_audit.log"
AUDIT_LOG_RETENTION_DAYS = 2555  # 7 years (standard PV reference period)

# ===== DATA SOURCE CONFIGURATION =====
DATA_SOURCE_OPTIONS = ["local", "demo_hf", "faers_live"]

# HuggingFace Demo Dataset (Synthetic ICSRs)
HF_REPO_ID = "PrashantRGore/synthetic-faers-1m-v3"
HF_DATA_FILE = "synthetic_faers_1m_v3.csv"
HF_CACHE_DIR = str(DATA_DIR / "hf_cache")

# FAERS Live / Public Data
FAERS_BASE_URL = "https://fis.fda.gov/content/Exports"
FAERS_RAW_DIR = DATA_DIR / "faers_raw"
FAERS_PROCESSED_DIR = DATA_DIR / "faers_processed"

# ===== ML CONFIGURATION =====
NUM_SHARDS = 10  # For SISA architecture
RANDOM_SEED = 42
TEST_SIZE = 0.2

# Model Versioning & Directory
MODEL_VERSION = "3.0.0-prototype"
MODEL_DIR = DATA_DIR / "models"
MODEL_DIR.mkdir(parents=True, exist_ok=True)
MODEL_REGISTRY_DIR = MODELS_DIR / "registry"

# ===== MLFLOW CONFIGURATION =====
_mlruns_path = str(MLRUNS_DIR.absolute()).replace("\\", "/")
MLFLOW_TRACKING_URI = f"file:///{_mlruns_path}"
MLFLOW_EXPERIMENT_NAME = "PV_Signal_Detection"
MLFLOW_ENABLE_SYSTEM_METRICS = False
MLFLOW_ARTIFACT_LOCATION = str(MLRUNS_DIR / "artifacts")

# ===== RAG / LLM CONFIGURATION =====
OLLAMA_BASE_URL = "http://localhost:11434"
DEFAULT_LLM_MODEL = "llama3.2"  # Governed by Meta Llama Community License Agreement
LLM_TEMPERATURE = 0.7
LLM_MAX_TOKENS = 2000

# ===== LOGGING CONFIGURATION =====
LOG_LEVEL = "INFO"
LOG_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
LOG_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"
ENABLE_CONSOLE_LOGGING = True
ENABLE_FILE_LOGGING = True
LOG_FILE = LOGS_DIR / "pv_signal.log"
LOG_FILE_MAX_BYTES = 10 * 1024 * 1024  # 10 MB
LOG_FILE_BACKUP_COUNT = 5

# ===== CANONICAL SCHEMA CONFIGURATION =====
CANONICAL_SCHEMA_VERSION = "1.0"
TIER1_REQUIRED_COLUMNS = [
    "drug_name",
    "event_term",
    "case_count",
    "prr",
    "ror",
    "chi2",
    "label",
]

# ===== SIGNAL DETECTION THRESHOLDS =====
MIN_CASE_COUNT = 3
PRR_THRESHOLD = 2.0
CHI2_THRESHOLD = 4.0

# ===== SISA MACHINE UNLEARNING CONFIGURATION =====
# SISA = Sharded, Isolated, Sliced, Aggregated (Bourtoule et al., 2021)
# Note: In regulatory pharmacovigilance, statutory public health monitoring obligations
# (e.g. GDPR Art. 17(3)(c) and Art. 9(2)(i)) generally provide legal exceptions to erasure.
# SISA is implemented here as an algorithmic proof-of-concept for research and education.
ENABLE_RIGHT_TO_BE_FORGOTTEN = True
UNLEARNING_METHOD = "SISA"

# Audit trail event types
AUDIT_TRAIL_EVENTS = [
    "data_load",
    "signal_detection",
    "model_training",
    "model_prediction",
    "unlearning_request",
    "sar_generation",
    "data_export",
    "user_action",
]

# ===== PERFORMANCE CONFIGURATION =====
MAX_WORKERS = 4
CHUNK_SIZE = 10000
CACHE_ENABLED = True
CACHE_DIR = DATA_DIR / "cache"

# ===== REGULATORY DISCLAIMER =====
DISCLAIMER_TEXT = """
⚠️ **RESEARCH & EDUCATIONAL PROTOTYPE ONLY**

This platform is a proof-of-concept research implementation demonstrating pharmacovigilance
signal detection, statistical disproportionality, SISA machine unlearning, and explainability.

**Not a Production System:**
- NOT FDA 21 CFR Part 11 validated (no GAMP-5 computer system validation performed).
- NOT approved for real patient PHI/PII.
- Not a substitute for commercial pharmacovigilance safety databases (Argus, ArisGlobal, SAS).
- SISA machine unlearning is an algorithmic exploration; statutory PV safety mandates
  typically take precedence under GDPR Art. 17(3)(c). Consult a Data Protection Officer.
- All algorithmic signals require clinical review by qualified medical reviewers.
"""

# ===== FEATURE FLAGS =====
FEATURES = {
    "canonical_schema": True,
    "tier2_features": True,
    "audit_logging": ENABLE_AUDIT_LOGGING,
    "mlflow_tracking": True,
    "shap_explainability": True,
    "rag_sar_generation": True,
    "sisa_unlearning": ENABLE_RIGHT_TO_BE_FORGOTTEN,
    "drug_portfolio_filter": True,
}

# ===== SYSTEM VALIDATION STATUS =====
VALIDATION_STATUS = {
    "code_version": "3.0.0-prototype",
    "last_validation": "2026-01-11",
    "validation_level": "None (Research / Educational Prototype)",
    "deployment_status": "Research Prototype (Not Validated for Regulatory Submissions)",
    "regulatory_framework_alignment": ["EMA GVP Module IX (Educational)", "CIOMS XIV (Conceptual)"],
}

# ===== TIER ANALYSIS MODES =====
ANALYSIS_MODES = {
    "tier1": {
        "name": "Quantitative (Tier 1)",
        "features": ["case_count", "prr", "ror", "chi2"],
        "description": "Baseline disproportionality metrics for signal screening",
    },
    "tier2": {
        "name": "Qualitative (Tier 2)",
        "features": "all",
        "description": "Full engineered clinical & causal feature set",
    },
}

# ===== SERVICE HEALTH CHECK ENDPOINTS =====
SERVICE_ENDPOINTS = {
    "ollama": OLLAMA_BASE_URL,
    "mlflow": MLFLOW_TRACKING_URI,
}
