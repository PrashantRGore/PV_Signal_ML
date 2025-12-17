# PV-Signal-ML: Pharmacovigilance Signal Detection Research Prototype

**Research implementation demonstrating signal detection algorithms, ML-based triage, SISA machine unlearning, and regulatory compliance concepts for educational and portfolio purposes.**

---

## ⚠️ IMPORTANT REGULATORY DISCLAIMER

**This is a research/proof-of-concept project, NOT a validated production system.**

- NOT FDA 21 CFR Part 11 validated (no formal GAMP5 validation performed)
- NOT for use with real patient data (demonstration on aggregated/public data only)
- NOT a substitute for commercial PV systems (SAS, Snowflake, validated PV platforms)
- Educational value: illustrates how enterprise PV systems work algorithmically
- Portfolio project: showcases signal detection, ML, unlearning, and compliance concepts

For real-world use, a separate, validated implementation is required.

---

## 🎯 What Is This?

`pv-signal-ml` is a **research prototype** that demonstrates how enterprise pharmacovigilance (PV) systems can implement:

- Statistical signal detection (PRR, Chi-square) on aggregated safety data
- ML-based triage with XGBoost and SHAP explainability
- RAG-style narrative generation of Signal Assessment Reports (SARs) using a local LLM via Ollama
- GDPR-oriented governance: data lineage, audit logging, and machine unlearning (SISA)

Status: **RESEARCH PROTOTYPE** (non-production, non-validated).

---

## 🏗️ Enterprise Mapping

| Layer | Enterprise Standard | This Prototype | Why This Works |
| --- | --- | --- | --- |
| Data Lake | Snowflake / Databricks | SQLite + CSV + Parquet | Same relational schema and queries; scalability is an infra choice. |
| Statistics | SAS / R | Python (pandas, numpy, scipy) | PRR and Chi-square formulas are identical across tools. |
| ML Engine | SageMaker / Vertex AI | Local XGBoost + MLflow | Algorithm and tracking logic are equivalent; deployment differs. |
| Context (RAG) | GraphRAG / Neo4j | Direct Ollama API | Provides contextual SAR narratives without graph DB for MVP. |
| UI | React/Angular | Streamlit (app_enhanced.py) | UI shell differs; workflow and logic are the same. |
| Compliance | Full GxP stack | Governance docs + lineage + audit logs | Concepts are implemented at prototype level for learning. |

---

## 🛠️ Tech Stack (Actual)

- **Data & Stats:** pandas, numpy, scipy, SQLite, CSV/Parquet
- **ML & Explainability:** XGBoost, scikit-learn, SHAP, MLflow
- **Unlearning:** SISA ensemble (sharded XGBoost models)
- **LLM / RAG:** Ollama API (direct calls, no LangChain, no ChromaDB)
- **UI:** Streamlit single-page app (`app_enhanced.py`)
- **Compliance Tooling:** data lineage JSONs, GDPR deletion registry, audit logging

---

## 🚀 Quick Start

### Prerequisites

- Python 3.9+
- Git
- Ollama installed and a local model (for example `llama3.2:3b`)

### Installation

```bash
git clone <repo-url>
cd PV_Signal_ML

python -m venv venv
# Windows
venv\Scripts\activate
# macOS/Linux
# source venv/bin/activate

pip install -r requirements.txt

# Pull Ollama model (one-time)
ollama pull llama3.2:3b

Run the Streamlit App

streamlit run app_enhanced.py

Then open http://localhost:8501 in your browser.

🧮 Core Workflows
FAERS ingestion & signal computation: faers_build_signals.py

Statistical engine (PRR / Chi-square): stats_engine.py

ML triage + MLflow tracking: pv_signal_ml_pipeline.py

Explainability (SHAP): shap_analysis_simple.py and Tab 3 in app_enhanced.py

SAR generation with Ollama: sar_generator.py and Tab 4 in app_enhanced.py

GDPR deletion registry: gdpr_deletion_registry.py

Audit logging: audit_logging.py

The main user-facing workflow is through app_enhanced.py, which exposes tabs for signal detection, ML validation, explainability, SAR generation, MLflow runs, and SISA unlearning.

📁 Project Structure (Simplified)

PV_Signal_ML/
├── app_enhanced.py              # Main Streamlit app (all tabs, including SISA)
├── api.py                       # FastAPI service (optional)
│
├── faers_build_signals.py       # FAERS ingestion & signal computation
├── stats_engine.py              # PRR / Chi-square calculations
├── prepare_ml_features.py       # Feature engineering for ML
├── pv_signal_ml_pipeline.py     # XGBoost training + MLflow
├── shap_analysis_simple.py      # Standalone SHAP analysis
│
├── sar_generator.py             # Direct Ollama SAR generator
├── rag_signal_evidence.py       # Evidence retrieval (embeddings + PubMed)
├── signal_report_builder.py     # SAR / PSMF formatting
│
├── data_lineage.py              # Data lineage and provenance
├── gdpr_deletion_registry.py    # Right-to-be-forgotten registry
├── audit_logging.py             # Access/event logging
├── governance_dpia.md           # DPIA and governance notes
├── change_control.py            # Predetermined change control plan
│
├── templates/
│   └── signal_report_template.md
│
├── sar_reports/
│   └── reports/                 # Generated SAR markdown/JSON
├── ml_data/                     # ML feature matrices
├── lineage/                     # Lineage JSONs
├── audit_logs/                  # Log files
├── gdpr_registry/               # Deletion registry files
│
├── models/
│   └── sisa/                    # SISA shard models (shard_*.pkl)
│
├── Experimental/                # Experimental UIs / notebooks
└── requirements.txt


🔐 SISA Machine Unlearning (Tab 6)
SISA (Sharded, Isolated, Sliced, Aggregated) is used here to support machine unlearning: efficiently removing the influence of specific cases from an ensemble without full retraining.

Conceptual flow

graph TB
    A[Training data<br/>N cases] --> B[Split into k shards]
    B --> C0[Shard 0 model]
    B --> C1[Shard 1 model]
    B --> C2[Shard 2 model]
    B --> C3[...]
    B --> C9[Shard k-1 model]

    C0 --> D[Ensemble prediction]
    C1 --> D
    C2 --> D
    C3 --> D
    C9 --> D

    E[Unlearn request<br/>case_id] --> F{Find shard<br/>containing case}
    F --> G[Retrain that shard<br/>without case]
    G --> H[Replace shard model<br/>in ensemble]
    H --> D

Implementation in this repo
src/ml/sisa_trainer.py – SISATrainer class with train(...) and unlearn(case_id)

models/sisa/shard_*.pkl – one XGBoost model per shard

app_enhanced.py – Tab 2 ("ML Validation (SISA)") and Tab 6 ("Machine Unlearning") provide UI to train shards and submit unlearning requests

Example training usage (script or notebook):

from src.ml.sisa_trainer import SISATrainer

trainer = SISATrainer(model_dir="models/sisa")
results = trainer.train(signals_df, n_shards=10)
print(results["auc"], results["n_shards"])

Unlearning a specific case:

result = trainer.unlearn(case_id=5432)
print(result)


This retrains only the affected shard and updates the ensemble.

Basic smoke tests via running:

python faers_build_signals.py 2024-01-01 2024-03-31
python pv_signal_ml_pipeline.py
python gdpr_deletion_registry.py
python audit_logging.py


Optional tests (to be expanded):

python -m pytest tests/

Formulas (PRR, Chi-square) are aligned with standard pharmacovigilance guidelines, but the system as a whole is not validated for regulatory use.

📚 References
The implementation is conceptually aligned with:

EMA GVP Module IX (Signal Management)

CIOMS XIV (Practical Aspects of Signal Detection)

Standard PRR / Chi-square disproportionality methods

XGBoost and SHAP documentation for ML and explainability

This repository is intended purely for research and education.
