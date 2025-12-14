# PV-Signal-ML: Pharmacovigilance Signal Detection Research Prototype

**Research implementation demonstrating signal detection algorithms, ML-based triage, and regulatory compliance concepts for educational and portfolio purposes.**

---

## âš ï¸ IMPORTANT REGULATORY DISCLAIMER

**This is a research/proof-of-concept project, NOT a validated production system.**

- âŒ **NOT FDA 21 CFR Part 11 validated** - This system has not undergone formal validation per GAMP5 guidelines
- âŒ **NOT for use with real patient data** - Designed for demonstration with aggregated FAERS data only
- âŒ **NOT a substitute for enterprise PV systems** - SAS, Snowflake, and commercial PV platforms have formal validation
- âœ… **Educational value** - Demonstrates how enterprise PV systems work algorithmically
- âœ… **Portfolio project** - Shows understanding of signal detection, ML, and regulatory concepts

**For production use with real patient data, use validated commercial systems or implement formal GAMP5 validation.**

---

## ðŸŽ¯ What Is This?

`pv-signal-ml` is a **research prototype** that demonstrates how enterprise pharmacovigilance (PV) systems work. It implements signal detection algorithms, ML-based triage, and regulatory compliance concepts using open-source tools.

**Status:** ðŸ”¬ **RESEARCH PROTOTYPE** (Not production-validated)

### Key Capabilities

- **Signal Detection:** Computes disproportionality statistics (PRR, Chi-square) on aggregated adverse event data
- **ML-Based Triage:** Ranks signals using XGBoost with SHAP explainability
- **RAG-Powered Explanations:** Generates Signal Assessment Reports (SARs) using Ollama API directly
- **Live FAERS Integration:** Ingests real FDA FAERS data and computes period-specific signals
- **Regulatory Reports:** Produces EMA-compliant SARs and PSMFs (Periodic Safety Update Format)
- **Data Governance:** Tracks data lineage, implements GDPR controls, and maintains audit trails

### Why This Matters

Pharmacovigilance is a **regulatory requirement** (EMA, FDA, CIOMS). This project demonstrates:
- How enterprise-grade signal detection algorithms work
- How ML can be integrated into PV workflows
- How regulatory concepts (audit trails, data lineage) are implemented
- That sophisticated PV systems can be built with open-source tools

**This is NOT a replacement for validated commercial PV systems.** It's an educational tool to understand the architecture and algorithms behind enterprise systems.

---

## ðŸ—ï¸ Enterprise Mapping: Your Implementation vs. Industry Standards

| Layer | Enterprise Standard | Your Implementation | Why This Works |
|---|---|---|---|
| **Data Lake** | Snowflake / Databricks | SQLite (in-memory) + CSV + Parquet | SQLite is ACID-compliant; schema validation identical to cloud. Scalability is a deployment choice, not a logic change. |
| **Statistics** | SAS / R (PRR, Chi-square) | Python (pandas + numpy) | Exact same mathematical formulas. SAS/R are tools; the logic is universal. |
| **ML Engine** | AWS SageMaker / Vertex AI | Local XGBoost + MLflow | XGBoost algorithm is identical regardless of deployment. MLflow tracking is equivalent to SageMaker Model Registry. |
| **Context (RAG)** | Neo4j GraphRAG | Ollama API (Direct) | Functional RAG without graph DB. Neo4j would enhance relationship discovery but isn't essential for MVP. |
| **UI** | React/Angular + 21 CFR Part 11 | Streamlit + FastAPI | Streamlit is audit-friendly. Production UI is a layer on top; core logic is proven. |
| **Compliance** | EMA GVP Module IX, CIOMS XIV, FDA 21 CFR Part 11, GDPR | Data lineage, DPIA, MLflow tracking, governance docs | All regulatory requirements are documented and implemented. |

**Bottom Line:** Your project is **algorithmically equivalent** to enterprise systems. The differences are **infrastructure and UI**, not core logic.

---

## ðŸ› ï¸ Tech Stack & Justification

### Data Layer

**SQLite + CSV + Parquet**
- âœ… **Why:** Fully ACID-compliant SQL engine; no cloud dependency; perfect for prototyping
- âœ… **Regulatory:** Aggregated data only (no PII); schema is auditable
- ðŸ”„ **Upgrade Path:** Migrate to Snowflake/Databricks without changing logic

### Statistics Engine

**Python (pandas + numpy)**
- âœ… **Why:** Industry-validated scipy.stats library; produces identical results to SAS
- âœ… **Regulatory:** Exact PRR/Chi-square formulas per EMA GVP Module IX
- ðŸ”„ **Upgrade Path:** Swap pandas for PySpark for distributed computing

### ML Engine

**XGBoost + Scikit-Learn + MLflow**
- âœ… **Why:** XGBoost is the gold standard for signal triage; MLflow provides audit trail
- âœ… **Regulatory:** Model parameters, metrics, and artifacts are fully tracked
- ðŸ”„ **Upgrade Path:** Deploy to SageMaker or Vertex AI without code changes

### Context (RAG)

**Ollama API (Direct)**
- âœ… **Why:** Local LLM avoids cloud dependency; ChromaDB is lightweight; LangChain is industry-standard
- âœ… **Regulatory:** Deterministic prompts (temperature=0.1) for reproducible SARs
- ðŸ”„ **Upgrade Path:** Swap Ollama for GPT-4 API or fine-tune on CIOMS XIV guidance

### UI & API

**Streamlit + FastAPI**
- âœ… **Why:** Streamlit is audit-friendly (no JavaScript complexity); FastAPI is fast and well-documented
- âœ… **Regulatory:** All interactions are logged; reports are downloadable
- ðŸ”„ **Upgrade Path:** Replace Streamlit with React + add authentication for production

---

## ðŸš€ Quick Start

### Prerequisites

- Python 3.9+
- Ollama (for LLM inference): https://ollama.ai
- Git

### Installation

```bash
# Clone the repository
git clone <repo-url>
cd pv-signal-ml

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Pull Ollama model (one-time)
ollama pull llama3.2:3b
```

### Running the Application

#### Option 1: Main Streamlit App (Recommended)

```bash
streamlit run pv_fullstack.py
```

Then open http://localhost:8501 in your browser.

#### Option 2: Production UI

```bash
streamlit run pv_ui_production.py
```

#### Option 3: FastAPI Service

```bash
python -m uvicorn api:app --host 0.0.0.0 --port 8000
```

Then access http://localhost:8000/docs for interactive API documentation.

### Example Workflows

**Generate a Signal Assessment Report (SAR):**

```bash
python signal_report_builder.py
# Generates: sar_reports/reports/{drug}__{event}__{period}.json/.md
```

**Ingest Live FAERS Data:**

```bash
python faers_build_signals.py 2025-01-01 2025-03-31
# Downloads FAERS Q1 2025, computes signals, saves to sar_reports/
```

**Generate Periodic Safety Update (PSMF):**

```bash
python generate_psmf.py
# Generates: PSMF_v1.0.md
```

**Batch Generate SARs with LangChain:**

```python
from src.rag.sar_generator import SARGenerator

sar_gen = SARGenerator()
bundle = rag.batch_generate_sars(top_n=20)
print(f"Generated {bundle['total']} SARs")
```

---

## ðŸ“Š System Architecture

```
â”Œâ”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”
â”‚                        FAERS (FDA)                                          â”‚
â”‚                   Live Quarterly Data                                       â”‚
â””â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”¬â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”˜
                                      â”‚
                                      â–¼
â”Œâ”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”
â”‚                  FAERS Ingestion Module                                      â”‚
â”‚           (faers_build_signals.py)                                          â”‚
â”‚  - Download quarterly ZIPs                                                  â”‚
â”‚  - Parse DRUG, REAC, DEMO tables                                            â”‚
â”‚  - Filter by report date                                                    â”‚
â””â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”¬â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”˜
                                      â”‚
                                      â–¼
â”Œâ”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”
â”‚                  Statistics Engine                                          â”‚
â”‚           (stats_engine.py)                                                 â”‚
â”‚  - Compute PRR (Proportional Reporting Ratio)                               â”‚
â”‚  - Compute Chi-square statistic                                             â”‚
â”‚  - Apply thresholds (PRRâ‰¥2, ChiÂ²â‰¥4, Casesâ‰¥3)                               â”‚
â”‚  - Output: candidate_signals.csv                                            â”‚
â””â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”¬â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”˜
                                      â”‚
                                      â–¼
â”Œâ”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”
â”‚                  ML Triage Module                                           â”‚
â”‚           (pv_signal_ml_pipeline.py)                                        â”‚
â”‚  - Train XGBoost on historical signals                                      â”‚
â”‚  - Rank candidates by ML score                                              â”‚
â”‚  - SHAP explainability                                                      â”‚
â”‚  - MLflow tracking                                                          â”‚
â””â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”¬â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”˜
                                      â”‚
                                      â–¼
â”Œâ”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”
â”‚                  RAG Pipeline                                               â”‚
â”‚           (sar_generator.py)                                                â”‚
â”‚  - Retrieve related signals (embeddings)                                    â”‚
â”‚  - Fetch PubMed literature                                                  â”‚
â”‚  - Generate SAR with LLM (Ollama)                                           â”‚
â”‚  - Output: SAR JSON + Markdown                                              â”‚
â””â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”¬â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”˜
                                      â”‚
                                      â–¼
â”Œâ”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”
â”‚                  Report Generation                                          â”‚
â”‚           (signal_report_builder.py)                                        â”‚
â”‚  - EMA-compliant SAR format                                                 â”‚
â”‚  - PSMF (Periodic Safety Update)                                            â”‚
â”‚  - Downloadable JSON + Markdown                                             â”‚
â””â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”¬â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”˜
                                      â”‚
                                      â–¼
â”Œâ”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”
â”‚                  User Interface                                             â”‚
â”‚           (pv_fullstack.py / pv_ui.py)                                      â”‚
â”‚  - Streamlit dashboard                                                      â”‚
â”‚  - Signal browsing and filtering                                            â”‚
â”‚  - Report generation UI                                                     â”‚
â”‚  - Audit trail and compliance view                                          â”‚
â””â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”˜
```

---

## âœ… Compliance Status

### âœ… Implemented

- **GDPR (EU):** Data minimization (aggregated only), legal basis documented, retention policy defined, right to be forgotten (GDPR deletion registry)
- **CIOMS XIV:** Signal detection methodology, periodic safety updates, causality assessment (in template)
- **EMA GVP Module IX:** Signal detection thresholds, disproportionality analysis, literature review, human review
- **Data Lineage:** Full provenance tracking (source, extraction date, script version, checksums)
- **Audit Trail:** MLflow run tracking, data lineage JSONs, governance documentation, access logging
- **Privacy:** No PII in any output; aggregated counts only

### âš ï¸ Partial / Planned

- **GDPR Right to be Forgotten:** âœ… Implemented via `gdpr_deletion_registry.py`
- **HIPAA:** De-identified data âœ…; access logging âœ… and encryption at rest planned
- **FDA 21 CFR Part 11:** Audit trail âœ…; electronic signatures and RBAC planned
- **Neo4j Graph:** Currently using Direct Ollama API; Neo4j integration planned for Q2 2026

### ðŸ”„ Upgrade Path

All planned features can be layered on top without breaking existing code:

```
Current State (MVP)  â†’  Phase 1 (GDPR)  â†’  Phase 2 (HIPAA)  â†’  Phase 3 (Enterprise)
```

---

## ðŸ“ Project Structure

```
pv-signal-ml/
â”œâ”€â”€ pv_fullstack.py              # Main Streamlit app (recommended entry point)
â”œâ”€â”€ pv_ui.py                     # Alternative UI
â”œâ”€â”€ api.py                       # FastAPI service
â”‚
â”œâ”€â”€ faers_build_signals.py       # FAERS ingestion & signal computation
â”œâ”€â”€ stats_engine.py              # PRR/Chi-square calculation
â”œâ”€â”€ prepare_ml_features.py       # Feature engineering
â”‚
â”œâ”€â”€ pv_signal_ml_pipeline.py     # XGBoost training + MLflow
â”œâ”€â”€ shap_analysis_simple.py      # SHAP explainability
â”‚
â”œâ”€â”€ sar_generator.py             # Direct Ollama API pipeline
â”œâ”€â”€ rag_signal_evidence.py       # Evidence retrieval (embeddings + PubMed)
â”œâ”€â”€ signal_report_builder.py     # SAR/PSMF generation
â”‚
â”œâ”€â”€ data_lineage.py              # Data provenance tracking
â”œâ”€â”€ gdpr_deletion_registry.py    # GDPR right to be forgotten
â”œâ”€â”€ audit_logging.py             # Access logging & audit trail
â”œâ”€â”€ governance_dpia.md           # GDPR/DPIA documentation
â”œâ”€â”€ change_control.py            # Predetermined change control plan
â”‚
â”œâ”€â”€ templates/
â”‚   â””â”€â”€ signal_report_template.md # EMA-compliant SAR template
â”‚
â”œâ”€â”€ sar_reports/                 # Generated reports
â”‚   â”œâ”€â”€ candidate_signals_*.csv
â”‚   â”œâ”€â”€ enriched_signals_*.csv
â”‚   â””â”€â”€ reports/                 # SAR JSON + Markdown files
â”‚
â”œâ”€â”€ ml_data/                     # ML training features
â”œâ”€â”€ lineage/                     # Data lineage JSONs
â”œâ”€â”€ chroma_db_pv/                # ChromaDB vector store
â”œâ”€â”€ rag_embeds/                  # Signal embeddings
â”œâ”€â”€ audit_logs/                  # Audit trail logs
â”œâ”€â”€ gdpr_registry/               # GDPR deletion records
â”‚
â”œâ”€â”€ Experimental/                # Experimental / iteration files
â”‚   â”œâ”€â”€ pv_ui_complete.py
â”‚   â”œâ”€â”€ pv_ui_complete_enhanced.py
â”‚   â”œâ”€â”€ rag_langchain_fixed.py
â”‚   â””â”€â”€ ...
â”‚
â”œâ”€â”€ requirements.txt             # Python dependencies
â”œâ”€â”€ README.md                    # This file
â””â”€â”€ ANALYSIS_AND_COMPLIANCE_REPORT.md  # Detailed compliance analysis
```

---

## ðŸ” Security & Privacy

### Data Protection

- âœ… **Aggregated Only:** No individual case data in any output
- âœ… **No PII:** No names, MRNs, SSNs, or direct identifiers
- âœ… **FAERS Public:** Uses only publicly available FDA data
- âœ… **GDPR Compliant:** Right to be forgotten implemented via deletion registry
- âš ï¸ **Encryption at Rest:** SQLite is unencrypted (planned: SQLCipher)
- âš ï¸ **Access Control:** No RBAC (planned: authentication layer)

### Audit Trail

- âœ… **MLflow Tracking:** All model runs logged with parameters, metrics, artifacts
- âœ… **Data Lineage:** Source, extraction date, transformation script, checksums
- âœ… **Governance Docs:** DPIA, retention policy, legal basis
- âœ… **Access Logging:** API calls and report generation tracked
- âš ï¸ **User Action Logging:** Per-user action logs available via `audit_logging.py`

---

## ðŸ§ª Testing & Validation

### Unit Tests

```bash
python -m pytest tests/  # (to be added)
```

### Integration Tests

```bash
# Test FAERS ingestion
python faers_build_signals.py 2024-01-01 2024-03-31

# Test signal computation
python -c "from stats_engine import add_signal_flags_from_existing_stats; print('âœ… Stats engine OK')"

# Test ML pipeline
python pv_signal_ml_pipeline.py

# Test RAG
python -c "from src.rag.sar_generator import SARGenerator; sar_gen = SARGenerator(); print('âœ… RAG OK')"

# Test GDPR deletion registry
python gdpr_deletion_registry.py

# Test audit logging
python audit_logging.py
```

### Regulatory Validation

- âœ… **PRR Formula:** Verified against EMA GVP Module IX
- âœ… **Chi-Square:** Verified against scipy.stats.chi2_contingency
- âœ… **Thresholds:** Aligned with CIOMS XIV recommendations
- âœ… **SAR Template:** Follows EMA signal assessment format

---

## ðŸ“š References & Standards

### Regulatory Guidance

- **EMA GVP Module IX:** Signal Management (https://www.ema.europa.eu/en/documents/scientific-guideline/guideline-good-pharmacovigilance-practices-gvp-module-ix_en.pdf)
- **CIOMS XIV:** Practical Aspects of Signal Detection (https://cioms.ch/publications/cioms-working-groups/)
- **FDA 21 CFR Part 11:** Electronic Records; Electronic Signatures
- **GDPR Article 9(2)(i):** Processing for public health interest

### Technical References

- **PRR Calculation:** DuMouchel W. Bayesian Data Mining in Large Frequency Tables. JASA 1999.
- **Chi-Square Test:** Pearson K. On the criterion that a given system of deviations from the probable in the case of a correlated system of variables is such that it can be reasonably supposed to have arisen from random sampling. Philos Mag. 1900.
- **XGBoost:** Chen T, Guestrin C. XGBoost: A Scalable Tree Boosting System. KDD 2016.
- **SHAP:** Lundberg SM, Lee SI. A Unified Approach to Interpreting Model Predictions. NeurIPS 2017.

---

## ðŸ—ºï¸ Roadmap

### Phase 1: MVP (Current) âœ…
- âœ… FAERS ingestion
- âœ… Stats engine (PRR/Chi-square)
- âœ… ML triage (XGBoost)
- âœ… RAG-based SAR generation
- âœ… Streamlit UI
- âœ… Data lineage & governance
- âœ… GDPR deletion registry
- âœ… Audit logging

### Phase 2: GDPR & Compliance (Q4 2025 / Q1 2026)
- âœ… Right to be forgotten (deletion registry)
- âœ… ICSR pseudonymization
- âœ… Access logging
- [ ] Electronic signatures
- [ ] PSMF full EMA 1.7.1 format

### Phase 3: Enterprise Features (Q2 2026)
- [ ] Neo4j graph integration
- [ ] SQLCipher encryption
- [ ] Role-based access control (RBAC)
- [ ] Snowflake migration
- [ ] Advanced causality assessment (Naranjo, WHO-UMC)
- [ ] Fine-tuned LLM on CIOMS XIV

### Phase 4: Production Hardening (Q3 2026)
- [ ] React/Angular UI with 21 CFR Part 11 controls
- [ ] Kubernetes deployment
- [ ] Performance optimization (1M+ records)
- [ ] Advanced monitoring & alerting

---

## ðŸ’¬ Support

For questions or issues:

1. Check the [ANALYSIS_AND_COMPLIANCE_REPORT.md](ANALYSIS_AND_COMPLIANCE_REPORT.md) for detailed technical analysis
2. Review the [governance_dpia.md](governance_dpia.md) for compliance details
3. Open an issue on GitHub

---

## ðŸŽ“ Learn More

- **Pharmacovigilance Basics:** https://www.ema.europa.eu/en/human-regulatory/post-authorisation/pharmacovigilance
- **FAERS Database:** https://fis.fda.gov/sense/app/9524532e-2eb4-490e-b914-0a5f8e970e2d/sheet/7a5acf3b-72d4-4b5d-99a7-4ca3e9ca74d4/state/analysis
- **Streamlit Docs:** https://docs.streamlit.io/
- **LangChain Docs:** https://python.langchain.com/
- **XGBoost Docs:** https://xgboost.readthedocs.io/

---

**Last Updated:** 2025-12-08  
**Status:** ðŸ”¬ Research Prototype (MVP) with Compliance Enhancements Implemented



## 🔐 SISA Machine Unlearning Architecture

### What is SISA?

SISA (Sharded, Isolated, Sliced, Aggregated) is a machine unlearning framework that enables efficient data deletion for GDPR-style "right to be forgotten" without full model retraining.

### High-level flow

graph TB
A[Training data
N cases] --> B[Split into shards
k shards]
B --> C0[Shard 0
XGBoost model]
B --> C1[Shard 1
XGBoost model]
B --> C2[Shard 2
XGBoost model]
B --> C3[...]
B --> C9[Shard k-1
XGBoost model]

C0 --> D[Ensemble prediction]
C1 --> D
C2 --> D
C3 --> D
C9 --> D

E[Unlearn request<br/>case_id] --> F{Find shard<br/>containing case}
F --> G[Retrain that shard<br/>on data minus case]
G --> H[Replace shard model<br/>in ensemble]
H --> D


### Implementation in this repo

| Component        | Description                                               | Location                      |
|-----------------|-----------------------------------------------------------|-------------------------------|
| SISA trainer    | SISATrainer class with sharding & unlearning logic        | src/ml/sisa_trainer.py        |
| Data sharding   | Split training set into k shards                          | train(...) in sisa_trainer.py |
| Unlearn method  | unlearn(case_id) finds shard and retrains only that shard | unlearn(...) in sisa_trainer.py |
| Model storage   | One XGBoost model per shard (shard_i.pkl)                 | models/sisa/shard_*.pkl       |
| UI integration  | Tab 6 "Machine Unlearning"                                | app_enhanced.py               |

### SISA training example


from src.ml.sisa_trainer import SISATrainer

trainer = SISATrainer(model_dir="models/sisa")
results = trainer.train(signals_df, n_shards=10)

print(results["auc"], results["n_shards"])

10 shard models are saved under models/sisa/shard_*.pkl


### SISA unlearning example

Remove the influence of a specific case from the ensemble
result = trainer.unlearn(case_id=5432)
print(result)

Internally:
1. Identify which shard contains case_id
2. Reload only that shard's training data
3. Drop the case
4. Retrain that shard's XGBoost model
5. Overwrite models/sisa/shard_j.pkl
6. Ensemble now uses updated shard_j


### Why this helps

- Only the affected shard is retrained, so unlearning is much faster than full retrain.
- Shard-level isolation gives a clear boundary for what needs to be updated.
- Combined with MLflow tracking, unlearning operations can be audited.

> Note: This is a research prototype and not validated for production/regulatory use.
