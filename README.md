# PV-Signal-ML: Open-Source Pharmacovigilance Signal Detection & Machine Learning Platform

**Zero-Cost, Legally Compliant, Local-First Pharmacovigilance Machine Learning System for Signal Detection, Statistical Disproportionality, Causal Inference, and Machine Unlearning (SISA).**

---

## ⚖️ Zero-Cost & Open-Access Legal Architecture

This project is strictly engineered to operate **free of commercial licensing fees and cloud API charges**:

1. **No MedDRA License Required**:
   Adverse event and indication terms are standardized using **openFDA Preferred Terms (PT)**, **NIH NLM RxNorm**, and **Medical Subject Headings (MeSH)** / ICD-10 cross-walks—all 100% public domain and freely accessible without MedDRA MSSO subscription fees.
2. **Zero Cloud LLM Dependencies**:
   Safety Assessment Report (SAR) generation runs entirely on local open-weight models via **Ollama** (`llama3.2`, `mistral`), with an automatic offline fallback to a **Deterministic Clinical Regulatory Template Engine** compliant with EMA GVP Module IX.
3. **No Target Leakage (OMOP Ground Truth Benchmark)**:
   Model evaluation and signal validation decouple ground-truth labels from arbitrary disproportionality thresholds ($PRR \ge 2 \land \chi^2 \ge 4$), utilizing the open-access **OMOP / OHDSI Reference Set** (*Ryan et al., 2013*) comprising curated positive and negative control pairs for AMI, UGIB, AKI, and DILI.

---

## 🎯 Key Capabilities

- **Vectorized Disproportionality Engine**:
  Computes Proportional Reporting Ratio (PRR), Reporting Odds Ratio (ROR) with 95% Confidence Intervals, and Yates' continuity-corrected Chi-Square ($\chi^2$) in milliseconds using NumPy matrix operations.
- **Multivariate Tabular Imputation (MICE & PyTorch DAE)**:
  Mathematically sound missing data imputation using **Multivariate Imputation by Chained Equations (MICE)** via `IterativeImputer` and an authentic **PyTorch Tabular Denoising Autoencoder** with corruption noise and reconstruction loss.
- **Full FAERS Adverse Event Extraction**:
  Preserves all reported adverse reactions per ICSR case (eliminating `.first()` truncation) and extracts concomitant medications directly from non-primary suspect drug records.
- **Algorithmic Regulatory Causality Assessment**:
  Evidence-based implementations of the **Naranjo Adverse Drug Reaction Probability Scale** (10 clinical items evaluated from case fields) and the **WHO-UMC Standardised Causality System** (temporal sequence, dechallenge, rechallenge, confounding risk).
- **SISA Machine Unlearning (GDPR Right to Be Forgotten)**:
  Sharded, Isolated, Sliced, Aggregated (SISA) ensemble training with XGBoost. When an unlearning request (`unlearn(case_id)`) arrives, only the single affected shard is retrained, avoiding expensive full-model re-training.
- **Explainability (SHAP TreeExplainer)**:
  Model-agnostic and tree-based feature attribution revealing individual risk drivers and global feature importance.

---

## 🏗️ Architecture Mapping

| Layer | Enterprise Equivalent | This Platform | Zero-Cost Implementation |
| :--- | :--- | :--- | :--- |
| **Data Ingestion** | SAS / Snowflake | FAERS Enriched Extractor + DuckDB / SQLite | Public FDA ASCII Quarterly Zips & Parquet storage |
| **Medical Terminology** | MedDRA MSSO | openFDA + NIH RxNorm + MeSH | Free US National Library of Medicine (NLM) REST APIs |
| **Disproportionality** | SAS Macros / PhV Tools | Vectorized NumPy Engine | Evans PRR, van Puijenbroek ROR, Yates' $\chi^2$ |
| **Ground Truth** | Proprietary Drug Safety DB | OMOP / OHDSI Reference Benchmark | 398 curated positive & negative controls (Ryan et al. 2013) |
| **Missing Data Imputer** | Commercial Statistical Tool | MICE (`IterativeImputer`) + PyTorch DAE | van Buuren MICE algorithm + PyTorch neural autoencoder |
| **ML Engine** | AWS SageMaker / Vertex AI | XGBoost + MLflow Tracking | Local multicore CPU/GPU execution |
| **Narrative RAG** | Azure OpenAI / Bedrock | Local Ollama + Deterministic Fallback | `llama3.2` / `mistral` local execution + GVP Module IX engine |
| **Compliance & Unlearning** | Manual Database Purge | SISA Trainer Ensemble | Retrains only 1 shard per GDPR unlearning request |

---

## 🛠️ Tech Stack & Requirements

- **Python**: 3.9+ (tested on Python 3.13)
- **Core ML & Numerics**: `numpy`, `pandas`, `scipy`, `scikit-learn`, `xgboost`, `torch`
- **Explainability & Tracking**: `shap`, `mlflow`
- **Pharmacovigilance Standards**: `requests` (NIH RxNav API), `biopython` (NCBI Entrez)
- **Interactive Interfaces**: `streamlit` (UI), `networkx` (Knowledge Graph)
- **Testing**: `pytest`, `pytest-cov`

---

## 🚀 Quick Start

### 1. Clone & Environment Setup

```bash
git clone https://github.com/PrashantRGore/PV_Signal_ML.git
cd PV_Signal_ML

# Create virtual environment
python -m venv venv

# Windows
venv\Scripts\activate

# macOS / Linux
source venv/bin/activate

# Install dependencies
pip install -r requirements.txt
```

### 2. Run Automated Verification Tests

Verify the statistical engines, OMOP benchmark, imputer, causality assessor, and SISA unlearning:

```bash
python -m pytest
```
*Expected: 15 passed unit tests.*

### 3. Launch Applications

**Primary Production UI (10 Tabs):**
```bash
streamlit run app_production_v3.py
```

**Alternative Enhanced UI (6 Tabs):**
```bash
streamlit run app_enhanced.py
```

Open `http://localhost:8501` in your browser.

*(Optional)* To enable local LLM narrative generation, launch Ollama:
```bash
ollama serve
ollama pull llama3.2
```
*Note: If Ollama is offline or not installed, `SARGenerator` automatically falls back to the built-in deterministic regulatory template engine without crashing.*

---

## 📁 Repository Structure

```text
PV_Signal_ML/
├── app_production_v3.py        # Primary Streamlit application (10 comprehensive tabs)
├── app_enhanced.py             # Alternative Streamlit UI (6-tab focused workflow)
├── pv_fullstack.py             # Fullstack entrypoint script
├── config.py                   # Central configuration & directory paths
├── pytest.ini                  # Pytest test discovery configuration
├── requirements.txt            # Python dependencies
├── README.md                   # System documentation
│
├── src/
│   ├── causal/
│   │   ├── causal_features.py       # Causal plausibility feature engineer (deterministic SHA256)
│   │   ├── knowledge_graph.py       # NetworkX biomedical knowledge graph
│   │   └── kg_data/                 # Drug mechanisms, DDIs, and organ tropisms
│   │
│   ├── data/
│   │   ├── canonical_schema.py      # 36-column canonical schema validator
│   │   ├── icsr_to_signals.py       # Vectorized ICSR-to-signal converter with OMOP labeling
│   │   ├── faers_downloader.py      # Resilient atomic FAERS downloader with zip validation
│   │   ├── hf_loader.py             # HuggingFace & local cache dataset loader
│   │   ├── enrichment/
│   │   │   └── faers_enriched_extractor.py  # FAERS ASCII parser preserving all reactions
│   │   └── standards/
│   │       ├── omop_benchmark.py    # OMOP/OHDSI ground truth reference controls
│   │       └── rxnorm_client.py     # NIH RxNav RxNorm API client with SQLite cache
│   │
│   ├── explainability/
│   │   ├── shap_analysis.py         # SHAP TreeExplainer with numeric-only filtering
│   │   └── shap_analyzer_v2.py      # Extended explainability analyzer
│   │
│   ├── features/
│   │   └── build_features_2020.py   # Monthly rolling window aggregation
│   │
│   ├── literature/
│   │   └── pubmed_miner.py          # PubMed E-utilities search with rate limiting
│   │
│   ├── ml/
│   │   ├── sisa_trainer.py          # Functional SISA unlearning ensemble trainer
│   │   ├── causality_scorer.py      # Gradient boosted signal prioritization
│   │   ├── feature_engineering/
│   │   │   └── feature_pipeline.py  # Height-aware BMI & clean polypharmacy parsing
│   │   └── nan_handler/
│   │       └── denoising_autoencoder.py  # MICE (IterativeImputer) & PyTorch Tabular DAE
│   │
│   ├── mlflow_tracking/
│   │   └── tracker.py               # MLflow experiment tracking integration
│   │
│   ├── rag/
│   │   └── sar_generator.py         # Local Ollama + Deterministic GVP Module IX engine
│   │
│   ├── regulatory/
│   │   └── causality_assessor.py    # Algorithmic WHO-UMC & Naranjo scales
│   │
│   ├── stats_engine/
│   │   └── disproportionality.py    # Vectorized Evans PRR, ROR 95% CI & Yates' Chi²
│   │
│   └── utils/
│       ├── pipeline_state.py        # Safe session state & DataFrame truthiness handler
│       ├── logger.py                # Structured system logging
│       └── service_health.py        # Local Ollama & MLflow service probes
│
├── tests/                           # Complete automated unit test suite (15 tests)
│   ├── test_disproportionality.py   # Analytical verification of PRR, ROR & Chi²
│   ├── test_omop_benchmark.py       # Ground-truth positive/negative control verification
│   ├── test_nan_handler.py          # MICE & PyTorch DAE imputation verification
│   ├── test_causality_assessor.py   # WHO-UMC & Naranjo algorithmic scoring tests
│   ├── test_feature_pipeline.py     # BMI calculation & reaction explosion verification
│   ├── test_sisa_unlearning.py      # SISA shard retraining & unlearning tests
│   ├── test_pipeline_state.py       # Pipeline state DataFrame truthiness tests
│   └── test_sar_generator.py        # Offline SAR template generation tests
│
└── data/                            # Local datasets and cache directories
```

---

## 📚 Legal Citations & Scientific References

1. **OMOP Ground Truth Reference Standard**:
   Ryan PB, Schuemie MJ, Gruber S, et al. *Empirical performance evaluation of methods for risk identification in healthcare data: results from the Experiments of the Observational Medical Outcomes Partnership.* Statistics in Medicine 2013; 32(30): 5601–5613. DOI: [10.1002/sim.5869](https://doi.org/10.1002/sim.5869). *(Open-access CC-BY 4.0)*.
2. **Proportional Reporting Ratio (PRR)**:
   Evans SJ, Waller PC, Davis S. *Use of proportional reporting ratios (PRRs) for signal detection from spontaneous adverse drug reaction reports.* Pharmacoepidemiol Drug Saf 2001; 10(6): 483–486. DOI: [10.1002/pds.677](https://doi.org/10.1002/pds.677).
3. **Reporting Odds Ratio (ROR)**:
   van Puijenbroek EP, Bate A, Leufkens HG, et al. *A comparison of methods used for disproportionality analysis.* Pharmacoepidemiol Drug Saf 2002; 11(1): 3–10. DOI: [10.1002/pds.668](https://doi.org/10.1002/pds.668).
4. **Naranjo ADR Probability Scale**:
   Naranjo CA, Busto U, Sellers EM, et al. *A method for estimating the probability of adverse drug reactions.* Clin Pharmacol Ther 1981; 30(2): 239–245. DOI: [10.1038/clpt.1981.154](https://doi.org/10.1038/clpt.1981.154).
5. **WHO-UMC Causality Assessment**:
   The use of the WHO-UMC system for standardised case causality assessment. Uppsala Monitoring Centre. Available at: [who-umc.org](https://who-umc.org).
6. **MICE Imputation Algorithm**:
   van Buuren S, Groothuis-Oudshoorn K. *mice: Multivariate Imputation by Chained Equations in R.* Journal of Statistical Software 2011; 45(3): 1–67. DOI: [10.18637/jss.v045.i03](https://doi.org/10.18637/jss.v045.i03).
7. **SISA Machine Unlearning**:
   Bourtoule L, Chandrasekaran V, Choquette-Choo CA, et al. *Machine unlearning via SISA training.* IEEE Symposium on Security and Privacy (SP) 2021; 141–159. DOI: [10.1109/SP40001.2021.00019](https://doi.org/10.1109/SP40001.2021.00019).
8. **Good Pharmacovigilance Practices (GVP)**:
   European Medicines Agency (EMA). *Guideline on good pharmacovigilance practices (GVP) Module IX – Signal management.* EMA/827661/2011 Rev 1.
9. **US FDA Adverse Event Reporting System (FAERS) & RxNorm**:
   United States Food and Drug Administration (FDA) FAERS Public Database and National Library of Medicine (NLM) RxNav API (Public Domain, US Government).

---

## ⚠️ Regulatory Disclaimer

*This project is a research and educational platform demonstrating algorithmic concepts in pharmacovigilance, statistical disproportionality, causal assessment, and machine unlearning. It is not an FDA 21 CFR Part 11 validated production environment. For regulatory submissions to health authorities (FDA, EMA, PMDA, MHRA), formal GAMP-5 computer system validation (CSV) and audit qualification are required.*
