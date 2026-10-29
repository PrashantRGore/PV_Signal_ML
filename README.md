# PV-Signal-ML: Open-Source Pharmacovigilance Signal Detection & Machine Learning Platform

**Zero-Cost, Local-First Pharmacovigilance Machine Learning Research Prototype for Signal Detection, Statistical Disproportionality, Causal Assessment, and Machine Unlearning (SISA).**

---

## ⚖️ Open-Access Legal Architecture & Regulatory Clarifications

This platform is strictly designed to operate **free of commercial software license fees and paid cloud API charges**, while upholding transparency regarding intellectual property and regulatory standards:

1. **MedDRA® Terminology & FAERS Public Data Distinctions**:
   - The U.S. Food and Drug Administration (FDA) publicly distributes quarterly Adverse Event Reporting System (FAERS) ASCII files and openFDA APIs for academic, research, and educational purposes. These extracts contain standardized reaction strings corresponding to MedDRA Preferred Terms (`PT`).
   - However, **MedDRA®** (Medical Dictionary for Regulatory Activities) is proprietary intellectual property owned by the International Council for Harmonisation (ICH) and administered by the Maintenance and Support Services Organization (MSSO). In formal commercial pharmacovigilance, clinical safety operations, or regulatory submissions, organizations must maintain an official MedDRA subscription license from MSSO.
   - This research prototype operates directly on freely distributed FDA FAERS public research extracts, **NIH NLM RxNorm**, and open biomedical crosswalks without packaging or distributing proprietary MedDRA dictionary files.

2. **Local LLMs & Meta Llama Community License**:
   - Safety Assessment Report (SAR) generation utilizes local open-weight large language models via **Ollama** (`llama3.2`, `mistral`).
   - Llama 3.2 is governed by the **Meta Llama Community License Agreement** (free for research, commercial, and educational applications subject to Meta's Acceptable Use Policy).
   - If Ollama is unavailable or offline, the system automatically falls back to an offline **Deterministic Clinical Regulatory Template Engine** aligned with EMA GVP Module IX principles.

3. **GDPR Article 17 Statutory Exception & SISA Unlearning**:
   - Under the European Union General Data Protection Regulation (GDPR) **Article 17(3)(b) & (c)** and **Article 9(2)(i)**, statutory pharmacovigilance safety monitoring obligations legally supersede patient "right to erasure" (right to be forgotten) requests due to overriding public health interests.
   - The **SISA** (Sharded, Isolated, Sliced, Aggregated) machine unlearning implementation (*Bourtoule et al., 2021*) is provided as an academic machine learning research testbed exploring algorithmic data revocation in ML models, rather than a clinical compliance mechanism for statutory adverse event records.

4. **No Target Leakage (OMOP Reference Standard)**:
   - Ground truth training labels are strictly decoupled from disproportionality metrics ($PRR \ge 2 \land \chi^2 \ge 4$). Screening flags (`is_evans_signal`) remain separate from the target `label`, preventing circular target leakage.
   - Reference evaluation leverages the open-access **OMOP / OHDSI Reference Set** (*Ryan et al., 2012*), featuring **81 curated reference control pairs** (40 positive and 41 negative controls across Acute Myocardial Infarction [AMI], Upper Gastrointestinal Bleeding [UGIB], Acute Kidney Injury [AKI], and Drug-Induced Liver Injury [DILI]).

---

## 🎯 Key Capabilities

- **Vectorized Disproportionality Engine**:
  Computes Proportional Reporting Ratio (PRR), Reporting Odds Ratio (ROR) with 95% Confidence Intervals, and Yates' continuity-corrected Chi-Square ($\chi^2$) in milliseconds using NumPy matrix operations.
- **Multivariate Tabular Imputation (MICE & PyTorch DAE)**:
  Mathematically sound missing data imputation using **Multivariate Imputation by Chained Equations (MICE)** via `IterativeImputer` and a **PyTorch Tabular Denoising Autoencoder** with corruption noise and reconstruction loss.
- **Full FAERS Adverse Event Extraction**:
  Preserves all reported adverse reactions per ICSR case (eliminating truncation) and extracts concomitant medications directly from non-primary suspect drug records.
- **Algorithmic Regulatory Causality Assessment**:
  Evidence-based implementations of the **Naranjo Adverse Drug Reaction Probability Scale** (10 clinical items evaluated from case fields) and the **WHO-UMC Standardised Causality System** (temporal sequence, dechallenge, rechallenge, confounding risk).
- **SISA Machine Unlearning (Research Prototype)**:
  Sharded, Isolated, Sliced, Aggregated (SISA) ensemble training with XGBoost. When an unlearning request (`unlearn(case_id)`) arrives, only the single affected shard is retrained, avoiding expensive full-model retraining.
- **Explainability (SHAP TreeExplainer)**:
  Model-agnostic and tree-based feature attribution revealing individual risk drivers and global feature importance.

---

## 🏗️ Architecture Mapping

| Layer | Enterprise Equivalent | This Platform | Zero-Cost Implementation |
| :--- | :--- | :--- | :--- |
| **Data Ingestion** | SAS / Snowflake | FAERS Enriched Extractor + DuckDB / SQLite | Public FDA ASCII Quarterly Zips & Parquet storage |
| **Medical Terminology** | MedDRA MSSO | openFDA + NIH RxNorm + MeSH | Free US National Library of Medicine (NLM) REST APIs |
| **Disproportionality** | SAS Macros / PhV Tools | Vectorized NumPy Engine | Evans PRR, van Puijenbroek ROR, Yates' $\chi^2$ |
| **Ground Truth** | Proprietary Drug Safety DB | OMOP / OHDSI Reference Benchmark | 81 curated positive & negative controls (Ryan et al. 2012) |
| **Missing Data Imputer** | Commercial Statistical Tool | MICE (`IterativeImputer`) + PyTorch DAE | van Buuren MICE algorithm + PyTorch neural autoencoder |
| **ML Engine** | AWS SageMaker / Vertex AI | XGBoost + MLflow Tracking | Local multicore CPU/GPU execution (pinned `mlflow<3`) |
| **Narrative RAG** | Azure OpenAI / Bedrock | Local Ollama + Deterministic Fallback | `llama3.2` / `mistral` local execution + GVP Module IX engine |
| **Unlearning Exploration** | Manual Database Purge | SISA Trainer Ensemble | Retrains only 1 shard per algorithmic unlearning request |

---

## 🛠️ Tech Stack & Requirements

- **Python**: 3.9+ (tested on Python 3.13)
- **Core ML & Numerics**: `numpy`, `pandas`, `scipy`, `scikit-learn`, `xgboost`, `torch`
- **Explainability & Tracking**: `shap`, `mlflow>=2.8.0,<3.0.0`
- **Pharmacovigilance Standards & Mining**: `requests` (NIH RxNav API), `biopython>=1.83` (NCBI Entrez E-utilities)
- **Interactive Interface**: `streamlit` (UI), `networkx` (Biomedical Knowledge Graph)
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
*Expected: 15 passed unit tests (or 14 passed, 1 skipped if PyTorch is not installed).*

### 3. Launch the Application

```bash
streamlit run app.py
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
├── app.py                      # Canonical Streamlit application
├── config.py                   # Central configuration & directory paths
├── pytest.ini                  # Pytest test discovery configuration
├── requirements.txt            # Python dependencies (pinned mlflow, biopython)
├── README.md                   # System documentation & legal notices
├── VALIDATION_STATUS.md        # System validation matrix & prototype status
├── PSMF_v1.0.md                # Educational PSMF Annex template
├── ARCHITECTURE_v3.md          # Technical architecture documentation
├── DEPLOYMENT_GUIDE.md         # Deployment & verification guide
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
│   │       ├── omop_benchmark.py    # OMOP/OHDSI ground truth reference controls (81 pairs)
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
│   │   ├── sisa_trainer.py          # SISA unlearning ensemble trainer (XGBoost)
│   │   ├── causality_scorer.py      # Prioritized causality scoring
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
   Ryan PB, Madigan D, Stang PE, Overhage JM, Racoosin JA, Hartzema AG. *Empirical assessment of methods for risk identification in healthcare data: results from the experiments of the Observational Medical Outcomes Partnership.* Statistics in Medicine 2012; 31(30): 4401–4415. DOI: [10.1002/sim.5620](https://doi.org/10.1002/sim.5620). *(The OMOP benchmark reference dataset comprising 81 positive and negative control pairs is openly released by the OHDSI / OMOP research community).*
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
10. **Meta Llama Community License**:
    Meta Platforms, Inc. *Llama 3.2 Community License Agreement.* Available at: [llama.meta.com/license](https://llama.meta.com/license).
11. **MedDRA® Trademark Notice**:
    MedDRA® is a registered trademark of the International Federation of Pharmaceutical Manufacturers and Associations (IFPMA) on behalf of the International Council for Harmonisation (ICH).

---

## ⚠️ Regulatory Disclaimer

*This project is an open-source research and educational prototype demonstrating algorithmic concepts in pharmacovigilance, statistical disproportionality, causal assessment, and machine unlearning. It is not an FDA 21 CFR Part 11 or GAMP-5 validated production environment. For regulatory submissions to health authorities (FDA, EMA, PMDA, MHRA), formal computer system validation (CSV), installation qualification (IQ), operational qualification (OQ), and performance qualification (PQ) are required.*
