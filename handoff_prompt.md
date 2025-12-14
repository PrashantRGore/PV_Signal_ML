# PV Signal ML Project Handoff

## Project Context
I'm continuing development of **PV Signal ML - Pharmacovigilance Signal Detection**, a production-grade pharmacovigilance signal detection system with ML-based causality assessment.

**Last Update:** December 14, 2025 1:51 PM IST
**Working Directory:** `C:\Users\koreo\PV_Signal_ML`
**Environment:** Windows 10/11, Python 3.13, Conda base

---

## Architecture Overview

### Two-Stage Pipeline

**Stage 1: Signal Detection** ✅ Implemented and working
- Method: Statistical (PRR, Chi-square, ROR)
- Input: FAERS adverse event data
- Output: Drug-event signal pairs

**Stage 2: Causality Assessment** ✅ Integrated, needs model training
- Method: ML (SISA-trained model)
- Input: Detected signals + features
- Output: Causality probability scores

**Stage 3: Explainability** ⚠️ Partially implemented
- Method: SHAP values
- Input: ML predictions
- Output: Feature importance

### Key Components
- **data_ingestion**: DataSourceManager (Live FAERS + Demo dataset)
- **signal_detection**: DisproportionalityAnalysis (PRR, Chi², ROR)
- **causality_ml**: CausalityScorer + SISATrainer (GDPR-compliant)
- **drug_filtering**: Portfolio-based targeted surveillance
- **ui**: Streamlit multi-tab dashboard

---

## Current Implementation Status

### ✅ Completed
✅ Live FAERS quarterly data download (2025 Q1)
✅ Statistical signal detection (PRR, Chi², ROR)
✅ Drug portfolio filtering (4 drugs loaded)
✅ SISA model architecture for machine unlearning
✅ CausalityScorer module for ML integration
✅ Multi-tab Streamlit interface
✅ Data source manager (Live/Demo switch)
✅ CSV export functionality

### 🔄 In Progress
🔄 SHAP explainability tab (needs model fix)
🔄 Causality ML integration (code ready, needs testing)
🔄 Model training on live FAERS data

### ⏳ Pending
⏳ RAG integration for evidence generation
⏳ Literature mining (PubMed/Embase)
⏳ GDPR right-to-be-forgotten workflow
⏳ Regulatory report generation (PSMF format)
⏳ MLflow experiment tracking

---

## File Structure

### Core Application Files
- `app_enhanced.py`: Main Streamlit application (enhanced with ML)
- `src/data/data_source_manager.py`: Handles Live FAERS and Demo datasets
- `src/analysis/disproportionality_analysis.py`: Statistical signal detection
- `src/ml/sisa_trainer.py`: SISA sharding for GDPR compliance
- `src/ml/causality_scorer.py`: ML causality assessment (NEW)
- `src/utils/logger.py`: Logging configuration

### Data & Models
- `data/faers_quarterly/`: Downloaded FAERS Q1 2025 data
- `data/processed/`: Processed datasets
- `models/`: Trained ML models (sisa_model.pkl)
- `logs/`: Application and MLflow logs

---

## Technical Stack

**Core:** Python 3.13, Streamlit, Pandas, NumPy
**ML:** scikit-learn, XGBoost, SHAP, (Transformers for future BERT integration)
**Development:** Conda, Git, PowerShell

---

## Data Sources

### Live FAERS Data
- **Period:** 2025Q1 (Jan-Mar 2025)
- **Status:** Downloaded and processed
- **Date Range:** 2025/01/15 to 2025/04/16
- **Records:** ~50M total, filtered by date range

### Drug Portfolio
- **Drugs:** ASPIRIN, IBUPROFEN, ALLOPURINOL, METFORMIN
- **Count:** 4 drugs loaded

---

## ML Models

### SISA Model (Causality Prediction)
- **Path:** `models/sisa_model.pkl`
- **Architecture:** SISA (Sharded, Isolated, Sliced, Aggregated)
- **Performance:** ~0.85 accuracy (reported)
- **Status:** Trained on demo data, needs FAERS retraining
- **Features:** prr, chi2, case_count, ror

**Purpose:** GDPR-compliant machine unlearning for "right to be forgotten" requests

---

## Recent Session Summary

### Changes Made
- Added drug portfolio filtering (Excel upload)
- Integrated SISA machine unlearning architecture
- Created CausalityScorer for Stage 2 ML assessment
- Fixed column display with conditional causality columns
- Added High/Moderate/Low causality metrics
- Prepared SHAP explainability (needs debugging)

### Key Decisions
- Two-stage pipeline: Statistical → ML Causality
- SISA sharding for GDPR right-to-be-forgotten
- Conditional column display for gradual feature rollout
- Portfolio-based filtering for targeted surveillance

---

## Known Issues

### High Priority
- SHAP tab error: 'SISATrainer' object has no attribute 'shard_models'
- Causality columns not showing (needs model training first)

### Medium Priority
- Model training slow on large FAERS dataset
- Memory usage high during signal computation

### Low Priority
- UI styling improvements needed
- Add more comprehensive logging

---

## Next Steps (Priority Order)

### Immediate Tasks
1. Test causality ML integration with live FAERS data
2. Fix SHAP explainability tab (model attribute error)
3. Verify High/Moderate/Low causality metrics display
4. Train SISA model on actual FAERS dataset

### Short-Term Goals
5. Implement RAG for evidence generation
6. Add literature mining (PubMed API)
7. Create GDPR data deletion workflow
8. Generate regulatory report templates

### Long-Term Vision
9. Integrate Drug-Causality-BERT v2
10. Build automated monitoring dashboard
11. Deploy to production environment
12. Add user authentication and audit trails

---

## Where I Need Help

**Current Focus:** Testing the complete two-stage pipeline (Statistical Signal Detection → ML Causality Assessment) with live FAERS data.

**Specific Questions:**
1. How to verify causality scoring is working correctly?
2. Best approach to train SISA model on large FAERS dataset efficiently?
3. Fixing SHAP explainability tab attribute error

**Context:** The app is running locally, causality integration code is in place, but needs testing after restart.

---

## Quick Start Commands

Project directory
cd C:\Users\koreo\PV_Signal_ML

Activate environment
conda activate base

Run application
streamlit run app_enhanced.py

Key files to check
- app_enhanced.py (line 82-105: causality integration)
- src/ml/causality_scorer.py (new module)
- models/sisa_model.pkl (trained model)


---

## Additional Context

- **Regulatory Compliance:** FDA/EMA pharmacovigilance standards, GDPR for data privacy
- **Portfolio Drugs:** Targeted surveillance for company-specific drugs
- **Signal Detection:** Statistical methods (PRR ≥ 2, Chi² ≥ 4) are regulatory standard
- **ML Enhancement:** Causality assessment prioritizes signals for manual review
- **Machine Unlearning:** SISA sharding allows deletion of specific cases without full retraining

Ready to continue from here! Please help with testing the causality integration.
