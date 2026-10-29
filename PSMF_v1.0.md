# MOCK PSMF ANNEX: SIGNAL MANAGEMENT SYSTEM (RESEARCH PROTOTYPE)

> [!WARNING]
> **EDUCATIONAL & RESEARCH USE ONLY — NOT A REGULATORY SUBMISSION**
> This document is an illustrative educational template demonstrating how an automated signal detection pipeline could theoretically be documented within a Pharmacovigilance System Master File (PSMF) annex under EMA GVP Module I / CIOMS XIV guidelines. It does NOT represent a validated commercial or clinical pharmacovigilance system, and has NOT undergone formal FDA 21 CFR Part 11 validation or regulatory inspection.

---

## 1.7.1 Signal Detection & Management Prototype

**System Name**: `pv-signal-ml` (Academic & Open Research Prototype)  
**Status**: Research Prototype (Non-Production / Unvalidated)  
**Target Environment**: Local workstation / Offline research testbed

### Core Components
- **Statistical Disproportionality Engine**: Calculates PRR, ROR, $\chi^2$ (with Yates' correction), IC (BCPNN empirical Bayes), and EBGM/EB05 (MGPS empirical Bayes) on openFDA / public FAERS data files.
- **Machine Learning Layer**: SISA (Sharded, Isolated, Sliced, Aggregated) ensemble architecture for reproducible training and privacy-preserving machine unlearning experimentation (*Bourtoule et al., 2021*).
- **RAG & Synthesis Module**: Local LLM prompt templating using Ollama (Llama 3.2 under Meta Community License) and PubMed E-utilities literature mining. Does not require external proprietary vector databases or cloud APIs.
- **Signal Summary Reports**: Structured markdown generation inspired by CIOMS Working Group XIV recommendations for exploratory signal evaluation.

### Architecture Flow
```
Public FAERS Data → Disproportionality Statistics → SISA ML Classifier → PubMed Miner / Local LLM → Exploratory SAR Reports
```

### Regulatory & Compliance Context
- **EMA GVP Module I / CIOMS XIV**: Prototype designed to reflect the concepts and terminology of international pharmacovigilance standards for research and educational purposes.
- **FDA 21 CFR Part 11 & GxP**: Not formally validated. Production regulatory deployment requires formal audit trails, electronic signature controls, and system validation (IQ/OQ/PQ) conducted by qualified pharmacovigilance and quality assurance professionals.
- **GDPR Art. 17(3)(c) Note**: In statutory pharmacovigilance, statutory public health safety monitoring obligations override standard erasure requests. SISA unlearning is included as an academic machine learning demonstration.
