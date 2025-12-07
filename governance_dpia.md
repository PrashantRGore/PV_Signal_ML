# PV Signal Detection DPIA / Governance

**Generated:** 2025-12-06T02:03:04.157697

## Purpose
Pharmacovigilance signal detection digital twin

## Legal Basis
Public health interest (Art. 9(2)(i) GDPR)

## Data Categories
- Aggregated case counts
- Pseudonymised drug-event pairs
## Privacy Measures
- No patient identifiers retained
- Aggregated counts only (A,B,C,D 2x2 tables)
- No free text/narratives processed
- FAERS public domain data
## Retention Policy
FAERS quarterly data + 10 years post-product lifecycle

## Human In Loop
ML triage supports stats engine; humans validate signals

## Audit Trail
MLflow + data lineage JSONs + Git commit history

