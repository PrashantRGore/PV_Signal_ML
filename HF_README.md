---
license: cc0-1.0
task_categories:
- tabular-classification
tags:
- pharmacovigilance
- drug-safety
- synthetic-data
- FAERS
- medical
pretty_name: Synthetic FAERS-Style Pharmacovigilance Dataset (1M Records)
size_categories:
- 1M<n<10M
---

# Synthetic FAERS-Style Pharmacovigilance Dataset (1M Records)

## Dataset Description

This is a **fully synthetic** dataset designed for pharmacovigilance signal detection and ML model training. Contains 1,000,000 synthetic Individual Case Safety Reports (ICSRs).

### Features

- **100% Synthetic**: No real patient data - GDPR/HIPAA compliant
- **1M Records**: Large-scale dataset for ML training
- **35 Features**: Causality, temporal, lab values, demographics
- **700+ Drug-Event Signals**: Realistic signal patterns

### Files

- `synthetic_faers_1m_pv_signal_ml.csv` - PV Signal ML compatible (recommended)
- `synthetic_faers_1m_v2.csv` - Original column names

### Quick Start

from datasets import load_dataset

Load the dataset
dataset = load_dataset("koreo/synthetic-faers-1m-v2",
data_files="synthetic_faers_1m_pv_signal_ml.csv",
split="train")

Convert to pandas
df = dataset.to_pandas()
print(f"Loaded {len(df):,} records")


### Schema

**Drug Information:**
- drug_name, dose, route, therapy_duration

**Event Information (MedDRA):**
- llt, pt, hlt, hlgt, soc

**Causality Features:**
- dechallenge, rechallenge, causality, action_drug

**Temporal:**
- ttp_days (time to onset), event_duration

**Demographics:**
- age, sex, weight, country

**Lab Values:**
- alt, ast, bilirubin, creatinine, bun

**Other:**
- concomitant_drugs, medical_history, serious, outcome

### Privacy

✅ GDPR Compliant  
✅ HIPAA Compliant  
✅ 100% Synthetic  
✅ No PII/PHI

### License

CC0-1.0 (Public Domain)
