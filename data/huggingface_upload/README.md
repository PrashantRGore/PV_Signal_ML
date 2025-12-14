# PV-Signal-ML Dataset (1M Synthetic ICSR Records)

## 🎯 Description

This dataset contains **1 million synthetic Individual Case Safety Reports (ICSRs)** for pharmacovigilance signal detection.

**⚠️ NO REAL PATIENT DATA** - Fully GDPR/HIPAA compliant synthetic data.

## 📊 Dataset Statistics

- **Total Records**: 1,000,000
- **Drug-Event Pairs**: 1,000,000
- **Unique Drugs**: 30
- **Unique Events**: 15
- **Medication Errors**: 141,875

## 📦 Files

| File | Records | Description |
|------|---------|-------------|
| drug.parquet | 1,000,000 | Drug exposures (product name, role) |
| 
eac.parquet | 1,000,000 | Adverse events (MedDRA PT terms) |
| demo.parquet | 1,000,000 | Patient demographics (age, sex) |
| outc.parquet | 1,000,000 | Case outcomes |

## 🔑 Key Features

✅ **1M+ drug-event pairs** for signal detection  
✅ **141,875 medication error records**  
✅ **FAERS-compatible schema**  
✅ **MedDRA-coded adverse events**  
✅ **Production-quality data structure**

## 📋 Schema

### drug.parquet
- case_id: Unique case identifier
- drug_name: Product/drug name
- 
ole_cod: Drug role (PS=Primary Suspect)

### reac.parquet
- case_id: Unique case identifier
- pt: MedDRA Preferred Term (adverse event)

### demo.parquet
- case_id: Unique case identifier
- ge: Patient age
- sex: Patient gender

### outc.parquet
- case_id: Unique case identifier
- outcome: Case outcome

## 🚫 Legal Disclaimer

**This dataset is 100% synthetic and contains ZERO real patient data.**

✅ **GDPR Compliant**: No personal identifiable information  
✅ **HIPAA Compliant**: No protected health information  
✅ **MIT License**: Free for research and educational use  

⚠️ **NOT FOR PRODUCTION**: For training, research, and portfolio demonstration only.

## 🎓 Use Cases

- Machine learning model training for pharmacovigilance
- Disproportionality analysis (PRR, ROR, Chi-square)
- SISA (machine unlearning) demonstrations
- SHAP explainability for regulatory compliance
- Benchmarking signal detection algorithms

## 🔗 Related Project

📦 **GitHub**: [PrashantRGore/PV_Signal_ML](https://github.com/PrashantRGore/PV_Signal_ML)

## 📄 Citation

@misc{pv-signal-ml-2025,
author = {Prashant R. Gore},
title = {PV-Signal-ML Dataset: 1M Synthetic ICSR Records},
year = {2025},
publisher = {HuggingFace},
url = {https://huggingface.co/datasets/PrashantRGore/pv-signal-ml-data}
}


---

**Generated**: December 2025  
**License**: MIT  
**Contact**: GitHub Issues
