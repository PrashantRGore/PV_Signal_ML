# Synthetic FAERS-Style Pharmacovigilance Dataset (1M Records)

## Dataset Description

This is a **fully synthetic** dataset designed to mimic the structure and complexity of FDA Adverse Event Reporting System (FAERS) Individual Case Safety Reports (ICSRs). It contains 1,000,000 synthetic adverse event reports for training and testing pharmacovigilance signal detection and machine learning models.

### Key Features

- **100% Synthetic**: No real patient data. Fully compliant with GDPR, HIPAA, and EMA privacy regulations
- **71 Unique Event Terms**: Mapped across MedDRA-like hierarchy (LLT → PT → HLT → HLGT → SOC)
- **25 Fictional Drugs**: Company suspect drugs with realistic pharmacovigilance parameters
- **~700 Drug-Event Signals**: Complex signal patterns for ML model training
- **Rich Feature Set**: 36 columns including causality assessment, dechallenge/rechallenge, concomitant medications, medical history, lab values, and temporal relationships

### Privacy Compliance

✅ **GDPR Compliant**: No personal identifiable information (PII)  
✅ **HIPAA Compliant**: No protected health information (PHI)  
✅ **K-Anonymity**: Minimum group size >= 5 for quasi-identifiers  
✅ **Differential Privacy**: Statistical noise added to sensitive attributes  
✅ **EMA GVP Compliant**: All required masking applied per Module VI Addendum II

**What is NOT included (per privacy regulations):**
- No patient names, initials, or medical record numbers
- No reporter names, addresses, phone numbers, or emails
- No hospital names or specific facilities
- No exact geographic locations (city/zip code) - country level only
- No dates of birth - only age ranges
- No FAERS case numbers or real regulatory identifiers

### Dataset Schema

| Column | Description | Data Type | Missing % |
|--------|-------------|-----------|-----------|
| `case_id` | Unique synthetic case identifier (SHA256 hash) | String | 0% |
| `receive_date` | Synthetic report receipt date | Date | 0% |
| `country` | Country of occurrence (ISO 3-letter code) | String | 0% |
| `age` | Patient age in years (with DP noise) | Integer | 0% |
| `age_group` | Regulatory age category | String | 0% |
| `sex` | Patient sex | String | 0% |
| `weight_kg` | Patient weight in kg | Float | ~35% |
| `suspect_drug` | Fictional company suspect drug name | String | 0% |
| `indication` | Drug indication | String | 0% |
| `route` | Route of administration | String | 0% |
| `dose` | Dose amount | Integer | 0% |
| `dose_unit` | Dose unit (mg, mcg, etc.) | String | 0% |
| `dose_frequency` | Dosing frequency | String | ~20% |
| `treatment_duration_days` | Duration of treatment | Integer | ~30% |
| `event_llt` | Adverse event (Lowest Level Term) | String | 0% |
| `event_pt` | Adverse event (Preferred Term) | String | 0% |
| `event_hlt` | Adverse event (High Level Term) | String | 0% |
| `event_hlgt` | Adverse event (High Level Group Term) | String | 0% |
| `event_soc` | Adverse event (System Organ Class) | String | 0% |
| `time_to_onset_days` | Time from drug start to event onset | Integer | 0% |
| `event_duration_days` | Duration of adverse event | Integer | ~40% |
| `seriousness` | Seriousness criteria (ICH E2B) | String | 0% |
| `outcome` | Event outcome | String | 0% |
| `action_taken` | Action taken with suspect drug | String | 0% |
| `dechallenge` | Dechallenge result | String | 0% |
| `rechallenge` | Rechallenge result | String | 0% |
| `concomitant_medications` | List of concomitant medications | String | 0% |
| `medical_history` | Relevant medical history | String | 0% |
| `reporter_type` | Reporter qualification | String | 0% |
| `report_type` | Type of report source | String | 0% |
| `causality_assessment` | WHO-UMC style causality category | String | 0% |
| `alt_u_l` | ALT lab value (U/L) | Float | ~25% |
| `ast_u_l` | AST lab value (U/L) | Float | ~25% |
| `bilirubin_mg_dl` | Total bilirubin (mg/dL) | Float | ~25% |
| `creatinine_mg_dl` | Serum creatinine (mg/dL) | Float | ~25% |
| `bun_mg_dl` | Blood urea nitrogen (mg/dL) | Float | ~25% |

### Use Cases

- **ML Model Training**: Train drug safety signal detection models (PRR, BCPNN, MGPS, BERT-based)
- **Algorithm Development**: Test new pharmacovigilance algorithms without privacy concerns
- **Education**: Learn ICSR structure and regulatory requirements
- **Benchmarking**: Compare model performance on standardized dataset
- **RAG Systems**: Literature mining and knowledge graph integration

### Known Limitations

- **Not Real FAERS Data**: Patterns are synthetic and may not reflect actual drug-event associations
- **Simplified MedDRA**: Uses MedDRA-like structure, not official licensed dictionary
- **Perfect for Development**: Designed for software development, not regulatory submission

### Citation

If you use this dataset, please cite:
Synthetic FAERS-Style Pharmacovigilance Dataset v2.0 (2025)
Generated: 2025-12-17
Records: 1,000,000
Privacy: GDPR/HIPAA Compliant


### License

Public Domain (CC0) - Fully synthetic data with no privacy restrictions

### Contact

For questions or issues, please open a discussion on this dataset's page.

---
**Disclaimer**: This is entirely synthetic data created for machine learning and software development purposes. It does not contain any real patient information and should not be used for actual drug safety decisions or regulatory submissions.
