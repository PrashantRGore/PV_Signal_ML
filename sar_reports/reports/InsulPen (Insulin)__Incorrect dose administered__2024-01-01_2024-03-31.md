# Signal Assessment Report – InsulPen (Insulin) / Incorrect dose administered

**Period assessed:** 2024/01/01:2024/03/31  
**Data sources:** Internal spontaneous DB / FAERS aggregate

## 1. Background

- Product: InsulPen (Insulin)
- Adverse event (MedDRA PT): Incorrect dose administered
- Detected via: disproportionality (PRR / chi-square) on spontaneous reports.
- Statistical thresholds: at least 3 cases, PRR ≥ 2.0, chi-square ≥ 4.0.

## 2. Quantitative evidence – current cases


- Cases (current dataset): 14455
- PRR (current dataset): 7.69
- Chi-square (current dataset): 71544.26


Interpretation (current cases): statistical criteria for a signal are met.

## 3. Quantitative evidence – FAERS period (2024/01/01:2024/03/31)


- No FAERS period aggregate available for this drug–event in the selected window.


## 4. Additional evidence (qualitative)


- **Regulatory guidance** – EMA GVP Module IX – Signal management (public)  
  Signal detected using routine disproportionality; requires clinical evaluation and consideration of risk minimisation measures in line with GVP Module IX.

- **Related signal** – Internal disproportionality analysis  
  Related signal: InsulPen (Insulin) – Incorrect dose administered (PRR=7.69, CASES=14455).

- **Related signal** – Internal disproportionality analysis  
  Related signal: PainAway (Oxycodone) – Accidental Overdose (PRR=9.29, CASES=17612).

- **Related signal** – Internal disproportionality analysis  
  Related signal: CardioFlow (Atenolol) – Bradycardia (PRR=12.62, CASES=24024).


## 5. Assessment and recommendation

**Signal status:** validated statistical signal

**Proposed actions for the MAH:**

- Include this drug–event pair in routine signal review meetings.

- Monitor new spontaneous reports and literature for this association.

- No immediate label changes are proposed at this stage; further assessment should consider clinical relevance, seriousness, and alternative explanations.

**Justification:**  
PRR=7.69 with 14455 cases exceeds predefined thresholds; clinical relevance must be assessed.

## 6. Methodological note

This assessment follows the principles of GVP Module IX on signal management and uses routine disproportionality methods on spontaneous reports (e.g. FAERS and internal MAH databases) complemented by literature and regulatory guidance.[web:68][web:80]  
Methods parameters (thresholds, input files) are captured in the model run metadata and methods artifact (for example, `methods_f7217712a96947589800c2bf33438600.json` where available).

## 7. Limitations

- Spontaneous reporting systems such as FAERS and internal MAH safety databases are subject to under-reporting and various reporting biases.
- Disproportionality signals do not establish causality; confounding by indication, co-medications, and data quality issues may contribute.
- Quantitative strength should always be interpreted together with clinical review, case narratives, and external evidence (literature and regulatory assessments).

## 8. Privacy and data protection

- Individual case data are processed in de-identified form and aggregated before signal detection.
- No direct identifiers or full narratives are exposed to the reporting or RAG layer.
- Current model: XGBoost v1.0, DP enabled: False.
- MLflow tracking: run_id f7217712a96947589800c2bf33438600, experiment pv-signal-ml (local file-based tracking).[web:126][web:135]
- A future configuration will allow training with differentially private SGD (DP-SGD) in line with published guidance to further bound the influence of any single ICSR.[web:89][web:95][web:106][web:109][web:112]