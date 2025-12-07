# Signal Assessment Report – CardioFlow (Atenolol) / Bradycardia

**Period assessed:** FAERS-ALL  
**Data sources:** Internal spontaneous DB / FAERS aggregate

## 1. Background

- Product: CardioFlow (Atenolol)
- Adverse event (MedDRA PT): Bradycardia
- Detected via: disproportionality (PRR / chi-square) on spontaneous reports.
- Statistical thresholds: at least 3 cases, PRR ≥ 2.0, chi-square ≥ 4.0.

## 2. Initial evidence (quantitative)

- Cases (a): 24024
- Non-event cases on drug (b): 0
- Event cases on other drugs (c): 0
- Other combinations (d): 0
- PRR: 12.62
- Chi-square: 194822.79

Interpretation: statistical criteria for a signal are met.

## 3. Additional evidence (qualitative)


- **Regulatory guidance** – EMA GVP Module IX – Signal management (public)  
  Signal detected using routine disproportionality; requires clinical evaluation and consideration of risk minimisation measures in line with GVP Module IX.

- **Related signal** – Internal disproportionality analysis  
  Related signal: CardioFlow (Atenolol) – Bradycardia (PRR=12.62, CASES=24024).

- **Related signal** – Internal disproportionality analysis  
  Related signal: ChemoGold (Doxorubicin) – Cardiomyopathy (PRR=10.94, CASES=20792).

- **Related signal** – Internal disproportionality analysis  
  Related signal: StatShield (Atorvastatin) – Rhabdomyolysis (PRR=5.27, CASES=9866).


## 4. Assessment and recommendation

**Signal status:** validated statistical signal

**Proposed actions for the MAH:**

- Include this drug–event pair in routine signal review meetings.

- Monitor new spontaneous reports and literature for this association.


**Justification:**  
PRR=12.62 with 24024 cases exceeds predefined thresholds; clinical relevance must be assessed.

## 5. Methodological note

This assessment follows the principles of GVP Module IX on signal management and uses routine disproportionality methods on spontaneous reports complemented by literature and regulatory guidance.[web:68][web:80]

## 6. Privacy and data protection

- Individual case data are processed in de-identified form and aggregated before signal detection.
- No direct identifiers or full narratives are exposed to the reporting or RAG layer.
- Current model: XGBoost v1.0, DP enabled: False.
- A future configuration will allow training with differentially private SGD (DP-SGD) in line with published guidance to further bound the influence of any single ICSR.[web:89][web:95][web:106][web:109][web:112]