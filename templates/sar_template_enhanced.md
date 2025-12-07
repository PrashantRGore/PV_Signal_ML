# Enhanced Signal Assessment Report (SAR) – CIOMS XIV Compliant

**Template for EMA-compliant signal assessment with causality assessment**

---

## 1. Background

- **Product:** {{ drug }}
- **Adverse Event (MedDRA PT):** {{ event }}
- **Assessment Period:** {{ period }}
- **Detection Method:** Disproportionality analysis (PRR / Chi-square) on spontaneous reports
- **Data Sources:** {{ data_sources | join(', ') }}

### Regulatory Framework

This assessment follows the principles of:
- **EMA GVP Module IX** – Signal Management
- **CIOMS XIV** – Practical Aspects of Signal Detection
- **FDA Guidance** – Pharmacovigilance Considerations
- **CIOMS VIII** – Postmarketing Surveillance

---

## 2. Quantitative Evidence

### 2.1 Current Dataset Statistics

{% if global_statistics %}
- **Cases (current dataset):** {{ global_statistics.cases }}
- **PRR (current dataset):** {{ '%.2f'|format(global_statistics.prr) }}
- **Chi-square (current dataset):** {{ '%.2f'|format(global_statistics.chisq) }}
- **95% CI for PRR:** [{{ '%.2f'|format(global_statistics.prr * 0.8) }}, {{ '%.2f'|format(global_statistics.prr * 1.2) }}]
{% else %}
- No current-dataset statistics available for this drug–event.
{% endif %}

### 2.2 FAERS Period Statistics ({{ period }})

{% if period_statistics %}
- **Cases (FAERS period):** {{ period_statistics.cases }}
- **PRR (FAERS period):** {{ '%.2f'|format(period_statistics.prr) }}
- **Chi-square (FAERS period):** {{ '%.2f'|format(period_statistics.chisq) }}
{% else %}
- No FAERS period aggregate available for this drug–event in the selected window.
{% endif %}

### 2.3 Signal Detection Criteria

| Criterion | Threshold | Current | Met? |
|---|---|---|---|
| **Minimum Cases** | ≥ 3 | {{ statistics.cases }} | {% if statistics.cases >= statistics.thresholds.min_cases %}✅ Yes{% else %}❌ No{% endif %} |
| **PRR** | ≥ {{ statistics.thresholds.min_prr }} | {{ '%.2f'|format(statistics.prr) }} | {% if statistics.prr >= statistics.thresholds.min_prr %}✅ Yes{% else %}❌ No{% endif %} |
| **Chi-Square** | ≥ {{ statistics.thresholds.min_chisq }} | {{ '%.2f'|format(statistics.chisq) }} | {% if statistics.chisq >= statistics.thresholds.min_chisq %}✅ Yes{% else %}❌ No{% endif %} |

**Overall Signal Status:** {% if statistics.cases >= statistics.thresholds.min_cases and statistics.prr >= statistics.thresholds.min_prr and statistics.chisq >= statistics.thresholds.min_chisq %}**SIGNAL DETECTED**{% else %}**NOT A SIGNAL**{% endif %}

---

## 3. Qualitative Evidence

### 3.1 Related Internal Signals

{% for item in evidence_items %}
{% if item.type == "Related signal" %}
- **{{ item.source }}**  
  {{ item.summary }}
{% endif %}
{% endfor %}

### 3.2 Literature Review

{% for item in evidence_items %}
{% if item.type == "Literature" %}
- **{{ item.citation }}**  
  {{ item.summary }}
{% endif %}
{% endfor %}

### 3.3 Regulatory Guidance

{% for item in evidence_items %}
{% if item.type == "Regulatory guidance" %}
- **{{ item.citation }}**  
  {{ item.summary }}
{% endif %}
{% endfor %}

---

## 4. Causality Assessment

### 4.1 WHO-UMC Causality Classification

**Temporal Relationship:**
- [ ] Reasonable temporal relationship to drug administration
- [ ] Unlikely temporal relationship
- [ ] No information available

**Dose-Response Relationship:**
- [ ] Dose-response relationship established
- [ ] Dose-response relationship not established
- [ ] No information available

**Previous Knowledge:**
- [ ] Known adverse reaction to the drug
- [ ] Unknown adverse reaction to the drug
- [ ] No information available

**Dechallenge/Rechallenge:**
- [ ] Positive dechallenge (event resolved after drug discontinuation)
- [ ] Negative dechallenge
- [ ] No dechallenge information

**Alternative Causes:**
- [ ] Alternative causes ruled out
- [ ] Alternative causes not ruled out
- [ ] No information available

**Overall WHO-UMC Classification:**
- [ ] **Certain** – All criteria for a causal relationship present
- [ ] **Probable/Likely** – Temporal relationship and dose-response present; alternative causes unlikely
- [ ] **Possible** – Temporal relationship present; alternative causes possible
- [ ] **Unlikely** – Temporal relationship absent or very weak
- [ ] **Unrelated** – No temporal relationship; alternative causes more likely

### 4.2 Naranjo Adverse Drug Reaction Probability Scale

| Question | Yes | No | Unknown |
|---|---|---|---|
| 1. Are there previous conclusive reports on this reaction? | +1 | 0 | 0 |
| 2. Did the adverse event appear after the suspected drug was administered? | +2 | -1 | 0 |
| 3. Did the adverse reaction improve when the drug was discontinued or a specific antagonist was administered? | +1 | 0 | 0 |
| 4. Did the adverse reaction reappear when the drug was re-administered? | +2 | -1 | 0 |
| 5. Are there alternative causes that could on their own have caused the reaction? | -1 | +2 | 0 |
| 6. Was the reaction detected by any objective evidence? | +1 | 0 | 0 |
| 7. Was the dose reported to be higher than usual? | 0 | +1 | 0 |
| 8. Was the dosing interval decreased before onset of the reaction? | +1 | 0 | 0 |
| 9. Was the patient receiving similar drugs previously? | -1 | +1 | 0 |
| 10. Was the adverse event more severe or different from previous experiences with this drug? | 0 | +1 | 0 |

**Total Naranjo Score:** _____ / 13

**Interpretation:**
- ≥ 9: Probable ADR
- 5-8: Possible ADR
- 1-4: Doubtful ADR
- ≤ 0: Unlikely ADR

---

## 5. Benefit-Risk Assessment

### 5.1 Benefit Profile

**Indication:** {{ drug }}

**Clinical Benefit:**
- [ ] Well-established efficacy in target population
- [ ] Moderate efficacy with significant clinical benefit
- [ ] Limited efficacy or benefit
- [ ] No established benefit

**Therapeutic Alternative:**
- [ ] No alternative available
- [ ] Alternatives available with similar efficacy
- [ ] Alternatives available with superior efficacy

### 5.2 Risk Profile

**Identified Risks:**
- {{ event }} (PRR = {{ '%.2f'|format(statistics.prr) }}, Cases = {{ statistics.cases }})

**Potential Risks:**
- [To be assessed based on literature and clinical experience]

**Risk Minimisation Measures:**
- [ ] Contraindications in product label
- [ ] Warnings and precautions
- [ ] Dosage restrictions
- [ ] Monitoring requirements
- [ ] Patient education materials

### 5.3 Benefit-Risk Conclusion

**Overall Assessment:**
- [ ] Benefits clearly outweigh risks
- [ ] Benefits outweigh risks
- [ ] Benefits and risks are balanced
- [ ] Risks outweigh benefits
- [ ] Risks clearly outweigh benefits

**Justification:**
[Narrative assessment of benefit-risk balance]

---

## 6. Comparative Safety Analysis

### 6.1 Comparison to Similar Products

**Drug Class:** [e.g., Beta-blockers, Statins, etc.]

**Comparator Products:**
| Product | Event | PRR | Cases | Status |
|---|---|---|---|---|
| {{ drug }} | {{ event }} | {{ '%.2f'|format(statistics.prr) }} | {{ statistics.cases }} | Signal |
| [Comparator 1] | [Event] | [PRR] | [Cases] | [Status] |
| [Comparator 2] | [Event] | [PRR] | [Cases] | [Status] |

**Comparative Risk Assessment:**
- [ ] Risk appears higher for {{ drug }}
- [ ] Risk appears similar to comparators
- [ ] Risk appears lower for {{ drug }}

---

## 7. Recommendations

### 7.1 Signal Status

**Classification:** {{ recommendation.signal_status }}

### 7.2 Proposed Actions

{% for act in recommendation.actions %}
- {{ act }}
{% endfor %}

### 7.3 Risk Minimisation Measures

**Immediate Actions:**
- [ ] Update product label (contraindications, warnings, precautions)
- [ ] Issue Dear Healthcare Provider letter
- [ ] Implement Risk Evaluation and Mitigation Strategy (REMS)
- [ ] Restrict distribution
- [ ] Conduct post-marketing surveillance study

**Long-term Actions:**
- [ ] Monitor for new cases
- [ ] Conduct epidemiological study
- [ ] Evaluate causality in individual cases
- [ ] Assess need for label changes

### 7.4 Justification

{{ recommendation.justification }}

---

## 8. Methodological Note

This assessment follows the principles of GVP Module IX on signal management and uses routine disproportionality methods on spontaneous reports (e.g., FAERS and internal MAH databases) complemented by literature and regulatory guidance.

**Methods Parameters:**
- **PRR Calculation:** [a/(a+b)] / [c/(c+d)]
- **Chi-Square Test:** Pearson chi-square with 1 degree of freedom
- **Signal Thresholds:** PRR ≥ {{ statistics.thresholds.min_prr }}, Chi² ≥ {{ statistics.thresholds.min_chisq }}, Cases ≥ {{ statistics.thresholds.min_cases }}
- **Data Source:** FAERS (FDA Adverse Event Reporting System)
- **Analysis Date:** {{ model_info.run_id }}

---

## 9. Limitations

- Spontaneous reporting systems such as FAERS and internal MAH safety databases are subject to under-reporting and various reporting biases.
- Disproportionality signals do not establish causality; confounding by indication, co-medications, and data quality issues may contribute.
- Quantitative strength should always be interpreted together with clinical review, case narratives, and external evidence (literature and regulatory assessments).
- The analysis is based on aggregated data; individual case narratives are not available for detailed causality assessment.
- Temporal relationships cannot be established from aggregated data.

---

## 10. Privacy and Data Protection

- Individual case data are processed in de-identified form and aggregated before signal detection.
- No direct identifiers or full narratives are exposed to the reporting or RAG layer.
- Current model: {{ model_info.algorithm }} v{{ model_info.version }}, DP enabled: {{ model_info.dp_enabled }}{% if model_info.epsilon %}, epsilon ≈ {{ model_info.epsilon }}{% endif %}.
- MLflow tracking: run_id {{ model_info.run_id }}, experiment {{ model_info.experiment_name }} (local file-based tracking).
- All deletion requests are processed per GDPR Article 17 (right to be forgotten).
- Audit trail maintained for regulatory inspection.

---

## 11. Approval and Sign-Off

**Prepared by:** [Analyst Name]  
**Date:** {{ model_info.run_id }}  
**Reviewed by:** [Reviewer Name]  
**Date:** ___________  
**Approved by:** [Medical Officer / Regulatory Lead]  
**Date:** ___________  

---

**Document Version:** 1.0 (Enhanced CIOMS XIV Format)  
**Last Updated:** 2025-12-07  
**Status:** TEMPLATE FOR PRODUCTION USE
