# Signal Assessment Report – {{ drug }} / {{ event }}

**Period assessed:** {{ period }}  
**Data sources:** {{ data_sources | join(', ') }}

## 1. Background

- Product: {{ drug }}
- Adverse event (MedDRA PT): {{ event }}
- Detected via: disproportionality (PRR / chi-square) on spontaneous reports.
- Statistical thresholds: at least {{ statistics.thresholds.min_cases }} cases, PRR ≥ {{ statistics.thresholds.min_prr }}, chi-square ≥ {{ statistics.thresholds.min_chisq }}.

## 2. Quantitative evidence – current cases

{% if global_statistics %}
- Cases (current dataset): {{ global_statistics.cases }}
- PRR (current dataset): {{ '%.2f'|format(global_statistics.prr) }}
- Chi-square (current dataset): {{ '%.2f'|format(global_statistics.chisq) }}
{% else %}
- No current-dataset statistics available for this drug–event.
{% endif %}

Interpretation (current cases): {% if global_statistics and global_statistics.prr >= statistics.thresholds.min_prr and global_statistics.cases >= statistics.thresholds.min_cases %}statistical criteria for a signal are met.{% else %}statistical criteria are not fully met; treat as a topic for monitoring only.{% endif %}

## 3. Quantitative evidence – FAERS period ({{ period }})

{% if period_statistics %}
- Cases (FAERS period): {{ period_statistics.cases }}
- PRR (FAERS period): {{ '%.2f'|format(period_statistics.prr) }}
- Chi-square (FAERS period): {{ '%.2f'|format(period_statistics.chisq) }}
{% else %}
- No FAERS period aggregate available for this drug–event in the selected window.
{% endif %}

## 4. Additional evidence (qualitative)

{% for item in evidence_items %}
- **{{ item.type }}** – {{ item.source }}  
  {{ item.summary }}
{% endfor %}

## 5. Assessment and recommendation

**Signal status:** {{ recommendation.signal_status }}

**Proposed actions for the MAH:**
{% for act in recommendation.actions %}
- {{ act }}
{% endfor %}
- No immediate label changes are proposed at this stage; further assessment should consider clinical relevance, seriousness, and alternative explanations.

**Justification:**  
{{ recommendation.justification }}

## 6. Methodological note

This assessment follows the principles of GVP Module IX on signal management and uses routine disproportionality methods on spontaneous reports (e.g. FAERS and internal MAH databases) complemented by literature and regulatory guidance.[web:68][web:80]  
Methods parameters (thresholds, input files) are captured in the model run metadata and methods artifact (for example, `methods_{{ model_info.run_id }}.json` where available).

## 7. Limitations

- Spontaneous reporting systems such as FAERS and internal MAH safety databases are subject to under-reporting and various reporting biases.
- Disproportionality signals do not establish causality; confounding by indication, co-medications, and data quality issues may contribute.
- Quantitative strength should always be interpreted together with clinical review, case narratives, and external evidence (literature and regulatory assessments).

## 8. Privacy and data protection

- Individual case data are processed in de-identified form and aggregated before signal detection.
- No direct identifiers or full narratives are exposed to the reporting or RAG layer.
- Current model: {{ model_info.algorithm }} v{{ model_info.version }}, DP enabled: {{ model_info.dp_enabled }}{% if model_info.epsilon %}, epsilon ≈ {{ model_info.epsilon }}{% endif %}.
- MLflow tracking: run_id {{ model_info.run_id }}, experiment {{ model_info.experiment_name }} (local file-based tracking).[web:126][web:135]
- A future configuration will allow training with differentially private SGD (DP-SGD) in line with published guidance to further bound the influence of any single ICSR.[web:89][web:95][web:106][web:109][web:112]
