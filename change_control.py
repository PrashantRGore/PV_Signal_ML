# CIOMS Production Governance: Predetermined Change Control Plan
PCP_VERSION = "1.0"
CHANGE_CONTROL_THRESHOLDS = {
    "model_drift": 2.0,  # Z-score
    "fairness_violation": 0.2,  # Subgroup imbalance
    "prr_shift": 1.5,  # 50% PRR distribution shift
    "new_cases_threshold": 1000  # Quarterly case volume
}

CHANGE_ACTIONS = {
    "retrain": "Full retrain with new FAERS quarter + RAG context update",
    "review": "Human-in-loop validation of top 100 signals",
    "rollback": "Revert to previous MLflow run",
    "audit": "Full lineage + DPIA re-assessment"
}

# Auto-generate PCP compliance report
def generate_pcp_report(drift_score, fairness_report):
    issues = []
    if drift_score > CHANGE_CONTROL_THRESHOLDS["model_drift"]:
        issues.append("MODEL_DRIFT → RETRAIN")
    if not fairness_report['gender_equity']:
        issues.append("FAIRNESS → SUBGROUP AUDIT")
    
    with open('governance/pcp_compliance.md', 'w') as f:
        f.write(f"# Predetermined Change Control Plan v{PCP_VERSION}\n\n")
        f.write(f"**Generated:** {pd.Timestamp.now()}\n\n")
        f.write("## Required Actions:\n")
        for issue in issues:
            f.write(f"- {issue}\n")
        f.write("\n## Thresholds:\n``````")
    
    return issues
