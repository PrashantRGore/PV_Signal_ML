from pathlib import Path
import json
import zipfile
import sys

from signal_report_builder import build_signal_report, SAR_DIR
from train_with_mlflow import BASE_DIR

REPORTS_DIR = SAR_DIR / "reports"
BUNDLES_DIR = BASE_DIR / "bundles"
BUNDLES_DIR.mkdir(exist_ok=True, parents=True)


def export_bundle(drug: str, event: str, period: str = "FAERS-ALL"):
    """
    Build signal report (if needed) and create a ZIP bundle containing:
    - JSON report
    - Markdown report
    - methods_{run_id}.json (if available)
    """
    # Ensure report exists
    build_signal_report(drug, event, period)

    # Sanitize drug, event, and period for safe filenames
    safe_drug = str(drug).replace('/', '_').replace('\\', '_').replace(':', '_').replace('*', '_').replace('?', '_').replace('"', '_').replace('<', '_').replace('>', '_').replace('|', '_')
    safe_event = str(event).replace('/', '_').replace('\\', '_').replace(':', '_').replace('*', '_').replace('?', '_').replace('"', '_').replace('<', '_').replace('>', '_').replace('|', '_')
    safe_period = period.replace("/", "-").replace(":", "_")
    
    json_path = REPORTS_DIR / f"{safe_drug}__{safe_event}__{safe_period}.json"
    md_path = REPORTS_DIR / f"{safe_drug}__{safe_event}__{safe_period}.md"

    if not json_path.exists():
        raise FileNotFoundError(json_path)

    with open(json_path, "r", encoding="utf-8") as f:
        report = json.load(f)

    run_id = report.get("model_info", {}).get("run_id")
    methods_path = SAR_DIR / f"methods_{run_id}.json" if run_id else None

    bundle_name = f"{safe_drug}__{safe_event}__{safe_period}.zip"
    bundle_path = BUNDLES_DIR / bundle_name

    with zipfile.ZipFile(bundle_path, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.write(json_path, arcname=json_path.name)
        if md_path.exists():
            zf.write(md_path, arcname=md_path.name)
        if methods_path and methods_path.exists():
            zf.write(methods_path, arcname=methods_path.name)

    print(f"✅ Bundle created: {bundle_path}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print('Usage: python export_assessment_bundle.py "DRUG" "EVENT" [PERIOD]')
        sys.exit(1)
    drug = sys.argv[1]
    event = sys.argv[2]
    period = sys.argv[3] if len(sys.argv) > 3 else "FAERS-ALL"
    export_bundle(drug, event, period)
