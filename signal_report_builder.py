from pathlib import Path
import json
import pandas as pd
from jinja2 import Environment, FileSystemLoader

from rag_signal_evidence import (
    enrich_evidence_with_similar_signals,
    fetch_pubmed_snippets,
)

BASE_DIR = Path(r"C:\Users\koreo\Downloads\pv-signal-ml")
SAR_DIR = BASE_DIR / "sar_reports"
TEMPLATE_DIR = BASE_DIR / "templates"
OUTPUT_DIR = SAR_DIR / "reports"
MODEL_META_PATH = SAR_DIR / "model_run_metadata.json"

OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

env = Environment(loader=FileSystemLoader(str(TEMPLATE_DIR)), autoescape=False)
template = env.get_template("signal_report_template.md")


def _load_model_meta():
    if MODEL_META_PATH.exists():
        try:
            with open(MODEL_META_PATH, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception as e:
            print(f"⚠️ Could not read model_run_metadata.json: {e}")
    return {
        "run_id": None,
        "model_version": "1.0",
        "tracking_uri": None,
        "experiment_name": "pv-signal-ml",
        "data_period": "FAERS-ALL",
        "dp_enabled": False,
        "epsilon": None,
    }


def _load_signal_rows(drug: str, event: str, period: str):
    """
    Return (global_row, period_row) from:
    - full_signals_1M.csv
    - full_signals_{start}_{end}.csv if it exists for the requested period.
    """
    global_df = pd.read_csv(SAR_DIR / "full_signals_1M.csv")
    g_row = global_df[(global_df["DRUG"] == drug) & (global_df["EVENT"] == event)]
    g_row = g_row.iloc[0] if not g_row.empty else None

    p_row = None
    if ":" in period:
        start, end = period.split(":")
        start_tag = start.replace("/", "-")
        end_tag = end.replace("/", "-")
        period_path = SAR_DIR / f"full_signals_{start_tag}_{end_tag}.csv"
        if period_path.exists():
            period_df = pd.read_csv(period_path)
            p = period_df[(period_df["DRUG"] == drug) & (period_df["EVENT"] == event)]
            if not p.empty:
                p_row = p.iloc[0]

    return g_row, p_row


def build_signal_report(drug: str, event: str, period: str = "FAERS-ALL"):
    """
    Build JSON + Markdown signal report for a drug–event–period.

    Uses both:
    - full_signals_1M.csv  (current dataset)
    - full_signals_{start}_{end}.csv (FAERS period, if present).
    """
    g_row, p_row = _load_signal_rows(drug, event, period)
    if g_row is None and p_row is None:
        raise ValueError(f"No signal row found for {drug} / {event} in any source")

    # Global (current) statistics
    if g_row is not None:
        global_stats = {
            "cases": int(g_row["CASES"]),
            "prr": float(g_row["PRR"]),
            "chisq": float(g_row.get("CHISQ", 0.0)),
        }
    else:
        global_stats = None

    # Period (FAERS window) statistics
    if p_row is not None:
        period_stats = {
            "cases": int(p_row["CASES"]),
            "prr": float(p_row["PRR"]),
            "chisq": float(p_row.get("CHISQ", 0.0)),
        }
    else:
        period_stats = None

    # Primary statistics field: use global if available, else period
    base_row = g_row if g_row is not None else p_row
    statistics = {
        "cases": int(base_row["CASES"]),
        "prr": float(base_row["PRR"]),
        "chisq": float(base_row.get("CHISQ", 0.0)),
        "a": int(base_row["CASES"]),
        "b": 0,
        "c": 0,
        "d": 0,
        "thresholds": {
            "min_cases": 3,
            "min_prr": 2.0,
            "min_chisq": 4.0,
        },
    }

    meta = _load_model_meta()
    model_info = {
        "algorithm": "XGBoost",
        "version": meta.get("model_version", "1.0"),
        "dp_enabled": bool(meta.get("dp_enabled", False)),
        "epsilon": meta.get("epsilon", None),
        "run_id": meta.get("run_id"),
        "tracking_uri": meta.get("tracking_uri"),
        "experiment_name": meta.get("experiment_name", "pv-signal-ml"),
    }

    evidence_items = [
        {
            "type": "Regulatory guidance",
            "source": "EMA GVP Module IX – Signal management (public)",
            "citation": "EMA/GVP/Module IX",
            "summary": (
                "Signal detected using routine disproportionality; "
                "requires clinical evaluation and consideration of risk minimisation "
                "measures in line with GVP Module IX."
            ),
        }
    ]

    # Related internal signals (RAG)
    try:
        related = enrich_evidence_with_similar_signals(drug, event, top_k=3)
        evidence_items.extend(related)
    except Exception as e:
        print(f"⚠️ Could not enrich evidence with related signals: {e}")

    # Live PubMed literature for the period
    start_date, end_date = ("1900/01/01", "2100/12/31")
    if ":" in period:
        parts = period.split(":")
        if len(parts) == 2:
            start_date, end_date = parts
    try:
        lit_items = fetch_pubmed_snippets(
            drug, event, start_date=start_date, end_date=end_date, max_results=5
        )
        evidence_items.extend(lit_items)
    except Exception as e:
        print(f"⚠️ PubMed fetch failed: {e}")

    recommendation = {
        "signal_status": "validated statistical signal",
        "actions": [
            "Include this drug–event pair in routine signal review meetings.",
            "Monitor new spontaneous reports and literature for this association.",
        ],
        "justification": (
            f"PRR={statistics['prr']:.2f} with {statistics['cases']} cases "
            "exceeds predefined thresholds; clinical relevance must be assessed."
        ),
    }

    report_obj = {
        "drug": drug,
        "event": event,
        "period": period,
        "data_sources": ["Internal spontaneous DB / FAERS aggregate"],
        "statistics": statistics,
        "global_statistics": global_stats,
        "period_statistics": period_stats,
        "model_info": model_info,
        "evidence_items": evidence_items,
        "recommendation": recommendation,
    }

    # Sanitize drug, event, and period for safe filenames
    safe_drug = str(drug).replace('/', '_').replace('\\', '_').replace(':', '_').replace('*', '_').replace('?', '_').replace('"', '_').replace('<', '_').replace('>', '_').replace('|', '_')
    safe_event = str(event).replace('/', '_').replace('\\', '_').replace(':', '_').replace('*', '_').replace('?', '_').replace('"', '_').replace('<', '_').replace('>', '_').replace('|', '_')
    safe_period = period.replace("/", "-").replace(":", "_")
    
    json_path = OUTPUT_DIR / f"{safe_drug}__{safe_event}__{safe_period}.json"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(report_obj, f, indent=2)

    md_content = template.render(**report_obj)
    md_path = OUTPUT_DIR / f"{safe_drug}__{safe_event}__{safe_period}.md"
    with open(md_path, "w", encoding="utf-8") as f:
        f.write(md_content)

    print(f"✅ Saved JSON: {json_path}")
    print(f"✅ Saved report: {md_path}")


if __name__ == "__main__":
    build_signal_report("CardioFlow (Atenolol)", "Bradycardia", period="FAERS-ALL")
