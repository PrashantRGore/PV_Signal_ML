import os
import pickle
import subprocess
import time
import socket
from pathlib import Path
from datetime import date

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px
import networkx as nx
import requests

# --------------------------------------------------------------------
# Paths and constants
# --------------------------------------------------------------------
BASE_DIR = Path(r"C:\Users\koreo\Downloads\pv-signal-ml")
SAR_DIR = BASE_DIR / "sar_reports"
REPORTS_DIR = SAR_DIR / "reports"
GRAPH_PATH = SAR_DIR / "signal_graph.pkl"
SHAP_IMG_PATH = SAR_DIR / "shap_production.png"
SIGNALS_PATH = SAR_DIR / "full_signals_1M.csv"

API_URL = "http://127.0.0.1:8000/signal-report"  # FastAPI service
API_HOST = "127.0.0.1"
API_PORT = 8000

# --------------------------------------------------------------------
# Helper Functions
# --------------------------------------------------------------------

def is_port_open(host, port, timeout=1):
    """Check if API port is open."""
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.settimeout(timeout)
    result = sock.connect_ex((host, port))
    sock.close()
    return result == 0


def start_api_server():
    """Start FastAPI server in background."""
    if is_port_open(API_HOST, API_PORT):
        return True  # Already running
    
    try:
        # Start API in background
        subprocess.Popen(
            ["python", "-m", "uvicorn", "api:app", f"--host", API_HOST, f"--port", str(API_PORT)],
            cwd=str(BASE_DIR),
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )
        
        # Wait for API to start
        for _ in range(30):  # Wait up to 30 seconds
            time.sleep(1)
            if is_port_open(API_HOST, API_PORT):
                return True
        
        return False
    except Exception as e:
        st.error(f"Failed to start API: {e}")
        return False


# Initialize session state
if 'api_started' not in st.session_state:
    st.session_state.api_started = False

# Try to start API if not already running
if not st.session_state.api_started:
    if is_port_open(API_HOST, API_PORT):
        st.session_state.api_started = True
    else:
        with st.spinner("Starting API server..."):
            if start_api_server():
                st.session_state.api_started = True
                st.success("✅ API server started successfully!")
            else:
                st.warning("⚠️ Could not start API server. Please ensure port 8000 is available.")

# --------------------------------------------------------------------
# Streamlit page config
# --------------------------------------------------------------------
st.set_page_config(
    page_title="PV-SIGNAL-ML Full GxP Stack",
    layout="wide",
)

st.title("🚀 PV-SIGNAL-ML – Full GxP Stack")
tabs = st.tabs(["Dashboard", "PSMF Annex Generator", "Methods & Governance"])

# ========================= TAB 1: DASHBOARD =========================
with tabs[0]:
    st.markdown(
        "**1M ICSR → PRR → XGBoost → SHAP → GraphRAG → EMA‑style Signal Reports**"
    )

    # Load signals (current dataset)
    if not SIGNALS_PATH.exists():
        st.error(f"Signals file not found: {SIGNALS_PATH}")
        st.stop()

    signals = pd.read_csv(SIGNALS_PATH)

    # Basic KPIs
    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Signals (PRR>2.0)", len(signals))
    col2.metric("Top PRR", f"{signals.PRR.max():.2f}")
    col3.metric("Top signal drug", signals.DRUG.iloc[0])
    col4.metric("Top signal event", signals.EVENT.iloc[0])

    # Top signals table
    st.subheader("🏆 Top PRR Signals (current dataset)")
    st.dataframe(
        signals[["DRUG", "EVENT", "PRR", "CASES"]].head(20),
        use_container_width=True,
    )

    # PRR vs CHISQ scatter
    if "CHISQ" in signals.columns:
        st.subheader("Signal Prioritisation Matrix (PRR vs CHISQ)")
        fig_scatter = px.scatter(
            signals,
            x="PRR",
            y="CHISQ",
            size="CASES",
            hover_data=["DRUG", "EVENT"],
            title="PRR vs CHISQ (bubble size = CASES)",
        )
        st.plotly_chart(fig_scatter, use_container_width=True)

    # SHAP image
    if SHAP_IMG_PATH.exists():
        st.subheader("SHAP Explainability")
        st.image(str(SHAP_IMG_PATH), use_container_width=True)

    # GraphRAG network
    st.subheader("GraphRAG Network (Drug–Event Graph)")
    if GRAPH_PATH.exists():
        with open(GRAPH_PATH, "rb") as f:
            G = pickle.load(f)

        fig_graph, ax = plt.subplots(figsize=(8, 6))
        pos = nx.spring_layout(G, seed=42)
        nx.draw_networkx_nodes(G, pos, node_color="lightblue", node_size=800, ax=ax)
        nx.draw_networkx_edges(G, pos, alpha=0.4, ax=ax)
        nx.draw_networkx_labels(G, pos, font_size=7, ax=ax)
        ax.set_axis_off()
        st.pyplot(fig_graph)
    else:
        st.info("Graph file not found, skipping network plot.")

    # EMA-style signal report generator
    st.subheader("📄 EMA‑style Signal Assessment Report Generator")

    # Period selection
    col_p1, col_p2 = st.columns(2)
    default_start = date(2024, 1, 1)
    default_end = date(2024, 3, 31)

    start_date = col_p1.date_input("Start date", value=default_start)
    end_date = col_p2.date_input("End date", value=default_end)

    if start_date > end_date:
        st.error("Start date must be on or before end date.")
        st.stop()

    period_str = f"{start_date.strftime('%Y/%m/%d')}:{end_date.strftime('%Y/%m/%d')}"
    st.write(f"**Selected assessment period:** {period_str}")

    # Select a signal row
    selected_idx = st.number_input(
        "Row index from top signals table (0–9)",
        min_value=0,
        max_value=min(9, len(signals) - 1),
        value=0,
        step=1,
    )

    sel = signals.iloc[int(selected_idx)]
    st.write(f"**Selected signal:** {sel.DRUG} → {sel.EVENT}")

    # Generate report
    if st.button("Generate signal assessment report"):
        payload = {
            "drug": sel.DRUG,
            "event": sel.EVENT,
            "period": period_str,
        }
        st.write("Calling report API...")
        try:
            resp = requests.post(API_URL, json=payload, timeout=120)
        except Exception as e:
            st.error(f"API call failed: {e}")
        else:
            if resp.status_code == 200:
                report = resp.json()
                safe_period = period_str.replace("/", "-").replace(":", "_")

                st.success(
                    "Report generated and stored under "
                    f"`sar_reports/reports/{sel.DRUG}__{sel.EVENT}__{safe_period}.json/md`."
                )

                st.markdown("### JSON evidence")
                st.json(report)

                # Quantitative comparison: current dataset vs FAERS period
                gs = report.get("global_statistics")
                ps = report.get("period_statistics")
                if gs or ps:
                    st.markdown("### Quantitative comparison (current vs FAERS period)")
                    rows = []
                    if gs:
                        rows.append(
                            {
                                "Source": "Current dataset (full_signals_1M)",
                                "Cases": gs["cases"],
                                "PRR": gs["prr"],
                                "Chi-square": gs["chisq"],
                            }
                        )
                    if ps:
                        rows.append(
                            {
                                "Source": f"FAERS {period_str}",
                                "Cases": ps["cases"],
                                "PRR": ps["prr"],
                                "Chi-square": ps["chisq"],
                            }
                        )
                    if rows:
                        comp_df = pd.DataFrame(rows)
                        st.dataframe(comp_df, use_container_width=True)

                # Markdown report + download + bundle hint
                md_path = REPORTS_DIR / f"{sel.DRUG}__{sel.EVENT}__{safe_period}.md"
                if md_path.exists():
                    st.markdown("### Markdown report")
                    with open(md_path, "r", encoding="utf-8") as f:
                        md_text = f.read()
                    st.markdown(md_text)

                    md_bytes = md_text.encode("utf-8")
                    st.download_button(
                        "Download report (.md)",
                        data=md_bytes,
                        file_name=md_path.name,
                        mime="text/markdown",
                    )

                    st.info(
                        "For archival or sharing, you can also run in a terminal:\n\n"
                        f"`python export_assessment_bundle.py \"{sel.DRUG}\" "
                        f"\"{sel.EVENT}\" \"{period_str}\"`"
                    )
                else:
                    st.info("Markdown file not found on disk.")
            else:
                st.error(f"API error: {resp.status_code} – {resp.text}")

    st.markdown(
        """
---
✅ Data: SQLite + Parquet • ✅ Stats: PRR/chi-square • ✅ ML: XGBoost + SHAP  
✅ Context: GraphRAG + RAG Reports • ✅ UI: Streamlit (audit‑oriented)
"""
    )

# ====================== TAB 2: PSMF ANNEX GENERATOR =================
with tabs[1]:
    st.header("📋 PSMF Annex D Generator")
    st.markdown("""
    Generate regulatory-compliant PSMF Annex D documentation describing your signal detection system.
    
    **PSMF Annex D** includes:
    - Signal detection methodology (PRR, Chi-square, ML triage)
    - Computerized system description (architecture, components)
    - Products monitored (extracted from signals database)
    - Quality assurance procedures
    - Regulatory compliance mapping
    """)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("📄 Generate PSMF Annex D")
        st.write("Click the button below to generate a comprehensive PSMF Annex D document.")
        
        if st.button("🔄 Generate PSMF Annex D", key="generate_psmf"):
            with st.spinner("Generating PSMF Annex D..."):
                try:
                    from psmf_annex_generator import generate_psmf_annex_d
                    result = generate_psmf_annex_d()
                    
                    if result['status'] == 'success':
                        st.success("✅ PSMF Annex D generated successfully!")
                        st.json({
                            'File': result['file_path'],
                            'Products Monitored': result['products_monitored'],
                            'Generated At': result['timestamp'],
                            'Content Size': f"{result['content_length']} characters"
                        })
                        
                        # Read and display the generated file
                        with open(result['file_path'], 'r', encoding='utf-8') as f:
                            psmf_content = f.read()
                        
                        st.markdown("### 📖 Generated Document Preview")
                        st.markdown(psmf_content[:2000] + "\n\n...[truncated for display]...")
                        
                        # Download button
                        psmf_bytes = psmf_content.encode('utf-8')
                        st.download_button(
                            "📥 Download PSMF Annex D (.md)",
                            data=psmf_bytes,
                            file_name=Path(result['file_path']).name,
                            mime="text/markdown"
                        )
                    else:
                        st.error(f"Failed to generate PSMF Annex D: {result.get('error', 'Unknown error')}")
                        
                except Exception as e:
                    st.error(f"❌ Error generating PSMF Annex D: {str(e)}")
                    st.write("**Troubleshooting:**")
                    st.write("1. Ensure `psmf_annex_generator.py` is in the project directory")
                    st.write("2. Check that `full_signals_1M.csv` exists in `sar_reports/`")
                    st.write("3. Verify all dependencies are installed")
    
    with col2:
        st.subheader("ℹ️ About PSMF Annex D")
        st.markdown("""
        **PSMF** (Periodic Safety Update Format) is a regulatory document required for:
        - **EMA:** Annual safety reports
        - **FDA:** Periodic safety updates
        - **CIOMS XIV:** Signal management procedures
        
        **Annex D** specifically covers:
        - Signal detection methodology
        - Computerized systems
        - Products monitored
        - Quality assurance
        
        **Note:** This is a technical annex. A complete PSMF requires:
        - Company information (QPPV, contacts)
        - Organizational procedures
        - Audit trail documentation
        - Regulatory history
        """)
        
        st.info("💡 **Tip:** Generate this document quarterly or when your signal detection system changes.")

# ====================== TAB 3: METHODS & GOVERNANCE =================
with tabs[2]:
    st.header("Methods & Governance")

    st.markdown("""
### What this system does

- Calculates disproportionality statistics (PRR, chi-square) on aggregated spontaneous reports.
- Trains an XGBoost model on structured features, with runs tracked in MLflow (local SQLite backend).
- Generates EMA-style signal assessment reports (JSON + Markdown) and optional ZIP bundles for archival.

### How PRR is computed

- For each drug–event pair, counts:
  - a: cases with the drug and the event.
  - b: cases with the drug and other events.
  - c: cases with the event and other drugs.
  - d: all other cases.
- Computes the Proportional Reporting Ratio as  
  PRR = [a / (a + b)] / [c / (c + d)], and a chi-square statistic for strength of association.[web:20]
- Applies basic thresholds (for example, at least 3 cases, PRR ≥ 2.0, chi-square ≥ 4.0) before flagging signals.[web:19]

### How reports are generated

- The API accepts `drug`, `event`, and `period`.
- For the current dataset, it pulls statistics from `full_signals_1M.csv`.
- For a specific period (e.g. 2025/01/01:2025/03/31), it uses a period-tagged signals file such as `full_signals_2025-01-01_2025-03-31.csv` when available.
- It enriches the quantitative view with:
  - Related internal signals via sentence-transformer embeddings.
  - Live PubMed literature snippets constrained to the selected period.
  - A templated EMA-style narrative, including recommendations, limitations, and privacy statements.[web:30][web:29]

### Controls and governance

- Data are de-identified and aggregated before modeling; no direct identifiers or full narratives are passed into RAG or reporting layers.[web:28]
- Model versions and parameters are tracked in MLflow (`sqlite:///mlflow.db`), enabling full run reproducibility.
- Methods, thresholds, and input files per run are captured in `sar_reports/methods_<run_id>.json`.
- Signal reports can be exported as ZIP bundles (JSON + Markdown + methods) for audit, review meetings, or regulatory interactions.
""")