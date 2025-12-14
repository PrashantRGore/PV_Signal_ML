"""
PV Signal ML Enhanced System - PRODUCTION VERSION
Complete MLflow audit trail for FDA/EMA/PMDA compliance
Working state: 10:31 PM December 14, 2025
"""

import streamlit as st
import pandas as pd
import mlflow
from pathlib import Path
import config
from src.data.data_source_manager import DataSourceManager
from src.stats_engine.disproportionality import DisproportionalityAnalysis
from src.ml.sisa_trainer import SISATrainer
from src.explainability.shap_analysis import SHAPAnalyzer
from src.utils.logger import setup_logger

logger = setup_logger(__name__)

# ===== PAGE CONFIG =====
st.set_page_config(
    page_title="PV Signal ML - Regulatory Compliant",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ===== DISCLAIMER =====
st.warning(config.DISCLAIMER_TEXT)

# ===== TITLE =====
st.title("🔬 Pharmacovigilance Signal Detection with SISA & SHAP")
st.markdown("Enhanced System with Machine Unlearning, Explainable AI, RAG-based SAR, and MLflow")

# ===== SESSION STATE INITIALIZATION =====
if 'data_loaded' not in st.session_state:
    st.session_state.data_loaded = False
if 'signals' not in st.session_state:
    st.session_state.signals = None
if 'trained_model' not in st.session_state:
    st.session_state.trained_model = None
if 'current_sar' not in st.session_state:
    st.session_state.current_sar = None
if 'current_signal_data' not in st.session_state:
    st.session_state.current_signal_data = None
if 'current_model_choice' not in st.session_state:
    st.session_state.current_model_choice = None

# ===== SIDEBAR: SERVICE STATUS =====
st.sidebar.markdown("### 🔌 Service Status")

try:
    from src.rag.sar_generator import SARGenerator
    sar_gen = SARGenerator()
    ollama_ok, ollama_msg = sar_gen.check_service()
    if ollama_ok:
        st.sidebar.success(f"🟢 Ollama: {ollama_msg}")
    else:
        st.sidebar.error(f"🔴 Ollama: {ollama_msg}")
except:
    st.sidebar.warning("⚪ Ollama: Not checked")

try:
    mlflow.set_tracking_uri(config.MLFLOW_TRACKING_URI)
    experiments = mlflow.search_experiments()
    st.sidebar.success(f"🟢 MLflow: Connected ({len(experiments)} experiments)")
except:
    st.sidebar.error("🔴 MLflow: Connection failed")

st.sidebar.markdown("---")

# ===== DATA SOURCE MANAGER =====
data_manager = DataSourceManager()
source_type = data_manager.select_data_source()

data = None
if source_type == "local":
    data = data_manager.load_local_dataset()
elif source_type == "demo_hf":
    data = data_manager.load_demo_dataset()
elif source_type == "faers_live":
    data = data_manager.load_faers_dataset()

if data is not None:
    st.session_state.data_loaded = True
    st.session_state.raw_data = data
    st.session_state.metadata = data_manager.get_metadata()
    
    st.sidebar.success(f"✅ Data loaded: {len(data):,} records")
    with st.sidebar.expander("📋 Data Metadata"):
        st.json(data_manager.get_metadata())

# ===== MAIN TABS =====
if st.session_state.data_loaded:
    tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
        "🔍 Signal Detection",
        "🤖 ML Validation",
        "💡 Explainability",
        "📝 SAR Reports (RAG)",
        "📋 MLflow Tracking",
        "🗑️ Unlearning"
    ])
    
    # ===== TAB 1: SIGNAL DETECTION =====
    with tab1:
        st.header("Statistical Signal Detection")
        
        if st.button("🚀 Run Signal Detection", type="primary"):
            with st.spinner("Computing signals..."):
                stats_engine = DisproportionalityAnalysis()
                signals = stats_engine.compute_signals(st.session_state.raw_data)
                st.session_state.signals = signals
                
                st.success(f"✅ Detected {len(signals)} signals")
        
        if st.session_state.signals is not None:
            st.subheader("Top Signals")
            st.dataframe(st.session_state.signals.head(20), use_container_width=True)
    
    # ===== TAB 2: ML VALIDATION =====
    with tab2:
        st.header("ML-Based Signal Validation (SISA)")
        
        if st.session_state.signals is not None:
            num_shards = st.number_input("Number of Shards", 5, 20, 10)
            
            if st.button("🤖 Train ML Model", type="primary"):
                with st.spinner("Training..."):
                    trainer = SISATrainer(Path("models/sisa"))
                    results = trainer.train(st.session_state.signals, n_shards=num_shards)
                    st.session_state.trained_model = trainer
                    st.success("✅ Training complete!")
                    st.metric("AUC", f"{results['auc']:.4f}")
        else:
            st.warning("⚠️ Run signal detection first")
    
    # ===== TAB 3: EXPLAINABILITY =====
    with tab3:
        st.header("Model Explainability (SHAP)")
        if st.session_state.trained_model:
            st.info("SHAP analysis available after training")
        else:
            st.warning("Train model first")
    
    # ===== TAB 4: SAR REPORTS - WORKING VERSION =====
    with tab4:
        st.header("Signal Assessment Reports (RAG-based)")
        st.info("Generate comprehensive SAR reports with causality assessment and literature evidence")
        
        if st.session_state.signals is not None:
            top_signals = st.session_state.signals.head(50)
            signal_options = top_signals.apply(
                lambda x: f"{x['drug_name']} - {x['event_term']} (PRR: {x['prr']:.2f}, Cases: {x['case_count']})",
                axis=1
            ).tolist()
            
            selected_idx = st.selectbox("Select Signal for SAR", range(len(signal_options)), format_func=lambda i: signal_options[i])
            selected_signal = top_signals.iloc[selected_idx]
            
            col1, col2 = st.columns([2, 1])
            with col1:
                model_choice = st.selectbox(
                    "Select LLM Model",
                    ["llama3.2", "llama3.1", "mistral", "mixtral"],
                    help="Ollama model for SAR generation"
                )
            with col2:
                st.metric("Signal Strength", f"{selected_signal['prr']:.2f}")
            
            # GENERATE BUTTON
            if st.button("📝 Generate SAR Report", type="primary"):
                with st.spinner("Generating comprehensive SAR..."):
                    try:
                        from src.rag.sar_generator import SARGenerator
                        generator = SARGenerator(model_name=model_choice)
                        
                        service_ok, msg = generator.check_service()
                        if not service_ok:
                            st.error(f"❌ {msg}")
                            st.info("💡 Start Ollama: ollama serve")
                            st.stop()
                        
                        signal_data = {
                            'drug_name': selected_signal['drug_name'],
                            'event_name': selected_signal['event_term'],
                            'count': int(selected_signal['case_count']),
                            'prr': float(selected_signal['prr']),
                            'chi_square': float(selected_signal['chi2'])
                        }
                        
                        sar_result = generator.generate_sar(signal_data, model_choice)
                        
                        # Store in session state
                        st.session_state.current_sar = sar_result
                        st.session_state.current_signal_data = signal_data
                        st.session_state.current_model_choice = model_choice
                        
                        st.success("✅ SAR Generated Successfully!")
                        st.rerun()
                        
                    except Exception as e:
                        st.error(f"❌ SAR generation failed: {str(e)}")
                        with st.expander("Error Details"):
                            st.code(str(e))
            
            # Display SAR if exists
            if hasattr(st.session_state, 'current_sar') and st.session_state.current_sar is not None:
                st.markdown("---")
                
                sar_result = st.session_state.current_sar
                signal_data = st.session_state.current_signal_data
                model_choice = st.session_state.current_model_choice
                
                # Display Causality
                if sar_result.get('causality_assessment'):
                    with st.expander("🔬 Causality Assessment", expanded=True):
                        causality = sar_result['causality_assessment']
                        col1, col2, col3 = st.columns(3)
                        
                        if 'who_umc' in causality:
                            with col1:
                                st.metric("WHO-UMC", causality['who_umc']['category'])
                                st.caption(causality['who_umc']['rationale'])
                        
                        if 'naranjo' in causality:
                            with col2:
                                st.metric("Naranjo Score", f"{causality['naranjo']['score']}")
                                st.caption(causality['naranjo']['category'])
                        
                        if 'consensus' in causality:
                            with col3:
                                st.metric("Consensus", causality['consensus']['level'])
                
                # Display Literature
                if sar_result.get('literature_evidence') and sar_result['literature_evidence'].get('n_papers', 0) > 0:
                    with st.expander("📚 Literature Evidence", expanded=True):
                        lit = sar_result['literature_evidence']
                        st.markdown(lit['evidence_text'])
                
                # Display Report
                st.subheader("Generated SAR Report")
                st.markdown(sar_result['report'])
                
                with st.expander("📊 Evidence Sources"):
                    for i, source in enumerate(sar_result['evidence_sources'], 1):
                        st.write(f"{i}. {source}")
                
                # Download button
                from datetime import datetime
                timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
                clean_drug = signal_data['drug_name'].replace(' ', '_').replace('/', '_')[:50]
                clean_event = signal_data['event_name'].replace(' ', '_').replace('/', '_')[:50]
                report_filename = f"SAR_{clean_drug}_{clean_event}_{timestamp}.md"
                
                full_report = f"""# Signal Assessment Report (SAR)

## Drug-Event Pair
**Drug:** {signal_data['drug_name']}  
**Event:** {signal_data['event_name']}  
**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S IST')}

## Signal Metrics
- **PRR:** {signal_data['prr']:.2f}
- **Case Count:** {signal_data['count']}
- **Chi-Square:** {signal_data['chi_square']:.2f}
- **LLM Model:** {model_choice}

---

{sar_result['report']}

---

## Evidence Sources
"""
                for i, source in enumerate(sar_result['evidence_sources'], 1):
                    full_report += f"{i}. {source}\n"
                
                st.markdown("---")
                col1, col2 = st.columns([3, 1])
                
                with col1:
                    st.download_button(
                        label="📥 Download SAR Report",
                        data=full_report,
                        file_name=report_filename,
                        mime="text/markdown",
                        use_container_width=True
                    )
                
                with col2:
                    if st.button("🗑️ Clear", use_container_width=True):
                        st.session_state.current_sar = None
                        st.session_state.current_signal_data = None
                        st.session_state.current_model_choice = None
                        st.rerun()
        
        else:
            st.warning("⚠️ Run signal detection first (Tab 1)")
    
    # ===== TAB 5: MLFLOW TRACKING =====
    with tab5:
        st.header("MLflow Experiment Tracking")
        st.info("Track model training experiments, metrics, and artifacts")
        
        try:
            mlflow.set_tracking_uri(config.MLFLOW_TRACKING_URI)
            experiments = mlflow.search_experiments()
            
            if experiments:
                st.subheader("Experiments")
                for exp in experiments:
                    if exp.name != "Default":
                        with st.expander(f"📂 {exp.name}"):
                            st.write(f"Experiment ID: {exp.experiment_id}")
                            
                            runs = mlflow.search_runs(experiment_ids=[exp.experiment_id])
                            if not runs.empty:
                                st.write(f"Total Runs: {len(runs)}")
                                st.dataframe(runs[['run_id', 'start_time', 'status']].head(10))
                            else:
                                st.info("No runs found")
            else:
                st.info("No experiments found")
        
        except Exception as e:
            st.error(f"MLflow error: {str(e)}")
    
    # ===== TAB 6: UNLEARNING =====
    with tab6:
        st.header("Machine Unlearning (GDPR Article 17)")
        st.info("Remove specific data points and retrain affected model shards")
        
        if st.session_state.trained_model:
            case_id = st.text_input("Enter Case ID to Remove")
            
            if st.button("🗑️ Unlearn Case", type="primary"):
                if case_id:
                    with st.spinner("Unlearning..."):
                        try:
                            result = st.session_state.trained_model.unlearn_case(case_id)
                            st.success(f"✅ Case removed from shard {result['shard_id']}")
                            st.info(f"Retrained shard in {result['retrain_time']:.2f}s")
                        except Exception as e:
                            st.error(f"❌ Unlearning failed: {str(e)}")
                else:
                    st.warning("⚠️ Please enter a case ID")
        else:
            st.warning("⚠️ Train a model first (Tab 2)")

else:
    st.info("👈 Select a data source from the sidebar to begin")

# ===== SIDEBAR FOOTER =====
st.sidebar.markdown("---")
st.sidebar.markdown("### System Info")
st.sidebar.markdown(f"**Version:** {config.VALIDATION_STATUS['code_version']}")
st.sidebar.markdown(f"**MLflow URI:** {config.MLFLOW_TRACKING_URI}")
