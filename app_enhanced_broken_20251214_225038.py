"""
PV_Signal_ML Enhanced System - REGULATORY COMPLIANT VERSION
Complete MLflow audit trail for FDA/EMA/PMDA compliance
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
from src.utils.regulatory_tracker import tracker

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

# ===== SIDEBAR: SERVICE STATUS INDICATORS =====
st.sidebar.markdown("### 🔌 Service Status")

# Ollama Status
try:
    from src.rag.sar_generator import SARGenerator
    sar_gen = SARGenerator()
    ollama_ok, ollama_msg = sar_gen.check_service()
    if ollama_ok:
        st.sidebar.success(f"🟢 Ollama: {ollama_msg}")
    else:
        st.sidebar.error(f"🔴 Ollama: {ollama_msg}")
except Exception as e:
    st.sidebar.warning("⚪ Ollama: Not checked")

# MLflow Status
try:
    mlflow.set_tracking_uri(config.MLFLOW_TRACKING_URI)
    experiments = mlflow.search_experiments()
    st.sidebar.success(f"🟢 MLflow: Connected ({len(experiments)} experiments)")
except Exception as e:
    st.sidebar.error("🔴 MLflow: Connection failed")

st.sidebar.markdown("---")

# ===== DATA SOURCE MANAGER =====
data_manager = DataSourceManager()
source_type = data_manager.select_data_source()

# Load data based on source
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
    
    # Display metadata
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
    
    # ========== TAB 1: SIGNAL DETECTION WITH MLFLOW LOGGING ==========
    with tab1:
        st.header("Statistical Signal Detection")
        
        if st.button("🚀 Run Signal Detection", type="primary"):
            with st.spinner("Computing disproportionality signals..."):
                stats_engine = DisproportionalityAnalysis()
                signals = stats_engine.compute_signals(st.session_state.raw_data)
                st.session_state.signals = signals
                
                # Apply ML Causality Assessment
                st.info("🤖 Applying ML causality model to rank signals...")
                from src.ml.causality_scorer import CausalityScorer
                scorer = CausalityScorer()
                if scorer.load_model():
                    signals = scorer.score_signals(signals)
                    signals = signals.sort_values('causality_score', ascending=False)
                    st.success("✅ Causality assessment complete")
                    
                    col1, col2, col3 = st.columns(3)
                    high_causal = (signals['causality_score'] > 0.7).sum()
                    mod_causal = ((signals['causality_score'] >= 0.3) & (signals['causality_score'] <= 0.7)).sum()
                    low_causal = (signals['causality_score'] < 0.3).sum()
                    
                    col1.metric("High Causality", f"{high_causal:,}", delta="Priority")
                    col2.metric("Moderate Causality", f"{mod_causal:,}")
                    col3.metric("Low Causality", f"{low_causal:,}")
                else:
                    st.info("ℹ️ Using statistical signals only (ML model not trained)")
                
                st.success(f"✅ Detected {len(signals)} signals")
                
                # ===== LOG TO MLFLOW FOR REGULATORY COMPLIANCE =====
                try:
                    run_id = tracker.log_signal_detection(
                        n_signals=len(signals),
                        n_records=len(st.session_state.raw_data),
                        prr_threshold=2.0,
                        chi2_threshold=4.0,
                        top_signals=signals,
                        data_source=source_type,
                        date_range=st.session_state.metadata.get('date_range')
                    )
                    st.info(f"📊 Logged to MLflow: Run ID {run_id[:8]}")
                except Exception as e:
                    st.warning(f"⚠️ MLflow logging failed: {str(e)}")
        
        # Drug Portfolio Filter
        with st.sidebar.expander("🎯 Filter by Drug Portfolio"):
            st.write("**Targeted Surveillance**")
            st.write("Upload Excel with your portfolio drugs to focus analysis")
            
            uploaded_drug_list = st.file_uploader(
                "Upload drug list (Excel)",
                type=['xlsx', 'xls'],
                help="Excel file with drug names (any column)",
                key="drug_filter_upload"
            )
            
            if uploaded_drug_list:
                from src.utils.drug_filter import DrugFilter
                if 'drug_filter' not in st.session_state:
                    st.session_state.drug_filter = DrugFilter()
                
                drug_count = st.session_state.drug_filter.load_from_excel(uploaded_drug_list)
                if drug_count > 0:
                    st.success(f"✅ Loaded {drug_count} portfolio drugs")
                    with st.expander("View Portfolio"):
                        drugs = st.session_state.drug_filter.drug_list
                        st.write(drugs[:20])
                        if len(drugs) > 20:
                            st.write(f"... and {len(drugs)-20} more")
        
        # Apply drug filter if active
        if st.session_state.signals is not None and 'drug_filter' in st.session_state:
            if st.session_state.drug_filter.active:
                original_count = len(st.session_state.signals)
                filtered_signals = st.session_state.drug_filter.filter_signals(st.session_state.signals)
                st.info(f"🎯 **Portfolio Filter Active:** {len(filtered_signals):,} signals from {original_count:,} total")
                
                coverage = st.session_state.drug_filter.get_coverage_report(st.session_state.signals)
                col1, col2, col3 = st.columns(3)
                col1.metric("Portfolio Drugs", coverage['total_portfolio'])
                col2.metric("Drugs with Signals", coverage['drugs_with_signals'])
                col3.metric("Drugs without Signals", coverage['drugs_without_signals'])
                
                if coverage['drugs_without_signals'] > 0:
                    with st.expander(f"View {coverage['drugs_without_signals']} drugs without signals"):
                        st.write(coverage['missing_list'])
                
                st.session_state.signals = filtered_signals
        
        # Display signals
        if st.session_state.signals is not None:
            st.subheader("Top Signals")
            
            col1, col2 = st.columns(2)
            with col1:
                prr_threshold = st.slider("PRR Threshold", 1.0, 10.0, 2.0, 0.1)
            with col2:
                min_cases = st.slider("Min Case Count", 3, 50, 3)
            
            filtered = st.session_state.signals[
                (st.session_state.signals['prr'] >= prr_threshold) &
                (st.session_state.signals['case_count'] >= min_cases)
            ]
            
            st.dataframe(
                filtered[['drug_name', 'event_term', 'case_count', 'prr', 'chi2'] +
                        (['causality_score', 'causality_level'] if 'causality_score' in filtered.columns else []) +
                        ['is_signal_prr']],
                use_container_width=True
            )
            
            csv = filtered.to_csv(index=False)
            st.download_button("📥 Download Signals CSV", csv, "signals.csv", "text/csv")
    
    # ========== TAB 2: ML VALIDATION (Already has MLflow) ==========
    with tab2:
        st.header("ML-Based Signal Validation (SISA)")
        st.info("""
        **SISA (Sharded, Isolated, Sliced, Aggregated) Training**
        - Enables efficient machine unlearning for GDPR compliance
        - Trains ensemble of models on data shards
        - Allows removal of individual cases without full retraining
        """)
        
        if st.session_state.signals is not None:
            col1, col2 = st.columns(2)
            with col1:
                enable_sisa = st.checkbox("Enable SISA Training", value=config.ENABLE_SISA)
            with col2:
                num_shards = st.number_input("Number of Shards", 5, 20, config.NUM_SHARDS)
            
            if st.button("🤖 Train ML Model", type="primary"):
                with st.spinner("Training SISA ensemble..."):
                    signals_df = st.session_state.signals
                    
                    if enable_sisa:
                        model_dir = Path("models/sisa")
                        trainer = SISATrainer(model_dir)
                        
                        results = trainer.train(signals_df, n_shards=num_shards)
                        
                        st.session_state.trained_model = trainer
                        st.session_state.training_results = results
                        
                        st.success("✅ SISA Model Trained Successfully!")
                        
                        col1, col2, col3 = st.columns(3)
                        with col1:
                            st.metric("AUC Score", f"{results['auc']:.4f}")
                        with col2:
                            st.metric("Training Samples", results['n_train'])
                        with col3:
                            st.metric("Features", results['n_features'])
                        
                        with st.expander("📊 Training Details"):
                            st.json({
                                'auc': results['auc'],
                                'n_features': results['n_features'],
                                'feature_names': results['feature_names'],
                                'n_train': results['n_train'],
                                'n_test': results['n_test'],
                                'n_shards': results['n_shards']
                            })
                            st.code(results['report'])
                    else:
                        st.warning("SISA disabled - using standard training")
        else:
            st.warning("⚠️ Run signal detection first (Tab 1)")
    
    # ========== TAB 3: EXPLAINABILITY ==========
    with tab3:
        st.header("Model Explainability (SHAP)")
        
        if st.session_state.trained_model is not None:
            if st.button("💡 Generate SHAP Explanations"):
                if not hasattr(st.session_state.trained_model, 'feature_pipeline'):
                    st.error("❌ Feature pipeline not found. Please retrain the model.")
                    st.stop()
                
                with st.spinner("Computing SHAP values..."):
                    try:
                        signals_df = st.session_state.signals
                        feature_pipeline = st.session_state.trained_model.feature_pipeline
                        
                        X = feature_pipeline.transform(signals_df)
                        st.info(f"📊 Using {len(feature_pipeline.feature_names)} features: {feature_pipeline.feature_names}")
                        
                        analyzer = SHAPAnalyzer()
                        
                        if hasattr(st.session_state.trained_model, 'model'):
                            model = st.session_state.trained_model.model
                        elif hasattr(st.session_state.trained_model, 'shard_models'):
                            model = st.session_state.trained_model.shard_models[0]
                        else:
                            st.error("Model structure incompatible")
                            st.stop()
                        
                        analyzer.compute_shap_values(model, X)
                        summary = analyzer.generate_global_summary()
                        st.success("✅ SHAP analysis complete!")
                        
                        st.subheader("Global Feature Importance")
                        st.image(summary['summary_plot'])
                        
                        st.subheader("Top Features")
                        top_features = pd.DataFrame(summary['top_features'])
                        st.dataframe(top_features, use_container_width=True)
                        
                        st.info(f"Analysis based on {summary['n_samples_used']} samples")
                        
                    except Exception as e:
                        st.error(f"❌ SHAP computation failed: {str(e)}")
                        with st.expander("Error Details"):
                            st.code(str(e))
        else:
            st.warning("⚠️ Train ML model first (Tab 2)")
    
    # ========== TAB 4: SAR REPORTS WITH MLFLOW LOGGING ==========
    # ===== TAB 4: SAR REPORTS - FIXED =====
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
                        
                        # ✅ STORE IN SESSION STATE
                        st.session_state.current_sar = sar_result
                        st.session_state.current_signal_data = signal_data
                        st.session_state.current_model_choice = model_choice
                        
                        st.success("✅ SAR Generated Successfully!")
                        st.rerun()
                        
                    except Exception as e:
                        st.error(f"❌ SAR generation failed: {str(e)}")
                        with st.expander("Error Details"):
                            st.code(str(e))
            
            # ===== DISPLAY SAR (if exists in session_state) =====
            if st.session_state.current_sar is not None:
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
                literature_count = 0
                if sar_result.get('literature_evidence'):
                    literature_count = sar_result['literature_evidence'].get('n_papers', 0)
                    if literature_count > 0:
                        with st.expander("📚 Literature Evidence", expanded=True):
                            lit = sar_result['literature_evidence']
                            st.markdown(lit['evidence_text'])
                
                # Display Report
                st.subheader("Generated SAR Report")
                st.markdown(sar_result['report'])
                
                with st.expander("📊 Evidence Sources"):
                    for i, source in enumerate(sar_result['evidence_sources'], 1):
                        st.write(f"{i}. {source}")
                
                # ===== SAVE & DOWNLOAD =====
                st.markdown("---")
                from datetime import datetime
                
                timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
                clean_drug = signal_data['drug_name'].replace(' ', '_').replace('/', '_')[:50]
                clean_event = signal_data['event_name'].replace(' ', '_').replace('/', '_')[:50]
                report_filename = f"SAR_{clean_drug}_{clean_event}_{timestamp}.txt"
                
                # Format report
                full_report = f'''# Signal Assessment Report (SAR)

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
'''
                for i, source in enumerate(sar_result['evidence_sources'], 1):
                    full_report += f"{i}. {source}\n"
                
                # Save locally
                from pathlib import Path
                reports_dir = Path("outputs/sar_reports")
                reports_dir.mkdir(parents=True, exist_ok=True)
                local_path = reports_dir / report_filename
                
                try:
                    with open(local_path, 'w', encoding='utf-8') as f:
                        f.write(full_report)
                    save_ok = True
                except Exception as e:
                    save_ok = False
                    st.error(f"Local save failed: {str(e)}")
                
                # Buttons
                col1, col2, col3 = st.columns([2, 2, 1])
                
                with col1:
                    st.download_button(
                        label="📥 Download SAR Report",
                        data=full_report,
                        file_name=report_filename,
                        mime="text/plain",
                        use_container_width=True
                    )
                
                with col2:
                    if save_ok:
                        st.success(f"✅ Saved: {report_filename}")
                    else:
                        st.error("❌ Save failed")
                
                with col3:
                    if st.button("🗑️ Clear", use_container_width=True):
                        st.session_state.current_sar = None
                        st.session_state.current_signal_data = None
                        st.session_state.current_model_choice = None
                        st.rerun()
                
                # MLflow logging
                try:
                    from src.utils.regulatory_tracker import tracker
                    run_id = tracker.log_sar_generation(
                        drug_name=signal_data['drug_name'],
                        event_name=signal_data['event_name'],
                        prr=signal_data['prr'],
                        case_count=signal_data['count'],
                        chi_square=signal_data['chi_square'],
                        causality_assessment=sar_result.get('causality_assessment'),
                        literature_count=literature_count,
                        model_used=model_choice,
                        report_path=str(local_path) if save_ok else None
                    )
                    st.info(f"📊 MLflow: {run_id[:8]}")
                except Exception as e:
                    st.warning(f"⚠️ MLflow failed: {str(e)}")
        
        else:
            st.warning("⚠️ Run signal detection first (Tab 1)")
        # ===== TAB 5: MLFLOW TRACKING =====
    with tab5:
        st.header("MLflow Experiment Tracking")
        st.info("Track model training experiments, metrics, and artifacts")
        
        try:
            mlflow.set_tracking_uri(config.MLFLOW_TRACKING_URI)
            
            # Display MLflow UI link at the top
            st.markdown("---")
            col1, col2, col3 = st.columns([2, 2, 1])
            
            with col1:
                st.subheader("🔗 MLflow UI Access")
                
                # Determine MLflow UI URL
                if config.MLFLOW_TRACKING_URI.startswith("file:"):
                    mlflow_ui_url = "http://localhost:5000"
                    st.markdown(f"""
                    **Local MLflow Server**
                    
                    To view experiments in MLflow UI:
                    1. Open PowerShell in project directory
                    2. Run: `mlflow ui`
                    3. Open: [{mlflow_ui_url}]({mlflow_ui_url})
                    """)
                    
                    if st.button("📋 Copy MLflow Command", use_container_width=True):
                        st.code("mlflow ui", language="bash")
                        st.info("💡 Run this in PowerShell to start MLflow UI")
                else:
                    mlflow_ui_url = config.MLFLOW_TRACKING_URI.replace("http://", "").replace("https://", "")
                    st.markdown(f"**MLflow Server:** [{mlflow_ui_url}]({config.MLFLOW_TRACKING_URI})")
            
            with col2:
                st.subheader("📊 Quick Stats")
                experiments = mlflow.search_experiments()
                total_experiments = len([exp for exp in experiments if exp.name != "Default"])
                
                st.metric("Total Experiments", total_experiments)
                
                # Count total runs across all experiments
                total_runs = 0
                for exp in experiments:
                    runs = mlflow.search_runs(experiment_ids=[exp.experiment_id])
                    total_runs += len(runs)
                
                st.metric("Total Runs", total_runs)
            
            with col3:
                st.subheader("🔄 Actions")
                if st.button("🔄 Refresh", use_container_width=True):
                    st.rerun()
            
            st.markdown("---")
            
            # Display experiments
            if experiments:
                st.subheader("📂 Experiments")
                
                # Create experiment table with links
                exp_data = []
                for exp in experiments:
                    if exp.name != "Default":
                        runs = mlflow.search_runs(experiment_ids=[exp.experiment_id])
                        exp_data.append({
                            "Experiment": exp.name,
                            "ID": exp.experiment_id,
                            "Runs": len(runs),
                            "Last Updated": runs['start_time'].max().strftime('%Y-%m-%d %H:%M') if not runs.empty and 'start_time' in runs.columns else "N/A"
                        })
                
                exp_df = pd.DataFrame(exp_data)
                st.dataframe(exp_df, use_container_width=True)
                
                # Experiment selector
                exp_names = [exp.name for exp in experiments if exp.name != "Default"]
                
                if exp_names:
                    selected_exp = st.selectbox("🔍 Select Experiment to View Details", exp_names)
                    
                    # Get experiment ID
                    exp_id = [exp.experiment_id for exp in experiments if exp.name == selected_exp][0]
                    runs = mlflow.search_runs(experiment_ids=[exp_id])
                    
                    if not runs.empty:
                        st.subheader(f"📋 Runs for: {selected_exp}")
                        st.caption(f"Total Runs: {len(runs)}")
                        
                        # Display runs table
                        display_cols = ['run_id', 'start_time', 'status']
                        metric_cols = [col for col in runs.columns if col.startswith('metrics.')]
                        param_cols = [col for col in runs.columns if col.startswith('params.')]
                        
                        display_cols.extend(metric_cols[:5])  # Show first 5 metrics
                        display_cols.extend(param_cols[:5])   # Show first 5 params
                        
                        # Make run_id clickable
                        runs_display = runs[display_cols].copy()
                        runs_display['run_id'] = runs_display['run_id'].apply(lambda x: x[:8])
                        
                        st.dataframe(runs_display, use_container_width=True)
                        
                        # Run selector
                        if len(runs) > 0:
                            run_idx = st.selectbox(
                                "🔍 Select Run for Detailed View", 
                                range(len(runs)), 
                                format_func=lambda i: f"Run {i+1}: {runs.iloc[i]['run_id'][:8]} - {runs.iloc[i]['start_time']}"
                            )
                            
                            selected_run = runs.iloc[run_idx]
                            
                            with st.expander("📊 Run Details", expanded=True):
                                col1, col2 = st.columns(2)
                                
                                with col1:
                                    st.markdown("**⚙️ Parameters**")
                                    params = {
                                        k.replace('params.', ''): v 
                                        for k, v in selected_run.items() 
                                        if k.startswith('params.') and pd.notna(v)
                                    }
                                    if params:
                                        st.json(params)
                                    else:
                                        st.info("No parameters logged")
                                
                                with col2:
                                    st.markdown("**📈 Metrics**")
                                    metrics = {
                                        k.replace('metrics.', ''): v 
                                        for k, v in selected_run.items() 
                                        if k.startswith('metrics.') and pd.notna(v)
                                    }
                                    if metrics:
                                        st.json(metrics)
                                    else:
                                        st.info("No metrics logged")
                                
                                # Show run metadata
                                st.markdown("**ℹ️ Run Metadata**")
                                metadata = {
                                    "Run ID": selected_run['run_id'],
                                    "Start Time": selected_run['start_time'],
                                    "End Time": selected_run.get('end_time', 'N/A'),
                                    "Status": selected_run['status'],
                                    "Experiment ID": selected_run['experiment_id']
                                }
                                st.json(metadata)
                                
                                # MLflow UI direct link
                                run_id = selected_run['run_id']
                                mlflow_run_url = f"http://localhost:5000/#/experiments/{exp_id}/runs/{run_id}"
                                st.markdown(f"🔗 [Open in MLflow UI]({mlflow_run_url})")
                    else:
                        st.info("No runs found for this experiment")
            else:
                st.info("No experiments found. Train a model or run signal detection to create experiments.")
        
        except Exception as e:
            st.error(f"MLflow connection error: {str(e)}")
            st.info("Make sure MLflow tracking server is running or using local file store")
            
            with st.expander("🔧 Troubleshooting"):
                st.markdown("""
                **To start MLflow UI:**
                1. Open PowerShell in your project directory
                2. Run: `mlflow ui`
                3. Open http://localhost:5000 in your browser
                
                **To check tracking URI:**
                Current URI: `{config.MLFLOW_TRACKING_URI}`
                """)

    # ===== TAB 6: UNLEARNING =====
    else:
                    st.info("No runs found for this experiment")
            else:
                st.info("No experiments found. Train a model to create experiments.")
        
        except Exception as e:
            st.error(f"MLflow connection error: {str(e)}")
            st.info("Make sure MLflow tracking server is running or using local file store")
    
    # ========== TAB 6: UNLEARNING WITH PROGRESS + MLFLOW LOGGING ==========
    with tab6:
        st.header("🗑️ Machine Unlearning (GDPR/HIPAA Compliance)")
        st.markdown("*Remove specific patient cases and retrain affected model shards*")
        
        if 'trained_model' not in st.session_state or st.session_state.trained_model is None:
            st.warning("⚠️ No trained SISA model found. Please train a model in Tab 2 first.")
        else:
            # Safe model info display
            if hasattr(st.session_state.trained_model, 'shard_models'):
                n_shards = len(st.session_state.trained_model.shard_models)
            elif hasattr(st.session_state.trained_model, 'n_shards'):
                n_shards = st.session_state.trained_model.n_shards
            else:
                n_shards = 'Unknown'
            
            st.info(f"📊 Current model: {n_shards} shards trained")
            
            case_id_input = st.text_input(
                "Enter Case ID to Remove:",
                placeholder="e.g., 1501",
                help="Enter the primary ID of the case to remove from the system"
            )
            
            if st.button("🗑️ Unlearn Case", type="primary", use_container_width=True):
                if not case_id_input:
                    st.error("❌ Please enter a valid Case ID")
                else:
                    try:
                        case_id = int(case_id_input)
                        
                        # Progress tracking UI
                        progress_bar = st.progress(0.0)
                        status_text = st.empty()
                        
                        def update_progress(progress, message):
                            progress_bar.progress(progress)
                            status_text.info(f"🔄 {message}")
                        
                        # Call unlearn with callback
                        with st.spinner("Processing unlearning request..."):
                            result = st.session_state.trained_model.unlearn(
                                case_id, 
                                progress_callback=update_progress
                            )
                        
                        progress_bar.progress(1.0)
                        
                        # Display results
                        if result.get('success'):
                            status_text.success("✅ Unlearning complete!")
                            st.success(f"✅ **Unlearning Successful**")
                            st.markdown(f"""
                            - **Case ID:** {case_id}
                            - **Shard:** {result.get('shard_id', 'N/A')}
                            - **Cases Removed:** {result.get('cases_removed', 1)}
                            - **Retraining Time:** {result.get('retrain_time', 0.0):.2f}s
                            """)
                            
                            # ===== LOG TO MLFLOW FOR REGULATORY COMPLIANCE =====
                            try:
                                run_id = tracker.log_unlearning_operation(
                                    case_id=case_id,
                                    shard_id=result.get('shard_id'),
                                    cases_removed=result.get('cases_removed', 1),
                                    retrain_time=result.get('retrain_time', 0.0),
                                    success=True,
                                    error_message=None,
                                    model_version_before=f"v1_shard_{result.get('shard_id')}",
                                    model_version_after=f"v1_shard_{result.get('shard_id')}_unlearned"
                                )
                                st.info(f"📊 Logged to MLflow: Run ID {run_id[:8]}")
                            except Exception as e:
                                st.warning(f"⚠️ MLflow logging failed: {str(e)}")
                        else:
                            status_text.error(f"❌ {result.get('message', 'Unknown error')}")
                            st.error(f"❌ Unlearning failed: {result.get('message', 'Unknown error')}")
                            
                            # Log failure to MLflow
                            try:
                                tracker.log_unlearning_operation(
                                    case_id=case_id,
                                    success=False,
                                    error_message=result.get('message', 'Unknown error')
                                )
                            except:
                                pass
                            
                    except ValueError:
                        st.error("❌ Case ID must be a valid number")
                    except Exception as e:
                        st.error(f"❌ Error during unlearning: {str(e)}")
                        progress_bar.empty()
                        status_text.empty()
                        
                        # Log error to MLflow
                        try:
                            tracker.log_unlearning_operation(
                                case_id=int(case_id_input) if case_id_input.isdigit() else 0,
                                success=False,
                                error_message=str(e)
                            )
                        except:
                            pass

else:
    st.info("👈 Select a data source from the sidebar to begin")

# ===== SIDEBAR FOOTER =====
st.sidebar.markdown("---")
st.sidebar.markdown("### 📊 System Info")
st.sidebar.markdown(f"**Version:** {config.VALIDATION_STATUS['code_version']}")
st.sidebar.markdown(f"**MLflow URI:** {config.MLFLOW_TRACKING_URI}")

try:
    from src.rag.sar_generator import SARGenerator
    sar_gen = SARGenerator()
    model_ok, model_msg = sar_gen.check_model()
    st.sidebar.markdown(f"**Ollama Model:** {sar_gen.model_name}")
except:
    pass





