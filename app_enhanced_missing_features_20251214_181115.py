"""
PV_Signal_ML Enhanced System
Main Streamlit Application with SISA, SHAP, SAR (RAG), MLflow, and Unlearning
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
    page_title="PV Signal ML - Enhanced",
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

# ===== MAIN TABS (CORRECTED - 6 TABS) =====
if st.session_state.data_loaded:
    
    # CORRECTED TAB STRUCTURE - 6 TABS WITHOUT REPORTS
    tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
        "🔍 Signal Detection",
        "🤖 ML Validation",
        "💡 Explainability",
        "📝 SAR Reports (RAG)",
        "📋 MLflow Tracking",
        "🗑️ Unlearning"
    ])
    
    # ========== TAB 1: SIGNAL DETECTION ==========
    with tab1:
        st.header("Statistical Signal Detection")
        
        if st.button("🚀 Run Signal Detection", type="primary"):
            with st.spinner("Computing disproportionality signals..."):
                stats_engine = DisproportionalityAnalysis()
                signals = stats_engine.compute_signals(st.session_state.raw_data)
                st.session_state.signals = signals
                
                # Stage 2: Apply ML Causality Assessment
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
        
        # Display signals if available
        if st.session_state.signals is not None:
            st.subheader("Top Signals")
            
            # Filter controls
            col1, col2 = st.columns(2)
            with col1:
                prr_threshold = st.slider("PRR Threshold", 1.0, 10.0, 2.0, 0.1)
            with col2:
                min_cases = st.slider("Min Case Count", 3, 50, 3)
            
            # Filter signals
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
            
            # Download button
            csv = filtered.to_csv(index=False)
            st.download_button(
                "📥 Download Signals CSV",
                csv,
                "signals.csv",
                "text/csv"
            )
    
    # ========== TAB 2: ML VALIDATION ==========
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
                        
                        # Train SISA model
                        results = trainer.train(signals_df, n_shards=num_shards)
                        st.session_state.trained_model = trainer
                        st.session_state.training_results = results
                        
                        st.success("✅ SISA Model Trained Successfully!")
                        
                        # Display metrics
                        col1, col2, col3 = st.columns(3)
                        with col1:
                            st.metric("AUC Score", f"{results['auc']:.4f}")
                        with col2:
                            st.metric("Training Samples", results['n_train'])
                        with col3:
                            st.metric("Features", results['n_features'])
                        
                        # Show detailed results
                        with st.expander("📊 Training Details"):
                            st.json({
                                'auc': results['auc'],
                                'n_features': results['n_features'],
                                'feature_names': results['feature_names'],
                                'n_train': results['n_train'],
                                'n_test': results['n_test']
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
                # Check feature pipeline exists
                if not hasattr(st.session_state.trained_model, 'feature_pipeline'):
                    st.error("❌ Feature pipeline not found. Please retrain the model.")
                    st.stop()
                
                with st.spinner("Computing SHAP values..."):
                    try:
                        signals_df = st.session_state.signals
                        feature_pipeline = st.session_state.trained_model.feature_pipeline
                        
                        # CRITICAL: Use feature pipeline
                        X = feature_pipeline.transform(signals_df)
                        
                        st.info(f"📊 Using {len(feature_pipeline.feature_names)} features: {feature_pipeline.feature_names}")
                        
                        analyzer = SHAPAnalyzer()
                        
                        # Get model
                        if hasattr(st.session_state.trained_model, 'model'):
                            model = st.session_state.trained_model.model
                        elif hasattr(st.session_state.trained_model, 'shard_models'):
                            model = st.session_state.trained_model.shard_models[0]
                        else:
                            st.error("Model structure incompatible")
                            st.stop()
                        
                        # SHAP computation
                        analyzer.compute_shap_values(model, X)
                        
                        # Generate global summary
                        summary = analyzer.generate_global_summary()
                        
                        st.success("✅ SHAP analysis complete!")
                        
                        # Display global summary
                        st.subheader("Global Feature Importance")
                        st.image(summary['summary_plot'])
                        
                        # Display top features
                        st.subheader("Top Features")
                        top_features = pd.DataFrame(summary['top_features'])
                        st.dataframe(top_features, use_container_width=True)
                        
                        st.info(f"Analysis based on {summary['n_samples_used']} samples")
                        
                    except Exception as e:
                        st.error(f"❌ SHAP computation failed: {str(e)}")
                        st.code(str(e))
        else:
            st.warning("⚠️ Train ML model first (Tab 2)")
    
    # ========== TAB 4: SAR REPORTS (RAG) ==========
    with tab4:
        st.header("Signal Assessment Reports (RAG-based)")
        st.info("Generate comprehensive SAR reports with causality assessment and literature evidence")
        
        if st.session_state.signals is not None:
            # Signal selection
            top_signals = st.session_state.signals.head(50)
            signal_options = top_signals.apply(
                lambda x: f"{x['drug_name']} - {x['event_term']} (PRR: {x['prr']:.2f}, Cases: {x['case_count']})",
                axis=1
            ).tolist()
            
            selected_idx = st.selectbox("Select Signal for SAR", range(len(signal_options)), format_func=lambda i: signal_options[i])
            selected_signal = top_signals.iloc[selected_idx]
            
            # Model selection
            col1, col2 = st.columns([2, 1])
            with col1:
                model_choice = st.selectbox(
                    "Select LLM Model",
                    ["llama3.2", "llama3.1", "mistral", "mixtral"],
                    help="Ollama model for SAR generation"
                )
            with col2:
                st.metric("Signal Strength", f"{selected_signal['prr']:.2f}")
            
            if st.button("📝 Generate SAR Report", type="primary"):
                with st.spinner("Generating comprehensive SAR..."):
                    try:
                        from src.rag.sar_generator import SARGenerator
                        
                        generator = SARGenerator(model_name=model_choice)
                        
                        # Check service
                        service_ok, msg = generator.check_service()
                        if not service_ok:
                            st.error(f"❌ {msg}")
                            st.info("Start Ollama: ollama serve")
                            st.stop()
                        
                        # Prepare signal data
                        signal_data = {
                            'drug_name': selected_signal['drug_name'],
                            'event_name': selected_signal['event_term'],
                            'count': int(selected_signal['case_count']),
                            'prr': float(selected_signal['prr']),
                            'chi_square': float(selected_signal['chi2'])
                        }
                        
                        # Generate SAR
                        sar_result = generator.generate_sar(signal_data, model_choice)
                        
                        st.success("✅ SAR Generated Successfully!")
                        
                        # Display causality if available
                        if sar_result.get('causality_assessment'):
                            with st.expander("🔬 Causality Assessment", expanded=True):
                                causality = sar_result['causality_assessment']
                                
                                col1, col2, col3 = st.columns(3)
                                if 'who_umc' in causality:
                                    with col1:
                                        st.metric("WHO-UMC", causality['who_umc']['category'])
                                if 'naranjo' in causality:
                                    with col2:
                                        st.metric("Naranjo Score", f"{causality['naranjo']['score']}")
                                if 'consensus' in causality:
                                    with col3:
                                        st.metric("Consensus", causality['consensus']['level'])
                        
                        # Display literature if available
                        if sar_result.get('literature_evidence') and sar_result['literature_evidence'].get('n_papers', 0) > 0:
                            with st.expander("📚 Literature Evidence", expanded=True):
                                lit = sar_result['literature_evidence']
                                st.markdown(lit['evidence_text'])
                        
                        # Display generated report
                        st.subheader("Generated SAR Report")
                        st.markdown(sar_result['report'])
                        
                        # Evidence sources
                        with st.expander("📊 Evidence Sources"):
                            for i, source in enumerate(sar_result['evidence_sources'], 1):
                                st.write(f"{i}. {source}")
                        
                        # Save option
                        if st.button("💾 Save Report"):
                            filepath = generator.save_report(sar_result)
                            st.success(f"✅ Report saved: {filepath}")
                        
                    except Exception as e:
                        st.error(f"❌ SAR generation failed: {str(e)}")
                        st.code(str(e))
        else:
            st.warning("⚠️ Run signal detection first (Tab 1)")
    
    # ========== TAB 5: MLFLOW TRACKING ==========
    with tab5:
        st.header("MLflow Experiment Tracking")
        st.info("Track model training experiments, metrics, and artifacts")
        
        try:
            # Set MLflow tracking URI
            mlflow.set_tracking_uri(config.MLFLOW_TRACKING_URI)
            
            # Get all experiments
            experiments = mlflow.search_experiments()
            
            if experiments:
                st.subheader("Experiments")
                exp_df = pd.DataFrame([{
                    'Name': exp.name,
                    'ID': exp.experiment_id,
                    'Artifact Location': exp.artifact_location
                } for exp in experiments])
                st.dataframe(exp_df, use_container_width=True)
                
                # Select experiment
                exp_names = [exp.name for exp in experiments]
                selected_exp = st.selectbox("Select Experiment", exp_names)
                
                # Get runs for selected experiment
                exp_id = [exp.experiment_id for exp in experiments if exp.name == selected_exp][0]
                runs = mlflow.search_runs(experiment_ids=[exp_id])
                
                if not runs.empty:
                    st.subheader(f"Runs for {selected_exp}")
                    
                    # Display runs
                    display_cols = ['run_id', 'start_time', 'status'] + [col for col in runs.columns if col.startswith('metrics.') or col.startswith('params.')]
                    st.dataframe(runs[display_cols], use_container_width=True)
                    
                    # Run details
                    if len(runs) > 0:
                        run_idx = st.selectbox("Select Run", range(len(runs)), format_func=lambda i: f"Run {i+1}: {runs.iloc[i]['run_id'][:8]}")
                        selected_run = runs.iloc[run_idx]
                        
                        with st.expander("Run Details", expanded=True):
                            col1, col2 = st.columns(2)
                            with col1:
                                st.write("**Parameters:**")
                                params = {k.replace('params.', ''): v for k, v in selected_run.items() if k.startswith('params.') and pd.notna(v)}
                                st.json(params)
                            with col2:
                                st.write("**Metrics:**")
                                metrics = {k.replace('metrics.', ''): v for k, v in selected_run.items() if k.startswith('metrics.') and pd.notna(v)}
                                st.json(metrics)
                else:
                    st.info("No runs found for this experiment")
            else:
                st.info("No experiments found. Train a model to create experiments.")
        
        except Exception as e:
            st.error(f"MLflow connection error: {str(e)}")
            st.info("Make sure MLflow tracking server is running or using local file store")
    
    # ========== TAB 6: UNLEARNING ==========
    with tab6:
        st.header("Machine Unlearning (GDPR Right to be Forgotten)")
        
        st.markdown("""
        **SISA Unlearning Process:**
        1. Enter Case ID to remove
        2. System identifies affected data shard
        3. Retrain only affected shard from checkpoint
        4. Update ensemble without full retraining
        """)
        
        if st.session_state.trained_model is not None:
            st.subheader("Current Model Version")
            col1, col2 = st.columns(2)
            with col1:
                st.metric("Model Version", "1")
            with col2:
                if hasattr(st.session_state.trained_model, 'shard_models'):
                    st.metric("Active Shards", len(st.session_state.trained_model.shard_models))
            
            st.subheader("Unlearn Case")
            case_id = st.text_input("Enter Case ID to Unlearn")
            
            if st.button("🗑️ Unlearn Case", type="primary"):
                if case_id:
                    with st.spinner(f"Unlearning case {case_id}..."):
                        try:
                            result = st.session_state.trained_model.unlearn(case_id)
                            
                            if result['status'] == 'unlearned':
                                st.success(f"✅ Case {case_id} successfully unlearned!")
                                st.json(result)
                            else:
                                st.warning(f"Case {case_id} not found in training data")
                        except Exception as e:
                            st.error(f"Unlearning failed: {str(e)}")
                else:
                    st.error("Please enter a case ID")
        else:
            st.warning("⚠️ Train SISA model first (Tab 2)")

else:
    st.info("👈 Select a data source from the sidebar to begin")

# ===== FOOTER =====
st.sidebar.markdown("---")
st.sidebar.markdown(f"**Version:** {config.VALIDATION_STATUS['code_version']}")
st.sidebar.markdown(f"**MLflow:** {config.MLFLOW_TRACKING_URI}")
