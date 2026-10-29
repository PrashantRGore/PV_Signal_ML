"""
PV Signal ML - Enhanced System
Version: 2.0.0-canonical (FIXED: No re-downloading)
Canonical Signal Schema with Automatic Cache Detection
"""

import streamlit as st
import pandas as pd
import mlflow
from pathlib import Path

import config
from src.data.data_source_manager import DataSourceManager
from src.data.canonical_schema import canonicalize_signals, CanonicalSchema
from src.utils.data_prep import prepare_data_for_signal_detection
from src.stats_engine.disproportionality import DisproportionalityAnalysis
from src.ml.sisa_trainer import SISATrainer
from src.explainability.shap_analysis import SHAPAnalyzer
from src.utils.logger import setup_logger

logger = setup_logger(__name__)

# ===== PAGE CONFIG =====
st.set_page_config(
    page_title="PV Signal ML - Enhanced",
    page_icon="💊",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ===== DISCLAIMER =====
st.warning(config.DISCLAIMER_TEXT)

# ===== TITLE =====
st.title("💊 Pharmacovigilance Signal Detection with SISA & SHAP")
st.markdown("Enhanced System with Machine Unlearning, Explainable AI, RAG-based SAR, and MLflow")

# ===== SESSION STATE INITIALIZATION =====
if "data_loaded" not in st.session_state:
    st.session_state.data_loaded = False
if "signals" not in st.session_state:
    st.session_state.signals = None
if "raw_data" not in st.session_state:
    st.session_state.raw_data = None
if "trained_model" not in st.session_state:
    st.session_state.trained_model = None
if "training_results" not in st.session_state:
    st.session_state.training_results = None
if "analysis_mode" not in st.session_state:
    st.session_state.analysis_mode = "tier1"
if "current_sar" not in st.session_state:
    st.session_state.current_sar = None
if "current_signal_data" not in st.session_state:
    st.session_state.current_signal_data = None
if "current_model_choice" not in st.session_state:
    st.session_state.current_model_choice = None
if "metadata" not in st.session_state:
    st.session_state.metadata = {}

# ===== SIDEBAR: SERVICE STATUS =====
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
except Exception:
    st.sidebar.warning("⚠️ Ollama: Not checked")

# MLflow Status
try:
    mlflow.set_tracking_uri(config.MLFLOW_TRACKING_URI)
    experiments = mlflow.search_experiments()
    st.sidebar.success(f"🟢 MLflow: Connected ({len(experiments)} experiments)")
except Exception:
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

# FIXED: Auto-detect signals and cache them
if data is not None:
    st.session_state.data_loaded = True
    st.session_state.raw_data = data
    st.session_state.metadata = data_manager.get_metadata()
    
    # Check if data already contains signal columns (Demo/FAERS sources)
    has_signal_columns = all(col in data.columns for col in ["prr", "ror", "chi2"])
    
    if has_signal_columns:
        # Data already contains computed signals - cache them immediately
        logger.info(f"✅ Detected pre-computed signals in {source_type} data")
        st.session_state.signals = canonicalize_signals(
            data, 
            source=source_type,
            preserve_tier2=True
        )
        st.sidebar.success(f"✅ Data loaded: {len(data):,} signals (cached)")
    else:
        # Raw ICSR data - needs signal computation
        logger.info(f"📊 Loaded raw ICSR data from {source_type}")
        st.sidebar.success(f"✅ Data loaded: {len(data):,} records")
    
    # Show metadata
    with st.sidebar.expander("📋 Data Metadata"):
        st.json(data_manager.get_metadata())

# ===== MAIN TABS =====
if st.session_state.data_loaded:
    tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
        "📊 Signal Detection",
        "🤖 ML Validation",
        "💡 Explainability",
        "📝 SAR Reports (RAG)",
        "📈 MLflow Tracking",
        "🗑️ Unlearning"
    ])
    
    # ===== TAB 1: SIGNAL DETECTION =====
    with tab1:
        st.header("📊 Statistical Signal Detection")
        
        # Determine if signals are already cached
        signals_cached = st.session_state.signals is not None
        
        if signals_cached:
            # Signals already loaded - display immediately
            st.success(f"✅ {len(st.session_state.signals):,} signals loaded from cache")
            
        else:
            # Raw ICSR data - show computation button
            st.info("📊 Raw ICSR data loaded. Click 'Run Signal Detection' to compute PRR/ROR/Chi²")
            
            # Local dataset mode selection
            if source_type == "local":
                st.subheader("Local Dataset Analysis Mode")
                col_mode1, col_mode2 = st.columns(2)
                with col_mode1:
                    if st.button("Quantitative (Tier 1)", help="Use minimal PRR/Chi² features only"):
                        st.session_state.analysis_mode = "tier1"
                with col_mode2:
                    if st.button("Qualitative (Tier 2)", help="Use full engineered feature set"):
                        st.session_state.analysis_mode = "tier2"
                st.info(f"Current mode: {'Tier 1 (Quantitative)' if st.session_state.analysis_mode == 'tier1' else 'Tier 2 (Qualitative)'}")
            
            # Signal detection button
            if st.button("🔬 Run Signal Detection", type="primary"):
                with st.spinner("Computing disproportionality signals..."):
                    try:
                        stats_engine = DisproportionalityAnalysis(
                            min_cases=config.MIN_CASE_COUNT
                        )
                        
                        # Compute signals from raw ICSR data
                        signals = stats_engine.compute_signals(st.session_state.raw_data)
                        
                        # Canonicalize
                        normalized_signals = canonicalize_signals(
                            signals,
                            source=source_type,
                            preserve_tier2=True
                        )
                        
                        # Cache in session state
                        st.session_state.signals = normalized_signals
                        
                        st.success(f"✅ Detected {len(normalized_signals):,} signals")
                        st.rerun()
                        
                    except Exception as e:
                        st.error(f"Signal detection failed: {str(e)}")
                        logger.error(f"Signal detection error: {e}", exc_info=True)
        
        # Display signals if available
        if st.session_state.signals is not None:
            df = st.session_state.signals
            
            # Metrics
            col1, col2, col3 = st.columns(3)
            total_pairs = len(df)
            signals_detected = int(df["label"].sum())
            
            col1.metric("Total Drug-Event Pairs", f"{total_pairs:,}")
            col2.metric("Signals Detected", f"{signals_detected:,}")
            if total_pairs > 0:
                rate = (signals_detected / total_pairs) * 100
                col3.metric("Signal Rate", f"{rate:.1f}%")
            
            # Drug Portfolio Filter
            with st.sidebar.expander("🎯 Filter by Drug Portfolio"):
                st.write("**Targeted Surveillance**")
                st.write("Upload Excel with your portfolio drugs to focus analysis")
                
                uploaded_drug_list = st.file_uploader(
                    "Upload drug list (Excel)",
                    type=["xlsx", "xls"],
                    help="Excel file with drug names (any column)",
                    key="drug_filter_upload"
                )
                
                if uploaded_drug_list:
                    from src.utils.drug_filter import DrugFilter
                    
                    if "drug_filter" not in st.session_state:
                        st.session_state.drug_filter = DrugFilter()
                    
                    drug_count = st.session_state.drug_filter.load_from_excel(uploaded_drug_list)
                    
                    if drug_count > 0:
                        st.success(f"✅ Loaded {drug_count} portfolio drugs")
                        
                        with st.expander("View Portfolio"):
                            drugs = st.session_state.drug_filter.drug_list
                            st.write(drugs[:20])
                            if len(drugs) > 20:
                                st.write(f"... and {len(drugs)-20} more")
                        
                        if st.session_state.signals is not None and "drug_filter" in st.session_state:
                            if st.session_state.drug_filter.active:
                                original_count = len(st.session_state.signals)
                                filtered_signals = st.session_state.drug_filter.filter_signals(
                                    st.session_state.signals
                                )
                                
                                st.info(f"🎯 Portfolio Filter Active: {len(filtered_signals):,} signals from {original_count:,} total")
                                
                                # Coverage report
                                coverage = st.session_state.drug_filter.get_coverage_report(
                                    st.session_state.signals
                                )
                                
                                col1, col2, col3 = st.columns(3)
                                col1.metric("Portfolio Drugs", coverage["total_portfolio"])
                                col2.metric("Drugs with Signals", coverage["drugs_with_signals"])
                                col3.metric("Drugs without Signals", coverage["drugs_without_signals"])
                                
                                if coverage["drugs_without_signals"] > 0:
                                    with st.expander(f"View {coverage['drugs_without_signals']} drugs without signals"):
                                        st.write(coverage["missing_list"])
                                
                                # Use filtered signals for display
                                st.session_state.signals = filtered_signals
            
            # Signal filtering controls
            st.subheader("Top Signals (Canonical Schema)")
            
            col1, col2 = st.columns(2)
            with col1:
                prr_threshold = st.slider("PRR Threshold", 1.0, 10.0, 2.0, 0.1)
            with col2:
                min_cases = st.slider("Min Case Count", 3, 50, 3)
            
            # Filter signals
            filtered = df[(df["prr"] >= prr_threshold) & (df["case_count"] >= min_cases)].copy()
            
            # Display columns
            display_cols = ["drug_name", "event_term", "case_count", "prr", "ror", "chi2", "label"]
            extra_cols = [c for c in filtered.columns if c not in display_cols]
            display_cols += extra_cols[:5]  # Show a few extra Tier 2 features
            
            st.dataframe(filtered[display_cols], use_container_width=True)
            
            # Download
            csv = filtered.to_csv(index=False)
            st.download_button(
                "📥 Download Signals CSV",
                csv,
                "signals_normalized.csv",
                "text/csv"
            )
    
    # ===== TAB 2: ML VALIDATION =====
    with tab2:
        st.header("🤖 ML-Based Signal Validation (SISA)")
        
        st.info("""
        **SISA (Sharded, Isolated, Sliced, Aggregated) Training**
        - Tier 1: Baseline model on canonical features (case_count, PRR, Chi²)
        - Tier 2: Extended model using all available engineered features
        """)
        
        if st.session_state.signals is not None:
            col1, col2 = st.columns(2)
            with col1:
                enable_sisa = st.checkbox("Enable SISA Training", value=True)
            with col2:
                num_shards = st.number_input(
                    "Number of Shards",
                    5, 20,
                    config.NUM_SHARDS
                )
            
            if st.button("🤖 Train ML Model", type="primary"):
                with st.spinner("Training SISA ensemble..."):
                    try:
                        signals_df = st.session_state.signals
                        model_dir = Path("models/sisa")
                        
                        trainer = SISATrainer(model_dir)
                        
                        # Determine feature mode
                        mode = st.session_state.analysis_mode if source_type == "local" else "auto"
                        
                        results = trainer.train(
                            signals_df,
                            n_shards=num_shards,
                            mode=mode
                        )
                        
                        st.session_state.trained_model = trainer
                        st.session_state.training_results = results
                        
                        st.success("✅ SISA Model Trained Successfully!")
                        
                        # Display metrics
                        col1, col2, col3 = st.columns(3)
                        with col1:
                            st.metric("AUC Score", f"{results['auc']:.4f}")
                        with col2:
                            st.metric("Training Samples", results["n_train"])
                        with col3:
                            st.metric("Features", results["n_features"])
                        
                        with st.expander("Training Details"):
                            st.json({
                                "auc": results["auc"],
                                "n_features": results["n_features"],
                                "feature_names": results["feature_names"],
                                "n_train": results["n_train"],
                                "n_test": results["n_test"],
                                "n_shards": results["n_shards"],
                            })
                            
                            summary_text = f"""
AUC: {results['auc']:.4f}
Features: {results['n_features']}
Train samples: {results['n_train']}
Test samples: {results['n_test']}
Shards: {results['n_shards']}
"""
                            st.code(summary_text)
                    
                    except Exception as e:
                        st.error(f"Training failed: {str(e)}")
                        logger.error(f"Training error: {e}", exc_info=True)
        else:
            st.warning("⚠️ Run signal detection first (Tab 1)")
    
    # ===== TAB 3: EXPLAINABILITY =====
    with tab3:
        st.header("💡 Model Explainability (SHAP)")
        
        if st.session_state.trained_model is not None:
            if st.button("💡 Generate SHAP Explanations"):
                with st.spinner("Computing SHAP values..."):
                    try:
                        signals_df = st.session_state.signals
                        trainer = st.session_state.trained_model
                        
                        # Get feature names from trainer
                        feature_names = trainer.meta.get("feature_names", [])
                        if not feature_names:
                            st.error("No feature names in trainer meta. Retrain the model.")
                            st.stop()
                        
                        # Build feature matrix
                        X = signals_df[feature_names].astype(float)
                        
                        st.info(f"Using {len(feature_names)} features: {feature_names}")
                        
                        # Initialize SHAP analyzer
                        analyzer = SHAPAnalyzer()
                        
                        # Get a representative model
                        if hasattr(trainer, "shard_models") and trainer.shard_models:
                            model = trainer.shard_models[0]
                        else:
                            st.error("No shard models available for SHAP.")
                            st.stop()
                        
                        # Compute SHAP values
                        analyzer.compute_shap_values(model, X)
                        summary = analyzer.generate_global_summary()
                        
                        st.success("✅ SHAP analysis complete!")
                        
                        # Display results
                        st.subheader("Global Feature Importance")
                        st.image(summary["summary_plot"])
                        
                        st.subheader("Top Features")
                        top_features = pd.DataFrame(summary["top_features"])
                        st.dataframe(top_features, use_container_width=True)
                        
                        st.info(f"Analysis based on {summary['n_samples_used']} samples")
                    
                    except Exception as e:
                        st.error(f"SHAP computation failed: {str(e)}")
                        with st.expander("Error Details"):
                            st.code(str(e))
        else:
            st.warning("⚠️ Train ML model first (Tab 2)")
    
    # ===== TAB 4: SAR REPORTS =====
    with tab4:
        st.header("📝 Signal Assessment Reports (RAG-based)")
        
        st.info("Generate SAR reports using canonical metrics (PRR, Chi², case counts)")
        
        if st.session_state.signals is not None:
            top_signals = st.session_state.signals.head(50).copy()
            
            signal_options = top_signals.apply(
                lambda x: f"{x['drug_name']} - {x['event_term']} (PRR: {x['prr']:.2f}, Cases: {int(x['case_count'])})",
                axis=1
            ).tolist()
            
            selected_idx = st.selectbox(
                "Select Signal for SAR",
                range(len(signal_options)),
                format_func=lambda i: signal_options[i]
            )
            
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
            
            if st.button("📝 Generate SAR Report", type="primary"):
                with st.spinner("Generating comprehensive SAR..."):
                    try:
                        from src.rag.sar_generator import SARGenerator
                        from datetime import datetime
                        
                        generator = SARGenerator(model_name=model_choice)
                        
                        # Check service
                        service_ok, msg = generator.check_service()
                        if not service_ok:
                            st.error(f"{msg}")
                            st.info("Start Ollama: `ollama serve`")
                            st.stop()
                        
                        # Prepare signal data
                        signal_data = {
                            "drug_name": selected_signal["drug_name"],
                            "event_name": selected_signal["event_term"],
                            "count": int(selected_signal["case_count"]),
                            "prr": float(selected_signal["prr"]),
                            "chi_square": float(selected_signal["chi2"]),
                        }
                        
                        # Generate SAR
                        sar_result = generator.generate_sar(signal_data, model_choice)
                        
                        # Cache results
                        st.session_state.current_sar = sar_result
                        st.session_state.current_signal_data = signal_data
                        st.session_state.current_model_choice = model_choice
                        
                        st.success("✅ SAR Generated Successfully!")
                        st.rerun()
                    
                    except Exception as e:
                        st.error(f"SAR generation failed: {str(e)}")
                        with st.expander("Error Details"):
                            st.code(str(e))
            
            # Display cached SAR
            if hasattr(st.session_state, "current_sar") and st.session_state.current_sar is not None:
                from datetime import datetime
                
                st.markdown("---")
                
                sar_result = st.session_state.current_sar
                signal_data = st.session_state.current_signal_data
                model_choice = st.session_state.current_model_choice
                
                st.subheader("Generated SAR Report")
                st.markdown(sar_result["report"])
                
                with st.expander("Evidence Sources"):
                    for i, source in enumerate(sar_result["evidence_sources"], 1):
                        st.write(f"{i}. {source}")
                
                # Download
                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                clean_drug = signal_data["drug_name"].replace(" ", "_").replace(",", "")[:50]
                clean_event = signal_data["event_name"].replace(" ", "_").replace(",", "")[:50]
                report_filename = f"SAR_{clean_drug}_{clean_event}_{timestamp}.txt"
                
                full_report = f"""Signal Assessment Report (SAR)
=====================================
Drug: {signal_data['drug_name']}
Event: {signal_data['event_name']}
Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")} IST

Signal Metrics:
- PRR: {signal_data['prr']:.2f}
- Case Count: {signal_data['count']}
- Chi-Square: {signal_data['chi_square']:.2f}
- LLM Model: {model_choice}

---

{sar_result['report']}

---
"""
                
                st.markdown("---")
                col1, col2 = st.columns([3, 1])
                with col1:
                    st.download_button(
                        label="📥 Download SAR Report",
                        data=full_report,
                        file_name=report_filename,
                        mime="text/plain",
                        use_container_width=True
                    )
                with col2:
                    if st.button("🗑️ Clear", use_container_width=True):
                        st.session_state.current_sar = None
                        st.session_state.current_signal_data = None
                        st.session_state.current_model_choice = None
                        st.rerun()
        
        elif st.session_state.signals is None:
            st.warning("⚠️ Run signal detection first (Tab 1)")
    
    # ===== TAB 5: MLFLOW TRACKING =====
    with tab5:
        st.header("📈 MLflow Experiment Tracking")
        st.info("Track model training experiments, metrics, and artifacts")
        
        try:
            mlflow.set_tracking_uri(config.MLFLOW_TRACKING_URI)
            experiments = mlflow.search_experiments()
            
            if experiments:
                st.subheader("Experiments")
                exp_df = pd.DataFrame([
                    {
                        "Name": exp.name,
                        "ID": exp.experiment_id,
                        "Artifact Location": exp.artifact_location,
                    }
                    for exp in experiments
                ])
                st.dataframe(exp_df, use_container_width=True)
                
                # Select experiment
                exp_names = [exp.name for exp in experiments]
                selected_exp = st.selectbox("Select Experiment", exp_names)
                exp_id = [exp.experiment_id for exp in experiments if exp.name == selected_exp][0]
                
                # Get runs
                runs = mlflow.search_runs(experiment_ids=[exp_id])
                
                if not runs.empty:
                    st.subheader(f"Runs for {selected_exp}")
                    
                    display_cols = ["run_id", "start_time", "status"] + [
                        col for col in runs.columns 
                        if col.startswith("metrics.") or col.startswith("params.")
                    ]
                    st.dataframe(runs[display_cols], use_container_width=True)
                    
                    # Run details
                    if len(runs) > 0:
                        run_idx = st.selectbox(
                            "Select Run",
                            range(len(runs)),
                            format_func=lambda i: f"Run {i+1}: {runs.iloc[i]['run_id'][:8]}"
                        )
                        
                        selected_run = runs.iloc[run_idx]
                        
                        with st.expander("Run Details", expanded=True):
                            col1, col2 = st.columns(2)
                            with col1:
                                st.write("**Parameters**")
                                params = {
                                    k.replace("params.", ""): v
                                    for k, v in selected_run.items()
                                    if k.startswith("params.") and pd.notna(v)
                                }
                                st.json(params)
                            with col2:
                                st.write("**Metrics**")
                                metrics = {
                                    k.replace("metrics.", ""): v
                                    for k, v in selected_run.items()
                                    if k.startswith("metrics.") and pd.notna(v)
                                }
                                st.json(metrics)
                else:
                    st.info("No runs found for this experiment")
            else:
                st.info("No experiments found. Train a model to create experiments.")
        
        except Exception as e:
            st.error(f"MLflow connection error: {str(e)}")
            st.info("Make sure MLflow tracking server is running or using local file store")
    
    # ===== TAB 6: UNLEARNING =====
    with tab6:
        st.header("🗑️ Machine Unlearning (GDPR Right to be Forgotten)")
        
        st.markdown("""
**SISA Unlearning Process:**
1. Enter Case ID to remove
2. System identifies affected data shard
3. Retrain only affected shard (from checkpoint)
4. Update ensemble without full retraining
""")
        
        if st.session_state.trained_model is not None:
            st.subheader("Current Model Version")
            col1, col2 = st.columns(2)
            with col1:
                st.metric("Model Version", "1")
            with col2:
                if hasattr(st.session_state.trained_model, "shard_models"):
                    st.metric("Active Shards", len(st.session_state.trained_model.shard_models))
            
            st.subheader("Unlearn Case")
            case_id = st.text_input("Enter Case ID to Unlearn")
            
            if st.button("🔥 Unlearn Case", type="primary"):
                if case_id:
                    with st.spinner(f"Unlearning case {case_id}..."):
                        try:
                            result = st.session_state.trained_model.unlearn(case_id)
                            
                            if result["status"] == "unlearned":
                                st.success(f"✅ Case {case_id} successfully unlearned!")
                                st.json(result)
                            else:
                                st.warning(f"⚠️ Case {case_id} not found in training data")
                        
                        except Exception as e:
                            st.error(f"Unlearning failed: {str(e)}")
                else:
                    st.error("Please enter a case ID")
        else:
            st.warning("⚠️ Train SISA model first (Tab 2)")

else:
    st.info("👈 Select a data source from the sidebar to begin")

# ===== SIDEBAR FOOTER =====
st.sidebar.markdown("---")
st.sidebar.markdown("### ℹ️ System Info")
st.sidebar.markdown(f"**Version:** {config.VALIDATION_STATUS['code_version']}")
st.sidebar.markdown(f"**MLflow URI:** {config.MLFLOW_TRACKING_URI}")

try:
    from src.rag.sar_generator import SARGenerator
    sar_gen = SARGenerator()
    model_ok, model_msg = sar_gen.check_model()
    st.sidebar.markdown(f"**Ollama Model:** {sar_gen.model_name}")
except Exception:
    pass
