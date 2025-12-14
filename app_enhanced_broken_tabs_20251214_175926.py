import streamlit as st
import pandas as pd
import webbrowser
from pathlib import Path
from datetime import datetime
import os
import config
from src.data.data_source_manager import DataSourceManager
from src.stats_engine.disproportionality import DisproportionalityAnalysis
from src.ml.sisa_trainer import SISATrainer
from src.explainability.shap_analysis import SHAPAnalyzer
from src.rag.sar_generator import SARGenerator
from src.mlflow_tracking.tracker import PVMLflowTracker
from src.utils.service_health import ServiceHealthChecker
from src.utils.logger import setup_logger

logger = setup_logger(__name__)

# ===== PAGE CONFIG =====
st.set_page_config(
    page_title='PV Signal ML - Enhanced',
    page_icon='🔬',
    layout='wide',
    initial_sidebar_state='expanded'
)

# ===== SESSION STATE INITIALIZATION =====
if 'data_loaded' not in st.session_state:
    st.session_state.data_loaded = False
if 'signals' not in st.session_state:
    st.session_state.signals = None
if 'trained_model' not in st.session_state:
    st.session_state.trained_model = None
if 'model_version' not in st.session_state:
    st.session_state.model_version = 1
if 'mlflow_tracker' not in st.session_state:
    st.session_state.mlflow_tracker = PVMLflowTracker()
if 'service_checker' not in st.session_state:
    st.session_state.service_checker = ServiceHealthChecker()

# ===== DISCLAIMER =====
st.warning(config.DISCLAIMER_TEXT)

# ===== TITLE =====
st.title('🔬 Pharmacovigilance Signal Detection with SISA & SHAP')
st.markdown('**Enhanced System with Machine Unlearning, Explainable AI, RAG-based SAR, and MLflow**')

# ===== SIDEBAR: SERVICE STATUS =====
st.sidebar.markdown('### 🔧 Service Status')
service_status = st.session_state.service_checker.get_service_status()

col1, col2 = st.sidebar.columns(2)
col1.markdown(f"**Ollama:** {'🟢' if service_status['ollama'] else '🔴'}")
col2.markdown(f"**MLflow:** {'🟢' if service_status['mlflow'] else '🔴'}")

if not service_status['ollama']:
    if st.sidebar.button('🚀 Start Ollama', use_container_width=True):
        with st.spinner('Starting Ollama...'):
            if st.session_state.service_checker.start_ollama():
                st.success('✅ Ollama started!')
                st.rerun()
            else:
                st.error('❌ Run manually: ollama serve')

if not service_status['mlflow']:
    if st.sidebar.button('🚀 Start MLflow UI', use_container_width=True):
        with st.spinner('Starting MLflow UI...'):
            if st.session_state.service_checker.start_mlflow_ui():
                st.success('✅ MLflow UI started')
                st.rerun()
            else:
                st.error('❌ Run manually: mlflow ui --port 5000')

st.sidebar.markdown('---')

# ===== SIDEBAR: PIPELINE STATUS WITH LABELS =====
st.sidebar.markdown('### 📊 Pipeline Status')
st.sidebar.metric('Model Version', st.session_state.model_version)

pipeline_stages = [
    ('Data Loaded', st.session_state.data_loaded, '📁'),
    ('Signals Detected', st.session_state.signals is not None, '🔍'),
    ('SISA Trained', st.session_state.trained_model is not None, '🤖'),
]

for label, status, icon in pipeline_stages:
    status_text = f"{icon} {label}: {'✅' if status else '⏸️'}"
    st.sidebar.markdown(status_text)

st.sidebar.markdown('---')

# ===== DATA SOURCE MANAGER =====
data_manager = DataSourceManager()
source_type = data_manager.select_data_source()

data = None
if source_type == 'local':
    data = data_manager.load_local_dataset()
elif source_type == 'demo_hf':
    data = data_manager.load_demo_dataset()
elif source_type == 'faers_live':
    data = data_manager.load_faers_dataset()

if data is not None:
    st.session_state.data_loaded = True
    st.session_state.raw_data = data
    st.session_state.metadata = data_manager.get_metadata()
    st.sidebar.success(f'✅ Data loaded: {len(data):,} records')
    with st.sidebar.expander('📋 Data Metadata'):
        st.json(data_manager.get_metadata())

# ===== MAIN TABS =====
if st.session_state.data_loaded:
    tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
        '🔍 Signal Detection',
        '🤖 ML Validation',
        '💡 Explainability',
        '📝 SAR Reports (RAG)',
        '📋 MLflow Tracking',
        '🗑️ Unlearning'
    ])
    
    # ========== TAB 1: SIGNAL DETECTION (COMPLETE IMPLEMENTATION) ==========
    with tab1:
        st.header('Statistical Signal Detection')
        
        if st.button('🚀 Run Signal Detection', type='primary'):
            with st.spinner('Computing disproportionality signals...'):
                try:
                    stats_engine = DisproportionalityAnalysis()
                    signals = stats_engine.compute_signals(st.session_state.raw_data)
                    st.session_state.signals = signals
                    
                    # Try to apply ML causality scoring if available
                    try:
                        from src.ml.causality_scorer import CausalityScorer
                        st.info('🤖 Applying ML causality model to rank signals...')
                        scorer = CausalityScorer()
                        
                        if scorer.load_model():
                            signals = scorer.score_signals(signals)
                            signals = signals.sort_values('causality_score', ascending=False)
                            st.session_state.signals = signals
                            st.success('✅ Causality assessment complete')
                            
                            col1, col2, col3 = st.columns(3)
                            high_causal = (signals['causality_score'] > 0.7).sum()
                            mod_causal = ((signals['causality_score'] >= 0.3) & (signals['causality_score'] <= 0.7)).sum()
                            low_causal = (signals['causality_score'] < 0.3).sum()
                            
                            col1.metric('High Causality', f'{high_causal:,}', delta='Priority')
                            col2.metric('Moderate Causality', f'{mod_causal:,}')
                            col3.metric('Low Causality', f'{low_causal:,}')
                    except:
                        st.info('ℹ️ Using statistical signals only (ML causality model not available)')
                    
                    # Log to MLflow
                    try:
                        st.session_state.mlflow_tracker.log_signal_detection(signals, {
                            'prr_threshold': 2.0,
                            'chi2_threshold': 4.0,
                            'min_cases': 3
                        })
                    except:
                        pass
                    
                    st.success(f'✅ Detected {len(signals)} signals')
                    
                    # Display summary metrics
                    col1, col2, col3 = st.columns(3)
                    col1.metric('Total Signals', f'{len(signals):,}')
                    col2.metric('Avg PRR', f"{signals['prr'].mean():.2f}")
                    col3.metric('Max Chi²', f"{signals['chi2'].max():.2f}")
                
                except Exception as e:
                    st.error(f'Error during signal detection: {str(e)}')
                    import traceback
                    st.code(traceback.format_exc())
        
        # Drug Portfolio Filter
        with st.sidebar.expander('🎯 Filter by Drug Portfolio'):
            st.write('**Targeted Surveillance**')
            st.write('Upload Excel with your portfolio drugs')
            
            uploaded_drug_list = st.file_uploader(
                'Upload drug list (Excel)',
                type=['xlsx', 'xls'],
                help='Excel file with drug names',
                key='drug_filter_upload'
            )
            
            if uploaded_drug_list:
                try:
                    from src.utils.drug_filter import DrugFilter
                    if 'drug_filter' not in st.session_state:
                        st.session_state.drug_filter = DrugFilter()
                    
                    drug_count = st.session_state.drug_filter.load_from_excel(uploaded_drug_list)
                    
                    if drug_count > 0:
                        st.success(f'✅ Loaded {drug_count} portfolio drugs')
                        with st.expander('View Portfolio'):
                            drugs = st.session_state.drug_filter.drug_list
                            st.write(drugs[:20])
                            if len(drugs) > 20:
                                st.write(f'... and {len(drugs)-20} more')
                except Exception as e:
                    st.error(f'Error loading drug filter: {str(e)}')
        
        # Apply drug filter if active
        if st.session_state.signals is not None:
            display_signals = st.session_state.signals.copy()
            
            if 'drug_filter' in st.session_state and st.session_state.drug_filter.active:
                original_count = len(display_signals)
                display_signals = st.session_state.drug_filter.filter_signals(display_signals)
                st.info(f'🎯 **Portfolio Filter Active:** {len(display_signals):,} signals from {original_count:,} total')
                
                # Coverage report
                coverage = st.session_state.drug_filter.get_coverage_report(st.session_state.signals)
                col1, col2, col3 = st.columns(3)
                col1.metric('Portfolio Drugs', coverage['total_portfolio'])
                col2.metric('Drugs with Signals', coverage['drugs_with_signals'])
                col3.metric('Drugs without Signals', coverage['drugs_without_signals'])
                
                if coverage['drugs_without_signals'] > 0:
                    with st.expander(f"View {coverage['drugs_without_signals']} drugs without signals"):
                        st.write(coverage['missing_list'])
            
            # Display signals table
            st.subheader('Detected Signals')
            
            # Filter controls
            col1, col2 = st.columns(2)
            with col1:
                prr_threshold = st.slider('PRR Threshold', 1.0, 10.0, 2.0, 0.1)
            with col2:
                min_cases = st.slider('Min Case Count', 3, 50, 3)
            
            # Apply filters
            filtered = display_signals[
                (display_signals['prr'] >= prr_threshold) &
                (display_signals['case_count'] >= min_cases)
            ]
            
            # Prepare display columns
            display_cols = ['drug_name', 'event_term', 'case_count', 'prr', 'chi2']
            if 'causality_score' in filtered.columns:
                display_cols.extend(['causality_score', 'causality_level'])
            if 'is_signal_prr' in filtered.columns:
                display_cols.append('is_signal_prr')
            
            st.dataframe(filtered[display_cols].head(100), use_container_width=True)
            
            # Download button
            csv = filtered.to_csv(index=False)
            st.download_button(
                '📥 Download Signals CSV',
                csv,
                'signals.csv',
                'text/csv',
                key='download_signals'
            )
        else:
            st.info('👆 Click "Run Signal Detection" to start analysis')
    
    # ========== TAB 2: ML VALIDATION (COMPLETE IMPLEMENTATION) ==========
    with tab2:
        st.header('ML-Based Signal Validation (SISA)')
        st.info('''
        **SISA (Sharded, Isolated, Sliced, Aggregated) Training**
        - Enables efficient machine unlearning for GDPR compliance
        - Trains ensemble of models on data shards
        - Allows removal of individual cases without full retraining
        ''')
        
        if st.session_state.signals is not None:
            col1, col2 = st.columns(2)
            with col1:
                enable_sisa = st.checkbox('Enable SISA Training', value=config.ENABLE_SISA)
            with col2:
                num_shards = st.number_input('Number of Shards', 5, 20, config.NUM_SHARDS)
            
            if st.button('🤖 Train ML Model', type='primary'):
                with st.spinner('Training SISA ensemble...'):
                    try:
                        signals_df = st.session_state.signals
                        
                        if enable_sisa:
                            model_dir = Path('models/sisa')
                            trainer = SISATrainer(model_dir)
                            
                            # Train SISA model
                            results = trainer.train(signals_df, n_shards=num_shards)
                            st.session_state.trained_model = trainer
                            st.session_state.training_results = results
                            
                            st.success('✅ SISA Model Trained Successfully!')
                            
                            # Log to MLflow
                            try:
                                st.session_state.mlflow_tracker.log_sisa_training(results, model=trainer)
                            except:
                                pass
                            
                            # Display metrics
                            col1, col2, col3 = st.columns(3)
                            with col1:
                                st.metric('AUC Score', f"{results['auc']:.4f}")
                            with col2:
                                st.metric('Training Samples', results['n_train'])
                            with col3:
                                st.metric('Features', results['n_features'])
                            
                            # Show detailed results
                            with st.expander('📊 Training Details'):
                                st.json({
                                    'auc': results['auc'],
                                    'n_features': results['n_features'],
                                    'n_train': results['n_train'],
                                    'n_test': results['n_test'],
                                    'n_shards': num_shards
                                })
                                st.code(results['report'])
                        else:
                            st.warning('SISA disabled - using standard training')
                    
                    except Exception as e:
                        st.error(f'Training error: {str(e)}')
                        import traceback
                        st.code(traceback.format_exc())
            
            # Show current model status
            if st.session_state.trained_model is not None:
                st.success('✅ Model is trained and ready')
                st.markdown('**Next Step:** Go to "Explainability" tab to generate SHAP explanations')
        else:
            st.warning('⚠️ Run signal detection first (Tab 1)')
    
    # ========== TAB 3: EXPLAINABILITY (COMPLETE IMPLEMENTATION) ==========
    with tab3:
        st.header('Model Explainability (SHAP)')
        
        if st.session_state.trained_model is not None:
            if st.button('💡 Generate SHAP Explanations'):
                with st.spinner('Computing SHAP values...'):
                    try:
                        signals_df = st.session_state.signals
                        X = signals_df[['case_count', 'prr', 'ror', 'chi2']].fillna(0)
                        feature_names = ['case_count', 'prr', 'ror', 'chi2']
                        
                        analyzer = SHAPAnalyzer()
                        
                        # Get model
                        if hasattr(st.session_state.trained_model, 'model'):
                            model = st.session_state.trained_model.model
                        elif hasattr(st.session_state.trained_model, 'shard_models'):
                            model = st.session_state.trained_model.shard_models[0]
                        else:
                            st.error('Model structure not compatible')
                            st.stop()
                        
                        # Compute SHAP
                        analyzer.compute_shap_values(model, X, feature_names=feature_names)
                        summary = analyzer.generate_global_summary()
                        
                        st.success('✅ SHAP analysis complete!')
                        
                        # Log to MLflow
                        try:
                            st.session_state.mlflow_tracker.log_shap_analysis(summary)
                        except:
                            pass
                        
                        # Display global summary
                        st.subheader('Global Feature Importance')
                        st.image(summary['summary_plot'])
                        
                        # Display top features
                        st.subheader('Top Features')
                        top_features_df = pd.DataFrame(summary['top_features'])
                        st.dataframe(top_features_df, use_container_width=True)
                    
                    except Exception as e:
                        st.error(f'SHAP error: {str(e)}')
                        st.info('Ensure model is XGBoost-compatible and data is properly formatted')
                        import traceback
                        st.code(traceback.format_exc())
        else:
            st.warning('⚠️ Train ML model first (Tab 2)')
    
    # ========== TAB 4: REPORTS (PLACEHOLDER) ==========
    with tab4:
        st.header('Signal Assessment Reports')
        
        if st.session_state.signals is not None:
            st.markdown('**Generate Summary Reports**')
            
            if st.button('📄 Generate Signal Summary'):
                st.info('Basic report generation')
                st.markdown('For advanced RAG-based SAR reports, use the next tab')
        else:
            st.warning('⚠️ No signals available')
    
    # ========== TAB 5: SAR REPORTS (COMPLETE IMPLEMENTATION) ==========
    with tab4:
        st.header('📝 Signal Assessment Report (SAR) Generation')
        st.markdown('*Generate regulatory-compliant SAR reports using RAG and Ollama LLMs*')
        
        # Service status check
        ollama_status = st.session_state.service_checker.check_ollama()
        
        if not ollama_status:
            st.error('🔴 Ollama service not running!')
            st.info('Click "Start Ollama" in the sidebar or run: ollama serve')
        else:
            st.success('🟢 Ollama service running')
        
        if st.session_state.signals is not None:
            signals_df = st.session_state.signals
            
            col1, col2 = st.columns([2, 1])
            
            with col1:
                signal_options = [
                    f"{row['drug_name']} → {row['event_term']} (PRR={row['prr']:.2f})"
                    for idx, row in signals_df.head(50).iterrows()
                ]
                selected_signal = st.selectbox('Select Signal for SAR:', signal_options)
            
            with col2:
                model_choice = st.selectbox(
                    'LLM Model:',
                    ['llama3.2', 'mistral', 'llama2'],
                    help='Ollama model for report generation'
                )
            
            st.info(f'📊 Using Model Version {st.session_state.model_version}')
            
            # Check if model is available
            if ollama_status:
                model_available = st.session_state.service_checker.check_ollama_model(model_choice)
                if not model_available:
                    st.warning(f'⚠️ Model {model_choice} not found')
                    if st.button(f'📥 Pull {model_choice} Model'):
                        with st.spinner(f'Pulling {model_choice}...'):
                            if st.session_state.service_checker.pull_ollama_model(model_choice):
                                st.success(f'✅ {model_choice} ready!')
                                st.rerun()
                            else:
                                st.error(f'❌ Failed. Run: ollama pull {model_choice}')
            
            if st.button('🔬 Generate SAR Report', type='primary', disabled=not ollama_status):
                idx = signal_options.index(selected_signal)
                signal_row = signals_df.iloc[idx]
                
                signal_data = {
                    'drug_name': signal_row['drug_name'],
                    'event_name': signal_row['event_term'],
                    'prr': signal_row['prr'],
                    'chi_square': signal_row.get('chi2', 0),
                    'count': signal_row['case_count']
                }
                
                with st.spinner(f'Generating SAR using {model_choice}...'):
                    try:
                        sar_gen = SARGenerator(model_name=model_choice)
                        sar_result = sar_gen.generate_sar(signal_data, model_choice)
                        filepath = sar_gen.save_report(sar_result)
                        
                        sar_result['model_version'] = st.session_state.model_version
                        
                        # Log to MLflow
                        try:
                            run_id = st.session_state.mlflow_tracker.log_sar_generation(sar_result)
                            st.success(f'✅ SAR Generated! (MLflow Run: {run_id[:8]}...)')
                        except:
                            st.success('✅ SAR Generated!')
                        
                        st.markdown('### Generated SAR Report')
                        st.markdown(f'**Model Version:** {st.session_state.model_version}')
                        st.text_area('Report Content:', sar_result['report'], height=400)
                        
                        if sar_result['evidence_sources']:
                            st.markdown('### Evidence Sources')
                            for i, source in enumerate(sar_result['evidence_sources'], 1):
                                st.markdown(f'{i}. {source}')
                        
                        with open(filepath, 'r', encoding='utf-8') as f:
                            st.download_button(
                                '📥 Download SAR Report',
                                f.read(),
                                file_name=os.path.basename(filepath),
                                mime='text/plain'
                            )
                    
                    except Exception as e:
                        st.error(f'Error generating SAR: {str(e)}')
                        if 'not running' in str(e):
                            st.info('Start Ollama using sidebar button')
                        elif 'not found' in str(e):
                            st.info(f'Pull model: ollama pull {model_choice}')
        else:
            st.info('⚠️ Run signal detection first')
    
    # ========== TAB 6: MLFLOW TRACKING (COMPLETE IMPLEMENTATION) ==========
    with tab4:
        st.header('📋 MLflow Experiment Tracking')
        st.markdown('*Audit trail for regulatory compliance*')
        
        mlflow_running = st.session_state.service_checker.check_mlflow()
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            if st.button('🔄 Refresh Runs'):
                st.rerun()
        
        with col2:
            if not mlflow_running:
                if st.button('🚀 Start MLflow UI'):
                    with st.spinner('Starting MLflow UI...'):
                        if st.session_state.service_checker.start_mlflow_ui():
                            st.success('✅ MLflow UI started!')
                            st.rerun()
                        else:
                            st.error('Run: mlflow ui --port 5000')
        
        with col3:
            if mlflow_running:
                if st.button('🌐 Open MLflow UI'):
                    webbrowser.open('http://localhost:5000')
                    st.success('Opening MLflow UI...')
            else:
                st.warning('MLflow UI not running')
        
        # Display runs
        try:
            runs = st.session_state.mlflow_tracker.get_all_runs()
            
            if len(runs) > 0:
                st.subheader('Pipeline Runs')
                runs_data = []
                for run in runs:
                    runs_data.append({
                        'Run ID': run.info.run_id[:8],
                        'Name': run.data.tags.get('mlflow.runName', 'N/A')[:40],
                        'Start Time': datetime.fromtimestamp(run.info.start_time/1000).strftime('%Y-%m-%d %H:%M'),
                        'Model Ver': run.data.params.get('model_version', 'N/A'),
                        'Unlearning': '🗑️' if run.data.params.get('unlearning_event') == 'true' else '',
                        'SHAP': '✅' if run.data.tags.get('shap_computed') == 'true' else '',
                        'SISA': '✅' if run.data.tags.get('sisa_trained') == 'true' else ''
                    })
                st.dataframe(pd.DataFrame(runs_data), use_container_width=True)
                
                st.subheader('Model Version History')
                version_history = st.session_state.mlflow_tracker.get_model_version_history()
                if version_history:
                    st.dataframe(pd.DataFrame(version_history), use_container_width=True)
            else:
                st.info('No experiment runs yet')
        
        except Exception as e:
            st.error(f'Error: {str(e)}')
    
    # ========== TAB 7: UNLEARNING (COMPLETE IMPLEMENTATION) ==========
    with tab5:
        st.header('Machine Unlearning (GDPR Right to be Forgotten)')
        st.markdown('''
        **SISA Unlearning Process:**
        1. Enter Case ID to remove
        2. System identifies affected data shard
        3. Retrain only affected shard from checkpoint
        4. Update ensemble without full retraining
        ''')
        
        if st.session_state.trained_model is not None:
            col1, col2 = st.columns([2, 1])
            with col1:
                case_id = st.text_input('Enter Case ID to Unlearn')
            with col2:
                st.metric('Current Model Version', st.session_state.model_version)
            
            if st.button('🗑️ Unlearn Case', type='primary'):
                if case_id:
                    with st.spinner(f'Unlearning case {case_id}...'):
                        try:
                            result = st.session_state.trained_model.unlearn(case_id)
                            
                            if result.get('status') == 'unlearned':
                                # Increment model version
                                st.session_state.model_version += 1
                                
                                # Log to MLflow
                                try:
                                    st.session_state.mlflow_tracker.log_unlearning_event({
                                        'case_id': case_id,
                                        'shards_retrained': result.get('shards', [])
                                    })
                                except:
                                    pass
                                
                                st.success(f'✅ Case {case_id} successfully unlearned!')
                                st.success(f'📊 New Model Version: {st.session_state.model_version}')
                                st.json(result)
                                st.warning('⚠️ Recompute SHAP and regenerate affected SAR reports')
                            else:
                                st.warning(f'Case {case_id} not found in training data')
                        
                        except Exception as e:
                            st.error(f'Unlearning error: {str(e)}')
                else:
                    st.error('Please enter a case ID')
        else:
            st.warning('⚠️ Train SISA model first (Tab 2)')

else:
    st.info('👈 Select a data source from the sidebar to begin')

# ===== FOOTER =====
st.sidebar.markdown('---')
st.sidebar.markdown(f"**Version:** {config.VALIDATION_STATUS['code_version']}")
st.sidebar.markdown(f'**MLflow:** {config.MLFLOW_TRACKING_URI}')


