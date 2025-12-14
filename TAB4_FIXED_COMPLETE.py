# ============================================
# TAB 4: SAR REPORTS - FIXED VERSION
# Save button now works correctly!
# ============================================

with tab4:
    st.header("Signal Assessment Reports (RAG-based)")
    st.info("Generate comprehensive SAR reports with causality assessment and literature evidence")
    
    # Initialize session state for SAR
    if 'current_sar' not in st.session_state:
        st.session_state.current_sar = None
    if 'current_signal_data' not in st.session_state:
        st.session_state.current_signal_data = None
    if 'current_model_choice' not in st.session_state:
        st.session_state.current_model_choice = None
    
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
        
        # GENERATE BUTTON (stores result in session_state)
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
                    st.rerun()  # Rerun to display the SAR
                    
                except Exception as e:
                    st.error(f"❌ SAR generation failed: {str(e)}")
                    with st.expander("Error Details"):
                        st.code(str(e))
        
        # ===== DISPLAY SAR (if available in session_state) =====
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
            
            # ===== SAVE & DOWNLOAD BUTTONS (OUTSIDE Generate block!) =====
            st.markdown("---")
            from datetime import datetime
            
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            clean_drug = signal_data['drug_name'].replace(' ', '_').replace('/', '_')[:50]
            clean_event = signal_data['event_name'].replace(' ', '_').replace('/', '_')[:50]
            report_filename = f"SAR_{clean_drug}_{clean_event}_{timestamp}.md"
            
            # Format complete report
            full_report = f"""# Signal Assessment Report (SAR)

## Drug-Event Pair
**Drug:** {signal_data['drug_name']}  
**Event:** {signal_data['event_name']}  
**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S IST')}

## Signal Strength Metrics
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
            
            # Save locally for backup
            from pathlib import Path
            reports_dir = Path("outputs/sar_reports")
            reports_dir.mkdir(parents=True, exist_ok=True)
            local_report_path = reports_dir / report_filename
            
            try:
                with open(local_report_path, 'w', encoding='utf-8') as f:
                    f.write(full_report)
                save_success = True
            except Exception as e:
                save_success = False
                st.error(f"Failed to save locally: {str(e)}")
            
            # Display buttons
            col1, col2, col3 = st.columns([2, 2, 1])
            
            with col1:
                st.download_button(
                    label="📥 Download SAR Report",
                    data=full_report,
                    file_name=report_filename,
                    mime="text/markdown",
                    use_container_width=True,
                    help="Download to your Downloads folder"
                )
            
            with col2:
                if save_success:
                    st.success(f"✅ Saved: {report_filename}")
                else:
                    st.error("❌ Local save failed")
            
            with col3:
                if st.button("🗑️ Clear", use_container_width=True):
                    st.session_state.current_sar = None
                    st.session_state.current_signal_data = None
                    st.session_state.current_model_choice = None
                    st.rerun()
            
            # Log to MLflow
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
                    report_path=str(local_report_path) if save_success else None
                )
                st.info(f"📊 Logged to MLflow: Run ID {run_id[:8]}")
            except Exception as e:
                st.warning(f"⚠️ MLflow logging failed: {str(e)}")
    
    else:
        st.warning("⚠️ Run signal detection first (Tab 1)")
