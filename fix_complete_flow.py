from pathlib import Path

dsm_file = Path('src/data/data_source_manager.py')
content = dsm_file.read_text(encoding='utf-8')

# Find and replace the entire FAERS section with session state management
old_faers_section = '''        # Initialize session state
        if 'dates_configured' not in st.session_state:
            st.session_state.dates_configured = False
        
        # Configure button
        if st.sidebar.button("📅 Configure Date Range", help="Click after selecting dates"):
            st.session_state.dates_configured = True
        
        if not st.session_state.dates_configured:
            st.sidebar.info("👆 Click 'Configure Date Range' to proceed")
            return None
        
        # Compute quarters after validation
        with st.spinner("🔍 Computing quarters..."):
            quarters = self._compute_quarters(start_date, end_date)
        
        st.sidebar.success(f"✅ Detected **{len(quarters)} quarters**: {', '.join([f'Q{q[1]}-{q[0]}' for q in quarters])}")
        
        selected_quarters = st.sidebar.multiselect(
            "Confirm quarters to download:",
            options=quarters,
            default=quarters,
            format_func=lambda q: f"Q{q[1]}-{q[0]}"
        )
        
        load_faers = st.sidebar.button("📡 Download & Process FAERS Data")
        
        st.sidebar.write(f"DEBUG: Button clicked = {load_faers}")
        st.sidebar.write(f"DEBUG: Selected quarters = {selected_quarters}")
        
        if load_faers and selected_quarters:
            st.sidebar.warning("🚀 DOWNLOAD STARTING...")
            try:
                downloader = FAERSDownloader()'''

new_faers_section = '''        # Initialize session state for date configuration
        if 'faers_dates_configured' not in st.session_state:
            st.session_state.faers_dates_configured = False
            st.session_state.faers_quarters = []
        
        # Configure button - store dates and quarters in session
        if st.sidebar.button("📅 Configure Date Range", help="Click after selecting dates"):
            st.session_state.faers_dates_configured = True
            st.session_state.faers_start_date = start_date
            st.session_state.faers_end_date = end_date
            # Compute quarters immediately
            st.session_state.faers_quarters = self._compute_quarters(start_date, end_date)
        
        # Show configuration status
        if not st.session_state.faers_dates_configured:
            st.sidebar.info("👆 Click 'Configure Date Range' to proceed")
            return None
        
        # Show detected quarters
        quarters = st.session_state.faers_quarters
        st.sidebar.success(f"✅ Detected **{len(quarters)} quarters**")
        st.sidebar.write(", ".join([f"Q{q[1]}-{q[0]}" for q in quarters]))
        
        # Quarter selection
        selected_quarters = st.sidebar.multiselect(
            "Select quarters to download:",
            options=quarters,
            default=quarters,
            format_func=lambda q: f"Q{q[1]}-{q[0]}"
        )
        
        # Download button
        if st.sidebar.button("📡 Download & Process FAERS Data", type="primary"):
            if not selected_quarters:
                st.sidebar.error("⚠️ Please select at least one quarter")
                return None
            
            st.sidebar.warning(f"🚀 STARTING DOWNLOAD OF {len(selected_quarters)} QUARTERS...")
            
            try:
                downloader = FAERSDownloader()'''

content = content.replace(old_faers_section, new_faers_section)
dsm_file.write_text(content, encoding='utf-8')
print('✅ Fixed with proper session state - quarters will persist!')
