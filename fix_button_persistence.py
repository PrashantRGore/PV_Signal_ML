from pathlib import Path

dsm_file = Path('src/data/data_source_manager.py')
content = dsm_file.read_text(encoding='utf-8')

# Fix the button persistence issue using session_state
old_config_section = '''        # Configure button
        configure_dates = st.sidebar.button("📅 Configure Date Range", help="Click after selecting dates")
        
        st.sidebar.write(f"DEBUG: Button clicked = {load_faers}")
        st.sidebar.write(f"DEBUG: Selected quarters = {selected_quarters}")
        
        if not configure_dates:
            st.sidebar.info("👆 Click 'Configure Date Range' to proceed")
            return None'''

new_config_section = '''        # Initialize session state
        if 'dates_configured' not in st.session_state:
            st.session_state.dates_configured = False
        
        # Configure button
        if st.sidebar.button("📅 Configure Date Range", help="Click after selecting dates"):
            st.session_state.dates_configured = True
        
        if not st.session_state.dates_configured:
            st.sidebar.info("👆 Click 'Configure Date Range' to proceed")
            return None'''

content = content.replace(old_config_section, new_config_section)
dsm_file.write_text(content, encoding='utf-8')
print('✅ Fixed button persistence with session_state!')
