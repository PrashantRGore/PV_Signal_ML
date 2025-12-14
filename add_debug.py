from pathlib import Path

# Add debugging to see if download is triggered
dsm_file = Path('src/data/data_source_manager.py')
content = dsm_file.read_text(encoding='utf-8')

# Find the download button section
old_button = '''        load_faers = st.sidebar.button("📡 Download & Process FAERS Data")
        
        if load_faers and selected_quarters:
            try:
                downloader = FAERSDownloader()'''

new_button = '''        load_faers = st.sidebar.button("📡 Download & Process FAERS Data")
        
        st.sidebar.write(f"DEBUG: Button clicked = {load_faers}")
        st.sidebar.write(f"DEBUG: Selected quarters = {selected_quarters}")
        
        if load_faers and selected_quarters:
            st.sidebar.warning("🚀 DOWNLOAD STARTING...")
            try:
                downloader = FAERSDownloader()'''

content = content.replace(old_button, new_button)
dsm_file.write_text(content, encoding='utf-8')
print('✅ Added debug output to see button state')
