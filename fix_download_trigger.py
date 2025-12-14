from pathlib import Path

dsm_file = Path('src/data/data_source_manager.py')
content = dsm_file.read_text(encoding='utf-8')

# Fix: Remove debug lines and ensure download persists
old_download = '''        # Download button
        if st.sidebar.button("📡 Download & Process FAERS Data", type="primary"):
            if not selected_quarters:
                st.sidebar.error("⚠️ Please select at least one quarter")
                return None
            
            st.sidebar.warning(f"🚀 STARTING DOWNLOAD OF {len(selected_quarters)} QUARTERS...")
            
            try:
                downloader = FAERSDownloader()
                with st.spinner(f"Downloading {len(selected_quarters)} quarters..."):
                    data = downloader.download_and_process(
                        quarters=selected_quarters,
                        start_date=start_date,
                        end_date=end_date
                    )'''

new_download = '''        # Download button - process immediately when clicked
        download_clicked = st.sidebar.button("📡 Download & Process FAERS Data", type="primary")
        
        if download_clicked:
            if not selected_quarters:
                st.sidebar.error("⚠️ Please select at least one quarter")
                return None
            
            st.sidebar.warning(f"🚀 STARTING DOWNLOAD OF {len(selected_quarters)} QUARTERS...")
            st.sidebar.write("📥 Connecting to FDA FAERS database...")
            
            try:
                # Create downloader and start download immediately
                from src.data.faers_downloader import FAERSDownloader
                downloader = FAERSDownloader()
                
                # Download with progress (this will show progress bars)
                data = downloader.download_and_process(
                    quarters=selected_quarters,
                    start_date=start_date,
                    end_date=end_date
                )'''

content = content.replace(old_download, new_download)
dsm_file.write_text(content, encoding='utf-8')
print('✅ Fixed download trigger - progress will now be visible!')
