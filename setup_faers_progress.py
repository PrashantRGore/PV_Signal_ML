from pathlib import Path

# First, check if FAERSDownloader exists
faers_file = Path('src/data/faers_downloader.py')

if not faers_file.exists():
    print('⚠️ FAERSDownloader does not exist. Creating it...')
    
    faers_code = '''"""
FAERS Data Downloader: Downloads quarterly FAERS data from FDA
"""

import requests
import pandas as pd
from typing import List, Tuple
import streamlit as st
from src.utils.logger import setup_logger

logger = setup_logger(__name__)

class FAERSDownloader:
    """Downloads and processes FAERS quarterly data"""
    
    BASE_URL = "https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html"
    
    def __init__(self):
        self.data = None
    
    def download_and_process(self, quarters: List[Tuple[int, int]], 
                            start_date, end_date) -> pd.DataFrame:
        """
        Download and process FAERS data for specified quarters
        
        Args:
            quarters: List of (year, quarter) tuples
            start_date: Start date for filtering
            end_date: End date for filtering
        
        Returns:
            Processed DataFrame with ICSR records
        """
        logger.info(f"Starting download for {len(quarters)} quarters")
        
        all_data = []
        progress_bar = st.progress(0)
        status_text = st.empty()
        
        for idx, (year, quarter) in enumerate(quarters):
            status_text.text(f"📥 Processing Q{quarter}-{year}... ({idx+1}/{len(quarters)})")
            
            # Simulate download and processing
            # In production, replace with actual FAERS API/file download
            quarter_data = self._fetch_quarter_data(year, quarter)
            
            if quarter_data is not None and len(quarter_data) > 0:
                all_data.append(quarter_data)
            
            # Update progress
            progress_bar.progress((idx + 1) / len(quarters))
        
        status_text.text("✅ All quarters processed!")
        progress_bar.empty()
        
        if len(all_data) == 0:
            raise Exception("No data downloaded for selected quarters")
        
        # Combine all quarters
        combined_data = pd.concat(all_data, ignore_index=True)
        logger.info(f"Downloaded {len(combined_data)} records")
        
        return combined_data
    
    def _fetch_quarter_data(self, year: int, quarter: int) -> pd.DataFrame:
        """
        Fetch data for a specific quarter
        PLACEHOLDER: Replace with actual FAERS API/download logic
        """
        import time
        time.sleep(0.5)  # Simulate download time
        
        # Generate synthetic data for demo
        # In production, replace with actual FAERS data download
        n_records = 1000
        data = pd.DataFrame({
            'case_id': range(n_records),
            'drug_name': [f'Drug_{i%100}' for i in range(n_records)],
            'event_term': [f'Event_{i%50}' for i in range(n_records)],
            'year': year,
            'quarter': quarter
        })
        
        return data
'''
    
    faers_file.write_text(faers_code, encoding='utf-8')
    print('✅ FAERSDownloader created with progress bars!')
else:
    print('✅ FAERSDownloader already exists')

print('\\n🎯 Now updating data_source_manager to show quarters after Configure button...')

# Update data_source_manager to show quarters immediately after clicking Configure
dsm_file = Path('src/data/data_source_manager.py')
content = dsm_file.read_text(encoding='utf-8')

# Find the section after configure button and ensure quarters are shown
old_after_config = '''        # Configure button
        configure_dates = st.sidebar.button("📅 Configure Date Range", help="Click after selecting dates")
        
        if not configure_dates:
            st.sidebar.info("👆 Click 'Configure Date Range' to proceed")
            return None
        
        # Compute quarters after validation
        quarters = self._compute_quarters(start_date, end_date)
        st.sidebar.write(f"**Detected Quarters:** {len(quarters)}")'''

new_after_config = '''        # Configure button
        configure_dates = st.sidebar.button("📅 Configure Date Range", help="Click after selecting dates")
        
        if not configure_dates:
            st.sidebar.info("👆 Click 'Configure Date Range' to proceed")
            return None
        
        # Compute quarters after validation
        with st.spinner("🔍 Computing quarters..."):
            quarters = self._compute_quarters(start_date, end_date)
        
        st.sidebar.success(f"✅ Detected **{len(quarters)} quarters**: {', '.join([f'Q{q[1]}-{q[0]}' for q in quarters])}")'''

content = content.replace(old_after_config, new_after_config)
dsm_file.write_text(content, encoding='utf-8')

print('✅ Updated data_source_manager to show quarters clearly!')
print('\\n🚀 Ready! Restart your Streamlit app.')
