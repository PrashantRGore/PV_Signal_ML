from pathlib import Path

faers_file = Path('src/data/faers_downloader.py')

faers_code = '''"""
FAERS Data Downloader: Downloads REAL quarterly FAERS data from FDA
"""

import requests
import pandas as pd
from typing import List, Tuple
import streamlit as st
import zipfile
import io
from src.utils.logger import setup_logger

logger = setup_logger(__name__)

class FAERSDownloader:
    """Downloads and processes FAERS quarterly data from FDA"""
    
    # Real FDA quarterly data download URLs
    BASE_URL = "https://fis.fda.gov/content/Exports/faers_ascii_{year}Q{quarter}.zip"
    
    def __init__(self):
        self.data = None
    
    def download_and_process(self, quarters: List[Tuple[int, int]], 
                            start_date, end_date) -> pd.DataFrame:
        """
        Download and process REAL FAERS data for specified quarters
        
        Args:
            quarters: List of (year, quarter) tuples
            start_date: Start date for filtering
            end_date: End date for filtering
        
        Returns:
            Processed DataFrame with ICSR records
        """
        logger.info(f"Starting REAL FAERS download for {len(quarters)} quarters")
        st.sidebar.info(f"📥 Downloading {len(quarters)} quarters from FDA...")
        
        all_data = []
        
        # Create progress bar in sidebar
        progress_bar = st.sidebar.progress(0, text="Initializing download...")
        
        for idx, (year, quarter) in enumerate(quarters):
            # Update progress bar
            progress_pct = idx / len(quarters)
            progress_bar.progress(progress_pct, text=f"📥 Downloading Q{quarter}-{year}... ({idx+1}/{len(quarters)})")
            
            st.sidebar.write(f"🔄 Processing **Q{quarter}-{year}**...")
            
            try:
                # Download real FAERS data
                quarter_data = self._fetch_quarter_data_real(year, quarter)
                
                if quarter_data is not None and len(quarter_data) > 0:
                    all_data.append(quarter_data)
                    st.sidebar.success(f"✅ Q{quarter}-{year}: {len(quarter_data):,} records")
                else:
                    st.sidebar.warning(f"⚠️ Q{quarter}-{year}: No data found")
            
            except Exception as e:
                st.sidebar.error(f"❌ Q{quarter}-{year}: {str(e)}")
                logger.error(f"Failed to download Q{quarter}-{year}: {e}")
        
        # Final progress
        progress_bar.progress(1.0, text="✅ Download complete!")
        progress_bar.empty()
        
        if len(all_data) == 0:
            raise Exception("No data downloaded for selected quarters")
        
        # Combine all quarters
        combined_data = pd.concat(all_data, ignore_index=True)
        logger.info(f"Downloaded {len(combined_data)} REAL records from FDA")
        st.sidebar.success(f"🎉 Total: **{len(combined_data):,}** records from FDA FAERS")
        
        return combined_data
    
    def _fetch_quarter_data_real(self, year: int, quarter: int) -> pd.DataFrame:
        """
        Fetch REAL data for a specific quarter from FDA
        Downloads ZIP file, extracts, and parses DRUG and REAC files
        """
        # Construct download URL
        url = self.BASE_URL.format(year=year, quarter=quarter)
        
        st.sidebar.write(f"🌐 Downloading from: {url}")
        
        # Download ZIP file
        response = requests.get(url, timeout=60)
        
        if response.status_code != 200:
            raise Exception(f"Failed to download: HTTP {response.status_code}")
        
        # Extract ZIP in memory
        with zipfile.ZipFile(io.BytesIO(response.content)) as z:
            # Find DRUG and REAC files
            drug_file = [f for f in z.namelist() if 'DRUG' in f.upper() and f.endswith('.txt')][0]
            reac_file = [f for f in z.namelist() if 'REAC' in f.upper() and f.endswith('.txt')][0]
            
            # Read files
            drug_df = pd.read_csv(z.open(drug_file), sep='$', encoding='latin1', 
                                 low_memory=False, on_bad_lines='skip')
            reac_df = pd.read_csv(z.open(reac_file), sep='$', encoding='latin1',
                                 low_memory=False, on_bad_lines='skip')
        
        # Merge on case ID
        # FAERS uses primaryid or caseid depending on year
        id_col_drug = 'primaryid' if 'primaryid' in drug_df.columns else 'caseid'
        id_col_reac = 'primaryid' if 'primaryid' in reac_df.columns else 'caseid'
        
        merged = drug_df.merge(reac_df, 
                              left_on=id_col_drug, 
                              right_on=id_col_reac, 
                              how='inner')
        
        # Standardize columns for PV Signal Detection
        result = pd.DataFrame({
            'case_id': merged[id_col_drug],
            'drug_name': merged['drugname'] if 'drugname' in merged.columns else merged.get('drug_name', ''),
            'event_term': merged['pt'] if 'pt' in merged.columns else merged.get('event_term', ''),
            'year': year,
            'quarter': quarter
        })
        
        return result.dropna(subset=['drug_name', 'event_term'])
'''

faers_file.write_text(faers_code, encoding='utf-8')
print('✅ REAL FAERS downloader implemented!')
print('📥 Will download actual quarterly ZIP files from FDA')
print('🔗 URL format: https://fis.fda.gov/content/Exports/faers_ascii_YYYYQQ.zip')
