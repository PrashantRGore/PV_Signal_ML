from pathlib import Path

print('Fixing FAERSDownloader to actually USE cache...')
faers_file = Path('src/data/faers_downloader.py')

# Complete working version with cache
new_content = '''"""
FAERS Data Downloader with Smart Caching
"""

import requests
import pandas as pd
from typing import List, Tuple
import streamlit as st
import zipfile
import io
from src.utils.logger import setup_logger
from src.data.faers_cache import FAERSCache

logger = setup_logger(__name__)

class FAERSDownloader:
    BASE_URL = "https://fis.fda.gov/content/Exports/faers_ascii_{year}Q{quarter}.zip"
    
    def __init__(self):
        self.cache = FAERSCache()
        logger.info(f"FAERSDownloader initialized with cache at {self.cache.cache_dir}")
    
    def download_and_process(self, quarters: List[Tuple[int, int]], 
                            start_date, end_date) -> pd.DataFrame:
        """Download with smart caching"""
        logger.info(f"Processing {len(quarters)} quarters")
        all_data = []
        progress_bar = st.sidebar.progress(0)
        status = st.sidebar.empty()
        
        for idx, (year, quarter) in enumerate(quarters):
            status.write(f"📊 Quarter {idx+1}/{len(quarters)}: Q{quarter}-{year}")
            
            try:
                # CHECK CACHE FIRST
                if self.cache.is_cached(year, quarter):
                    info = self.cache.get_info(year, quarter)
                    st.sidebar.info(f"💾 Loading Q{quarter}-{year} from cache...")
                    data = self.cache.load(year, quarter)
                    st.sidebar.success(f"✅ Loaded from cache: {len(data):,} records ({info.get('size_mb', 0):.1f} MB)")
                else:
                    # Download from FDA
                    st.sidebar.info(f"📥 Downloading Q{quarter}-{year} from FDA...")
                    data = self._download_quarter(year, quarter, status)
                    
                    if data is not None and len(data) > 0:
                        # SAVE TO CACHE
                        cache_path = self.cache.save(year, quarter, data)
                        st.sidebar.success(f"💾 Saved to cache: {cache_path.name}")
                        logger.info(f"Cached Q{quarter}-{year} at {cache_path}")
                
                if data is not None and len(data) > 0:
                    all_data.append(data)
            
            except Exception as e:
                st.sidebar.error(f"❌ Q{quarter}-{year}: {str(e)}")
                logger.error(f"Error Q{quarter}-{year}: {e}", exc_info=True)
            
            progress_bar.progress((idx + 1) / len(quarters))
        
        status.empty()
        progress_bar.empty()
        
        if not all_data:
            raise Exception("No data loaded")
        
        combined = pd.concat(all_data, ignore_index=True)
        st.sidebar.success(f"🎉 Total: {len(combined):,} drug-event pairs")
        logger.info(f"Combined {len(combined)} records from {len(all_data)} quarters")
        return combined
    
    def _download_quarter(self, year, quarter, status):
        """Download and parse a quarter from FDA"""
        url = self.BASE_URL.format(year=year, quarter=quarter)
        status.write(f"🌐 Downloading: {url}")
        
        response = requests.get(url, timeout=120)
        if response.status_code != 200:
            raise Exception(f"HTTP {response.status_code}")
        
        status.write("📦 Extracting ZIP file...")
        with zipfile.ZipFile(io.BytesIO(response.content)) as z:
            drug_files = [f for f in z.namelist() if 'DRUG' in f.upper() and f.endswith('.txt')]
            reac_files = [f for f in z.namelist() if 'REAC' in f.upper() and f.endswith('.txt')]
            
            if not drug_files or not reac_files:
                raise Exception("Required files not found in ZIP")
            
            status.write("📊 Parsing DRUG file...")
            drug_df = pd.read_csv(z.open(drug_files[0]), sep='$', encoding='latin1', 
                                 low_memory=False, on_bad_lines='skip')
            
            status.write("📊 Parsing REAC file...")
            reac_df = pd.read_csv(z.open(reac_files[0]), sep='$', encoding='latin1',
                                 low_memory=False, on_bad_lines='skip')
        
        status.write("🔗 Merging drug-event pairs...")
        id_col_drug = 'primaryid' if 'primaryid' in drug_df.columns else 'caseid'
        id_col_reac = 'primaryid' if 'primaryid' in reac_df.columns else 'caseid'
        
        merged = drug_df.merge(reac_df, left_on=id_col_drug, right_on=id_col_reac, how='inner')
        
        result = pd.DataFrame({
            'case_id': merged[id_col_drug],
            'drug_name': merged.get('drugname', merged.get('drug_name', '')),
            'event_term': merged.get('pt', merged.get('event_term', '')),
            'year': year,
            'quarter': quarter
        })
        
        # Clean data
        result = result.dropna(subset=['drug_name', 'event_term'])
        result = result[result['drug_name'].str.strip() != '']
        result = result[result['event_term'].str.strip() != '']
        result['drug_name'] = result['drug_name'].str.upper().str.strip()
        result['event_term'] = result['event_term'].str.upper().str.strip()
        
        return result
'''

faers_file.write_text(new_content, encoding='utf-8')
print('✅ COMPLETE: FAERSDownloader now uses cache!')
print('✅ Cache checks: self.cache.is_cached()')
print('✅ Cache loads: self.cache.load()')
print('✅ Cache saves: self.cache.save()')
print()
print('🚀 Restart your app now!')
