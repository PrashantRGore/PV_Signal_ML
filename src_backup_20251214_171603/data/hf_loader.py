"""
Fixed HuggingFace Loader with correct column mapping
"""
import pandas as pd
from pathlib import Path
from huggingface_hub import hf_hub_download
import config
from src.utils.logger import setup_logger

logger = setup_logger(__name__)


class HuggingFaceLoader:
    def __init__(self):
        self.repo_id = "PrashantRGore/pv-signal-ml-data"
        self.cache_dir = config.HF_CACHE_DIR
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        
    def load_dataset(self) -> pd.DataFrame:
        logger.info(f"Loading demo dataset from {self.repo_id}")
        
        # Check cache
        cached_file = self.cache_dir / "demo_combined.parquet"
        if cached_file.exists():
            logger.info("Loading from local cache")
            data = pd.read_parquet(cached_file)
            # Ensure correct column names
            if 'pt' in data.columns and 'event_term' not in data.columns:
                data = data.rename(columns={'pt': 'event_term'})
            return data
        
        try:
            logger.info("Downloading from HuggingFace...")
            
            # Try standardized files first
            try:
                drug_path = hf_hub_download(repo_id=self.repo_id, filename="drug.parquet", 
                                           repo_type="dataset", cache_dir=str(self.cache_dir))
                reac_path = hf_hub_download(repo_id=self.repo_id, filename="reac.parquet",
                                           repo_type="dataset", cache_dir=str(self.cache_dir))
                demo_path = hf_hub_download(repo_id=self.repo_id, filename="demo.parquet",
                                           repo_type="dataset", cache_dir=str(self.cache_dir))
                outc_path = hf_hub_download(repo_id=self.repo_id, filename="outc.parquet",
                                           repo_type="dataset", cache_dir=str(self.cache_dir))
                
                drug_df = pd.read_parquet(drug_path)
                reac_df = pd.read_parquet(reac_path)
                demo_df = pd.read_parquet(demo_path)
                outc_df = pd.read_parquet(outc_path)
                
                logger.info("✅ Loaded standardized files from root")
                
            except Exception as e:
                logger.warning(f"Standardized files not found, trying original location: {e}")
                
                # Fallback to original location
                drug_path = hf_hub_download(repo_id=self.repo_id, 
                                           filename="data/clean\\ICSR_PRODUCT_1M.parquet",
                                           repo_type="dataset", cache_dir=str(self.cache_dir))
                reac_path = hf_hub_download(repo_id=self.repo_id,
                                           filename="data/clean\\ICSR_EVENT_1M.parquet",
                                           repo_type="dataset", cache_dir=str(self.cache_dir))
                demo_path = hf_hub_download(repo_id=self.repo_id,
                                           filename="data/clean\\ICSR_PATIENT_1M.parquet",
                                           repo_type="dataset", cache_dir=str(self.cache_dir))
                
                drug_df = pd.read_parquet(drug_path)
                reac_df = pd.read_parquet(reac_path)
                demo_df = pd.read_parquet(demo_path)
                outc_df = reac_df.copy()  # Use event data for outcomes
                
                logger.info("✅ Loaded original files from data/clean/")
            
            # Standardize column names
            # Find case_id column
            drug_case_col = next((c for c in drug_df.columns if 'case' in c.lower() and 'id' in c.lower()), 'CASE_ID')
            reac_case_col = next((c for c in reac_df.columns if 'case' in c.lower() and 'id' in c.lower()), 'CASE_ID')
            demo_case_col = next((c for c in demo_df.columns if 'case' in c.lower() and 'id' in c.lower()), 'CASE_ID')
            
            # Find drug column
            drug_col = next((c for c in drug_df.columns if any(x in c.lower() for x in ['product', 'drug'])), None)
            
            # Find event column
            event_col = next((c for c in reac_df.columns if any(x in c.lower() for x in ['event', 'pt', 'reaction'])), None)
            
            logger.info(f"Detected columns - Drug: {drug_col}, Event: {event_col}")
            
            # Rename to standard names
            drug_df = drug_df.rename(columns={drug_case_col: 'case_id', drug_col: 'drug_name'})
            reac_df = reac_df.rename(columns={reac_case_col: 'case_id', event_col: 'event_term'})
            demo_df = demo_df.rename(columns={demo_case_col: 'case_id'})
            
            # Find and rename age/sex columns in demo
            age_col = next((c for c in demo_df.columns if 'age' in c.lower()), None)
            sex_col = next((c for c in demo_df.columns if any(x in c.lower() for x in ['sex', 'gender'])), None)
            
            if age_col and age_col != 'age':
                demo_df = demo_df.rename(columns={age_col: 'age'})
            if sex_col and sex_col != 'sex':
                demo_df = demo_df.rename(columns={sex_col: 'sex'})
            
            # Merge tables
            combined = drug_df[['case_id', 'drug_name']].merge(
                reac_df[['case_id', 'event_term']],
                on='case_id',
                how='inner'
            )
            
            # Add demographics if available
            demo_cols = ['case_id']
            if 'age' in demo_df.columns:
                demo_cols.append('age')
            if 'sex' in demo_df.columns:
                demo_cols.append('sex')
            
            combined = combined.merge(demo_df[demo_cols], on='case_id', how='left')
            
            # Cache
            combined.to_parquet(cached_file, index=False)
            logger.info(f"✅ Loaded and merged {len(combined):,} records")
            logger.info(f"Columns: {list(combined.columns)}")
            
            return combined
            
        except Exception as e:
            logger.error(f"Failed to load from HuggingFace: {e}")
            
            # Fallback to local
            local_dir = Path(r'C:\Users\koreo\PV_Signal_ML\data\huggingface_upload')
            if local_dir.exists():
                drug_df = pd.read_parquet(local_dir / 'drug.parquet')
                reac_df = pd.read_parquet(local_dir / 'reac.parquet')
                demo_df = pd.read_parquet(local_dir / 'demo.parquet')
                
                # Rename pt to event_term
                if 'pt' in reac_df.columns:
                    reac_df = reac_df.rename(columns={'pt': 'event_term'})
                
                combined = drug_df.merge(reac_df, on='case_id', how='inner')
                combined = combined.merge(demo_df, on='case_id', how='left')
                
                logger.info(f"✅ Loaded {len(combined):,} from local files")
                return combined
            
            raise ValueError(f"Could not load dataset: {e}")
