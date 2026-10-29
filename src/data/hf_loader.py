"""
HuggingFace Dataset Loader - Complete Implementation
- Original HuggingFaceLoader class (7 methods) - PRESERVED
- Module-level functions for __init__.py exports - ADDED
- DuckDB support for efficient querying - ADDED
- v3 dataset with 70K+ signals - UPDATED
- Canonical column mapping - ADDED
"""

import os
import pandas as pd
import duckdb
from datasets import load_dataset
from typing import Optional, Dict
from pathlib import Path
import logging

import config
from src.data.canonical_schema import canonicalize_signals, CanonicalSchema
from src.utils.logger import setup_logger

logger = setup_logger(__name__)

# =============================================================================
# CONFIGURATION
# =============================================================================

DEFAULT_REPO_ID = "PrashantRGore/synthetic-faers-1m-v3"
DEFAULT_CSV_FILE = "synthetic_faers_1m_v3.csv"
# Dynamic cache directory with backward compatibility
_env_cache = os.environ.get("PV_SIGNAL_CACHE_DIR")
if _env_cache:
    CACHE_DIR = Path(_env_cache)
elif Path("C:/Users/koreo/PV_Signal_Cache").exists():
    CACHE_DIR = Path("C:/Users/koreo/PV_Signal_Cache")
else:
    CACHE_DIR = Path(getattr(config, "HF_CACHE_DIR", config.DATA_DIR / "hf_cache"))

DUCKDB_PATH = CACHE_DIR / "faers_demo.duckdb" 

# Canonical column mapping (v3 -> internal)
COLUMN_MAPPING = {
    'case_id': 'case_id',
    'receive_date': 'fda_dt',
    'country': 'country',
    'age': 'age',
    'age_group': 'age_group',
    'sex': 'sex',
    'weight_kg': 'weight',
    'suspect_drug': 'drug_name',
    'indication': 'drug_indication',
    'route': 'route_of_admin',
    'dose': 'dose',
    'dose_unit': 'dose_unit',
    'dose_frequency': 'dose_freq',
    'treatment_duration_days': 'therapy_duration',
    'event_llt': 'llt',
    'event_pt': 'event_term',  # Maps to canonical event_term
    'event_hlt': 'hlt',
    'event_hlgt': 'hlgt',
    'event_soc': 'soc',
    'time_to_onset_days': 'ttp_days',
    'event_duration_days': 'event_duration',
    'seriousness': 'serious',
    'outcome': 'outcome',
    'action_taken': 'action_drug',
    'causality_assessment': 'causality',
    'dechallenge': 'dechallenge',
    'rechallenge': 'rechallenge',
    'concomitant_medications': 'concomitant_drugs',
    'medical_history': 'medical_history',
    'reporter_type': 'reporter_qualification',
    'report_type': 'report_source',
    'alt_u_l': 'alt',
    'ast_u_l': 'ast',
    'bilirubin_mg_dl': 'bilirubin',
    'creatinine_mg_dl': 'creatinine',
    'bun_mg_dl': 'bun',
}


# =============================================================================
# ORIGINAL CLASS - PRESERVED FOR BACKWARDS COMPATIBILITY
# =============================================================================

class HuggingFaceLoader:
    """
    Original class-based loader (v2 API).
    Load synthetic 1M ICSR demo dataset from HuggingFace.
    PRESERVED for backwards compatibility.
    """

    def __init__(self, repo_id: Optional[str] = None):
        self.repo_id = repo_id or DEFAULT_REPO_ID
        self.cache_dir = config.HF_CACHE_DIR
        self.data_file = DEFAULT_CSV_FILE

    def load_dataset(self) -> pd.DataFrame:
        """
        Load Demo 1M dataset from HuggingFace, compute signals, and canonicalize.
        """
        logger.info(f"Loading raw ICSR data from HuggingFace: {self.repo_id}")
        logger.info(f"Using specific file: {self.data_file}")

        try:
            # Try local cache first
            cache_path = CACHE_DIR / self.data_file
            if cache_path.exists():
                logger.info(f"Loading from local cache: {cache_path}")
                df = pd.read_csv(cache_path)
            else:
                # Load from HuggingFace
                dataset = load_dataset(
                    self.repo_id,
                    data_files=self.data_file,
                    cache_dir=self.cache_dir,
                    split="train",
                )
                df = dataset.to_pandas()
            
            logger.info(f"Loaded {len(df):,} raw ICSR records")

            # Check if data is raw ICSR or pre-computed signals
            has_prr = "prr" in df.columns or "PRR" in df.columns

            if has_prr:
                # Pre-computed signals - just canonicalize
                canonical_df = canonicalize_signals(df, source="demo_hf", preserve_tier2=True)
            else:
                # Raw ICSR data - compute signals first
                icsr_df = self._map_to_icsr_format(df)
                from src.stats_engine.disproportionality import DisproportionalityAnalysis
                stats_engine = DisproportionalityAnalysis(min_cases=config.MIN_CASE_COUNT)
                signals_df = stats_engine.compute_signals(icsr_df)
                canonical_df = signals_df

            CanonicalSchema.validate(canonical_df)
            logger.info(f"✅ Demo 1M dataset ready: {len(canonical_df):,} signals")
            
            return canonical_df

        except Exception as e:
            logger.error(f"Failed to load dataset: {e}")
            raise

    def _map_to_icsr_format(self, df: pd.DataFrame) -> pd.DataFrame:
        """Map v3 columns to ICSR format."""
        return df.rename(columns=COLUMN_MAPPING)

    def validate_schema(self, df: pd.DataFrame) -> bool:
        """Validate schema has required columns."""
        has_case = 'case_id' in df.columns
        has_drug = any(c in df.columns for c in ['drug_name', 'suspect_drug'])
        has_event = any(c in df.columns for c in ['event_term', 'pt', 'event_pt'])
        return has_case and has_drug and has_event

    def get_sample(self, n: int = 10000) -> pd.DataFrame:
        """Get sample of n records."""
        df = self.load_dataset()
        return df.sample(n=min(n, len(df)), random_state=42)

    def check_availability(self) -> bool:
        """Check if dataset is available."""
        try:
            cache_path = CACHE_DIR / self.data_file
            return cache_path.exists()
        except:
            return False

    def get_metadata(self) -> dict:
        """Get dataset metadata."""
        return {
            'repo_id': self.repo_id,
            'csv_file': self.data_file,
            'version': 'v3',
            'injected_signals': '70,803 (7.1%)'
        }


# =============================================================================
# MODULE-LEVEL FUNCTIONS - REQUIRED BY __init__.py
# =============================================================================

def load_demo_dataset_duckdb(use_cache: bool = True) -> pd.DataFrame:
    """
    Load demo dataset using DuckDB for efficient querying.
    This is the main function used by data_source_manager.
    
    Args:
        use_cache: Use local cache if available
        
    Returns:
        DataFrame with canonical column names
    """
    try:
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        csv_path = CACHE_DIR / DEFAULT_CSV_FILE
        
        if not csv_path.exists() or not use_cache:
            # Download from HuggingFace
            logger.info(f"📥 Downloading from HuggingFace: {DEFAULT_REPO_ID}")
            print(f"📥 First-time download: {DEFAULT_REPO_ID}")
            print(f"📄 Loading: {DEFAULT_CSV_FILE} (50+ Tier 2 features)")
            
            dataset = load_dataset(
                DEFAULT_REPO_ID,
                data_files=DEFAULT_CSV_FILE,
                split="train"
            )
            df = dataset.to_pandas()
            df.to_csv(csv_path, index=False)
            logger.info(f"Cached to: {csv_path}")
        else:
            # Load from cache
            logger.info(f"📁 Loading from local cache: {csv_path}")
            print(f"📁 Loading from local cache")
            df = pd.read_csv(csv_path)
        
        # Apply canonical mapping
        df_canonical = df.rename(columns=COLUMN_MAPPING)
        
        # Create DuckDB database
        conn = duckdb.connect(str(DUCKDB_PATH))
        conn.execute("DROP TABLE IF EXISTS faers_demo")
        conn.execute("CREATE TABLE faers_demo AS SELECT * FROM df_canonical")
        
        logger.info(f"✅ Loaded {len(df_canonical):,} records into DuckDB")
        print(f"✅ Loaded {len(df_canonical):,} records with {len(df_canonical.columns)} features")
        
        conn.close()
        return df_canonical
        
    except Exception as e:
        logger.error(f"Failed to load demo dataset: {e}")
        raise


def get_cache_info() -> Dict:
    """
    Get information about cached datasets.
    Used by data_source_manager for cache management.
    """
    try:
        csv_path = CACHE_DIR / DEFAULT_CSV_FILE
        duckdb_path = DUCKDB_PATH
        
        info = {
            'cache_dir': str(CACHE_DIR),
            'csv_cached': csv_path.exists(),
            'duckdb_cached': duckdb_path.exists(),
            'repo_id': DEFAULT_REPO_ID,
            'version': 'v3'
        }
        
        if csv_path.exists():
            size_mb = csv_path.stat().st_size / (1024 * 1024)
            info['csv_size_mb'] = round(size_mb, 1)
        
        if duckdb_path.exists():
            size_mb = duckdb_path.stat().st_size / (1024 * 1024)
            info['duckdb_size_mb'] = round(size_mb, 1)
        
        return info
        
    except Exception as e:
        logger.error(f"Failed to get cache info: {e}")
        return {'error': str(e)}


def query_signals_duckdb(sql_query: str) -> pd.DataFrame:
    """
    Execute arbitrary SQL query on DuckDB faers_demo table.
    Used by data_source_manager for flexible querying.
    """
    try:
        if not DUCKDB_PATH.exists():
            logger.warning('DuckDB database not found. Run load_demo_dataset_duckdb() first.')
            return pd.DataFrame()
        
        conn = duckdb.connect(str(DUCKDB_PATH))
        result = conn.execute(sql_query).df()
        conn.close()
        
        logger.info(f'Query returned {len(result)} rows')
        return result
        
    except Exception as e:
        logger.error(f'Query failed: {e}')
        raise


def clear_cache() -> Dict:
    """
    Clear all cached datasets (CSV and DuckDB).
    Used for forcing fresh download or testing.
    """
    try:
        result = {
            'csv_deleted': False,
            'duckdb_deleted': False,
            'cache_dir': str(CACHE_DIR)
        }
        
        csv_path = CACHE_DIR / DEFAULT_CSV_FILE
        if csv_path.exists():
            csv_path.unlink()
            result['csv_deleted'] = True
            logger.info(f'Deleted CSV cache: {csv_path}')
        
        if DUCKDB_PATH.exists():
            DUCKDB_PATH.unlink()
            result['duckdb_deleted'] = True
            logger.info(f'Deleted DuckDB cache: {DUCKDB_PATH}')
        
        return result
        
    except Exception as e:
        logger.error(f'Failed to clear cache: {e}')
        return {'error': str(e)}
