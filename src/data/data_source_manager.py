"""
Data Source Manager: Handles Local, Demo (HuggingFace), and FAERS Live sources
NOW WITH: Canonical schema integration - all sources output canonical format
"""

import streamlit as st
import pandas as pd
from pathlib import Path
from datetime import datetime
from typing import Optional, Tuple, List
import json

import config
from src.data.hf_loader import load_demo_dataset_duckdb, get_cache_info
from src.data.faers_downloader import FAERSDownloader
from src.data.canonical_schema import canonicalize_signals, CanonicalSchema
from src.utils.logger import setup_logger

logger = setup_logger(__name__)

# ===== Export canonical functions for backward compatibility =====

def normalize_signals(df: pd.DataFrame, source: str) -> pd.DataFrame:
    """
    Legacy wrapper - redirects to canonical implementation.
    
    DEPRECATED: Use canonicalize_signals() directly.
    """
    logger.warning("normalize_signals() is deprecated. Use canonicalize_signals() instead.")
    return canonicalize_signals(df, source, preserve_tier2=True)


# ===== Export feature definitions from canonical schema =====

def get_tier1_features() -> List[str]:
    """Get Tier 1 feature columns for ML"""
    return CanonicalSchema.get_tier1_feature_columns()

def get_tier2_features(df: pd.DataFrame) -> List[str]:
    """Get Tier 2 feature columns for ML"""
    return CanonicalSchema.get_tier2_columns(df)

# Legacy constants for backward compatibility
TIER1_FEATURES = CanonicalSchema.get_tier1_feature_columns()
TIER2_FEATURES_1M = []  # Dynamically populated at runtime


class DataSourceManager:
    """Manages three mutually exclusive data sources for PV signal detection"""

    def __init__(self):
        self.current_source = None
        self.data = None
        self.metadata = {}

    def select_data_source(self) -> str:
        """
        Streamlit UI for data source selection
        Returns: Selected source type ('local', 'demo_hf', 'faers_live')
        """
        st.sidebar.header("📊 Data Source Selection")
        source_type = st.sidebar.radio(
            "Choose Data Source:",
            options=config.DATA_SOURCE_OPTIONS,
            format_func=lambda x: {
                "local": "🗂️ Local ICSR Dataset",
                "demo_hf": "🤗 Demo 1M Dataset (HuggingFace)",
                "faers_live": "📡 Live FAERS Quarterly Data",
            }[x],
            help="Select ONE data source. Sources are mutually exclusive.",
        )
        self.current_source = source_type
        return source_type

    def load_local_dataset(self) -> Optional[pd.DataFrame]:
        """Load user-provided local dataset"""
        st.sidebar.subheader("Local Dataset Configuration")
        
        input_method = st.sidebar.radio(
            "Input Method:",
            ["File Path", "File Upload"],
        )

        try:
            if input_method == "File Path":
                file_path = st.sidebar.text_input(
                    "Enter file path (CSV/Parquet/SQLite):",
                    placeholder=r"C:\path\to\data.csv",
                )
                if file_path and Path(file_path).exists():
                    data = self._read_file(file_path)
                    self._validate_icsr_schema(data)
                    self.metadata = {
                        "source": "local",
                        "file_path": file_path,
                        "load_time": datetime.now().isoformat(),
                        "record_count": len(data),
                    }
                    logger.info(f"Loaded local dataset: {len(data)} records")
                    return data
                elif file_path:
                    st.sidebar.error(f"❌ File not found: {file_path}")

            else:  # File Upload
                uploaded_file = st.sidebar.file_uploader(
                    "Upload ICSR file:",
                    type=["csv", "parquet", "db"],
                    help="Supported: CSV, Parquet, SQLite",
                )
                if uploaded_file:
                    data = self._read_uploaded_file(uploaded_file)
                    self._validate_icsr_schema(data)
                    self.metadata = {
                        "source": "local_upload",
                        "filename": uploaded_file.name,
                        "load_time": datetime.now().isoformat(),
                        "record_count": len(data),
                    }
                    logger.info(f"Uploaded dataset: {len(data)} records")
                    return data

        except Exception as e:
            st.sidebar.error(f"❌ Error loading local dataset: {str(e)}")
            logger.error(f"Local dataset load failed: {e}")
        
        return None

    def load_demo_dataset(self) -> Optional[pd.DataFrame]:
        """Load 1M synthetic demo dataset from HuggingFace"""
        st.sidebar.subheader("Demo Dataset (HuggingFace)")
        st.sidebar.info(
            """
📦 **1M Synthetic ICSR Records**
- NO REAL PATIENT DATA
- GDPR/HIPAA Compliant
- Educational Use Only
            """
        )

        use_demo = st.sidebar.checkbox(
            "Load Demo Dataset",
            value=False,
            help="1M synthetic records from HuggingFace",
        )

        if use_demo:
            try:
                with st.spinner("📥 Downloading demo dataset from HuggingFace..."):
                    data = load_demo_dataset_duckdb()
                    
                    self.metadata = {
                        "source": "demo_hf",
                        "repo_id": config.HF_REPO_ID,
                        "load_time": datetime.now().isoformat(),
                        "record_count": len(data),
                        "synthetic": True,
                        "legal_status": "GDPR_HIPAA_COMPLIANT",
                    }
                    
                    st.sidebar.success(f"✅ Loaded {len(data):,} demo records")
                    logger.info(f"Demo dataset loaded: {len(data)} records")
                    return data

            except Exception as e:
                st.sidebar.error(f"❌ Failed to load demo dataset: {str(e)}")
                logger.error(f"Demo dataset load failed: {e}")
        
        return None

    def load_faers_dataset(
        self, 
        start_date: datetime = None,  # ✅ ADD PARAMETERS
        end_date: datetime = None      # ✅ ADD PARAMETERS
    ) -> Optional[pd.DataFrame]:
        """
        Load FAERS data via precomputed ranked signals for selected date range.
    
        ✅ FIXED: Now accepts start_date and end_date parameters
        """
        st.sidebar.subheader("Live FAERS Quarterly Data")

        # ✅ Use provided dates if available
        if start_date is None:
            start_date = st.session_state.get('faers_start', datetime(2024, 1, 1))
        if end_date is None:
            end_date = st.session_state.get('faers_end', datetime.now())

        # Store in session state for later use
        if "faers_configured" not in st.session_state:
            st.session_state.faers_configured = False
            st.session_state.faers_quarters = []
            st.session_state.faers_start = start_date
            st.session_state.faers_end = end_date

        # Get today's date
        today = datetime.now()

        # Date inputs (only show if not already provided)
        if start_date == st.session_state.get('faers_start') and end_date == st.session_state.get('faers_end'):
            col1, col2 = st.sidebar.columns(2)
            with col1:
                start_date = st.date_input(
                    "Start Date",
                    value=start_date,
                    min_value=datetime(2004, 1, 1),
                    max_value=today,
                    help="FAERS data available from 2004 onwards",
                    key="faers_start_input"
                )
            with col2:
                end_date = st.date_input(
                    "End Date",
                    value=end_date,
                    min_value=datetime(2004, 1, 1),
                    max_value=today,
                    help="Select end date",
                    key="faers_end_input"
                )

        # Validate dates
        if start_date > end_date:
            st.sidebar.error("⚠️ End date must be >= Start date")
            return None

        # Configure button - SAVE TO SESSION STATE
        if st.sidebar.button("📅 Configure Date Range", key="faers_config_btn"):
            st.session_state.faers_configured = True
            st.session_state.faers_start = start_date
            st.session_state.faers_end = end_date
            st.session_state.faers_quarters = self._compute_quarters(
                start_date, end_date
            )
            st.rerun()

        # Check if configured
        if not st.session_state.faers_configured:
            st.sidebar.info("👆 Click 'Configure Date Range' to proceed")
            return None

        # Show quarters (from session state)
        quarters = st.session_state.faers_quarters
        st.sidebar.success(
            f"✅ Detected {len(quarters)} quarters: "
            f"{', '.join([f'Q{q[1]}-{q[0]}' for q in quarters])}"
        )

        # Quarter selection
        selected_quarters = st.sidebar.multiselect(
            "Select quarters to include:",
            options=quarters,
            default=quarters,
            format_func=lambda q: f"Q{q[1]}-{q[0]}",
            key="faers_quarter_select"
        )

        if not selected_quarters:
            st.sidebar.error("⚠️ Select at least one quarter")
            return None

        # Convert to labels used by FAERSSignalBuilder (e.g. "2025Q2")
        quarter_labels = [f"{y}Q{q}" for (y, q) in selected_quarters]

        # ✅ Check cache first (avoid re-downloading)
        raw_dir = Path("data/faers_raw")
        raw_dir.mkdir(parents=True, exist_ok=True)
    
        missing = []
        existing = []
        for q_label in quarter_labels:
            zip_path = raw_dir / f"{q_label}.zip"
            if zip_path.exists():
                existing.append(q_label)
            else:
                missing.append(q_label)

        # Show status
        if existing:
            st.sidebar.success(f"✅ Cached: {', '.join(existing)}")
        if missing:
            st.sidebar.warning(f"⚠️ Need to download: {', '.join(missing)}")

        # BUTTON 1: Download missing quarters (only show if any missing)
        if missing:
            if st.sidebar.button(
                f"📥 Download {len(missing)} Missing Quarter(s)", 
                type="secondary",
                key="faers_download_btn"
            ):
                try:
                    downloader = FAERSDownloader()
                    with st.spinner(f"Downloading {len(missing)} quarter(s)..."):
                        downloader.download_raw_zips(missing)
                    st.sidebar.success("✅ Download complete! Now click Build/Load below.")
                    logger.info(f"Downloaded raw zips: {missing}")
                    st.rerun()
                except Exception as e:
                    st.sidebar.error(f"❌ Download failed: {str(e)}")
                    logger.error("Raw zip download error", exc_info=True)
                    return None

            st.sidebar.info("👆 Download missing quarters first before building signals")
            return None

        # BUTTON 2: Build / Load FAERS Signals (only show when all quarters available)
        if st.sidebar.button("⚙️ Build / Load FAERS Signals", type="primary", key="faers_build_btn"):
            try:
                start_q = quarter_labels[0]
                end_q = quarter_labels[-1]
                key = f"{start_q}_{end_q}"

                # ✅ Check if already built (cache)
                output_dir = Path("data/outputs")
                output_dir.mkdir(parents=True, exist_ok=True)
                signals_path = output_dir / f"ranked_signals_{key}_predictions.csv"
            
                if signals_path.exists():
                    st.sidebar.info(f"✅ Loading cached signals: {signals_path.name}")
                    data = pd.read_csv(signals_path)
                else:
                    # Build signals via FAERSSignalBuilder
                    config_path = Path("data/faers_config.json")
                    config_payload = {"quarters": quarter_labels}
                    config_path.parent.mkdir(parents=True, exist_ok=True)
                    with open(config_path, "w") as f:
                        json.dump(config_payload, f, indent=2)

                    st.sidebar.info("⚙️ Building ranked signals for selected quarters...")
                    with st.spinner("Computing FAERS signals (offline engine)..."):
                        from faers_build_signals import FAERSSignalBuilder
                        builder = FAERSSignalBuilder()
                        ok = builder.process_live_data()
                    
                        if not ok:
                            st.sidebar.error("❌ FAERS signal building failed")
                            return None

                    # Load the newly created signals file
                    if not signals_path.exists():
                        st.sidebar.error(f"❌ Signals file not found: {signals_path}")
                        return None

                    data = pd.read_csv(signals_path)
            
                self.metadata = {
                    "source": "faers_live_signals",
                    "quarters": quarter_labels,
                    "date_range": f"{st.session_state.faers_start} to {st.session_state.faers_end}",
                    "load_time": datetime.now().isoformat(),
                    "record_count": len(data),
                    "signals_file": str(signals_path),
                }
            
                st.sidebar.success(
                    f"✅ Loaded {len(data):,} FAERS signals from {signals_path.name}"
                )
                logger.info(f"FAERS signals loaded: {len(data)} rows from {signals_path}")
            
                # Reset configuration for next selection
                st.session_state.faers_configured = False
                return data

            except Exception as e:
                st.sidebar.error(f"❌ FAERS signals flow failed: {str(e)}")
                logger.error("FAERS signals flow error", exc_info=True)
                return None

        # If button not clicked yet
        st.sidebar.info("👆 Click 'Build / Load FAERS Signals' to proceed")
        return None

    def _compute_quarters(self, start_date, end_date) -> List[Tuple[int, int]]:
        """
        Compute all quarters overlapping the date range
        Returns: List of (year, quarter) tuples
        """
        # Convert date objects to datetime for consistent comparison
        if not isinstance(start_date, datetime):
            start_date = datetime.combine(start_date, datetime.min.time())
        if not isinstance(end_date, datetime):
            end_date = datetime.combine(end_date, datetime.min.time())

        quarters: List[Tuple[int, int]] = []
        current = start_date

        while current <= end_date:
            year = current.year
            month = current.month
            quarter = (month - 1) // 3 + 1

            if (year, quarter) not in quarters:
                quarters.append((year, quarter))

            # Move to next quarter
            if quarter == 4:
                current = datetime(year + 1, 1, 1)
            else:
                current = datetime(year, quarter * 3 + 1, 1)

        return quarters

    def _read_file(self, file_path: str) -> pd.DataFrame:
        """Read CSV, Parquet, or SQLite file"""
        path = Path(file_path)
        if path.suffix == ".csv":
            return pd.read_csv(file_path)
        elif path.suffix == ".parquet":
            return pd.read_parquet(file_path)
        elif path.suffix == ".db":
            import sqlite3
            conn = sqlite3.connect(file_path)
            return pd.read_sql_query("SELECT * FROM icsr", conn)
        else:
            raise ValueError(f"Unsupported file type: {path.suffix}")

    def _read_uploaded_file(self, uploaded_file) -> pd.DataFrame:
        """Read uploaded file from Streamlit uploader"""
        if uploaded_file.name.endswith(".csv"):
            return pd.read_csv(uploaded_file)
        elif uploaded_file.name.endswith(".parquet"):
            return pd.read_parquet(uploaded_file)
        else:
            raise ValueError(f"Unsupported file type: {uploaded_file.name}")

    def _validate_icsr_schema(self, data: pd.DataFrame):
        """Basic validation of ICSR dataset schema"""
        required_cols = ["case_id", "drug_name", "event_term"]
        missing = [col for col in required_cols if col not in data.columns]
        
        if missing:
            raise ValueError(f"Missing required columns: {missing}")
        if len(data) == 0:
            raise ValueError("Dataset is empty")

    def get_metadata(self) -> dict:
        """Return metadata about loaded dataset"""
        return self.metadata

    def save_metadata(self):
        """Save metadata to disk for audit trail"""
        metadata_file = config.DATA_DIR / "last_load_metadata.json"
        with open(metadata_file, "w") as f:
            json.dump(self.metadata, f, indent=2)
        logger.info(f"Metadata saved to {metadata_file}")

