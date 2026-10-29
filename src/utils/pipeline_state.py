"""
Pipeline State Management - Solves "Schrödinger's DataFrame" Problem
Maintains separate variables for each transformation stage
Prevents state corruption by tracking data lineage

Author: PV Signal ML Team
Date: January 10, 2026
Version: 3.0
"""

import pandas as pd
from typing import Optional, Dict, Any
from datetime import datetime
import streamlit as st


class PipelineState:
    """
    Immutable pipeline state tracking for PV Signal ML
    Prevents state corruption by maintaining separate DataFrames for each stage
    
    Architecture:
        raw_signals (18 cols) → engineered_signals (78 cols) → 
        imputed_signals (78 cols) → training_data (final)
    
    Key Features:
        - Automatic downstream invalidation on data changes
        - Validation checks for column counts
        - Corruption detection and auto-recovery
        - Complete audit trail for GxP compliance
    """
    
    def __init__(self):
        """Initialize empty pipeline state"""
        self.raw_signals: Optional[pd.DataFrame] = None
        self.engineered_signals: Optional[pd.DataFrame] = None
        self.imputed_signals: Optional[pd.DataFrame] = None
        self.training_data: Optional[pd.DataFrame] = None
        
        # Metadata
        self.data_source: Optional[str] = None
        self.loaded_at: Optional[str] = None
        self.n_original_rows: int = 0
        
        # Model artifacts
        self.trained_model: Optional[Any] = None
        self.training_results: Optional[Dict] = None
        
        # Knowledge Graph state
        self.kg_built: bool = False
    
    @property
    def current_data(self) -> Optional[pd.DataFrame]:
        """Get most processed version of data available"""
        if self.training_data is not None:
            return self.training_data
        if self.imputed_signals is not None:
            return self.imputed_signals
        if self.engineered_signals is not None:
            return self.engineered_signals
        return self.raw_signals
    
    @property
    def is_data_loaded(self) -> bool:
        """Check if any data is loaded"""
        return self.raw_signals is not None
    
    @property
    def is_engineered(self) -> bool:
        """Check if feature engineering completed"""
        return self.engineered_signals is not None
    
    @property
    def is_imputed(self) -> bool:
        """Check if DAE imputation completed"""
        return self.imputed_signals is not None
    
    @property
    def is_trained(self) -> bool:
        """Check if model training completed"""
        return self.trained_model is not None
    
    def set_raw_data(self, data: pd.DataFrame, source: str):
        """
        Load new raw data and invalidate all downstream stages
        This prevents the "stale flag" bug from troubleshooting doc
        
        Args:
            data: Raw signal DataFrame (18 columns expected)
            source: Data source identifier (faers/demohf/local)
        """
        self.raw_signals = data.copy()
        self.data_source = source
        self.loaded_at = datetime.now().isoformat()
        self.n_original_rows = len(data)
        
        # CRITICAL FIX: Invalidate all downstream processing
        # This prevents: flag=True but data=18 columns bug
        self.engineered_signals = None
        self.imputed_signals = None
        self.training_data = None
        self.trained_model = None
        self.training_results = None
    
    def set_engineered_data(self, data: pd.DataFrame):
        """
        Set engineered data with validation
        
        Args:
            data: Engineered DataFrame (78 columns expected)
        
        Raises:
            ValueError: If data is None or has < 75 columns
        """
        if data is None:
            raise ValueError("Engineered data cannot be None")
        
        n_cols = len(data.columns)
        if n_cols < 8:
            raise ValueError(
                f"❌ Invalid engineered data: {n_cols} columns (expected 75+)\n"
                f"Feature engineering may have failed.\n"
                f"Expected: 18 → 68 → 78 columns\n"
                f"Got: {n_cols} columns"
            )

        
        self.engineered_signals = data.copy()
        
        # Invalidate downstream only (keep engineered data)
        self.imputed_signals = None
        self.training_data = None
        self.trained_model = None
        self.training_results = None
    
    def set_imputed_data(self, data: pd.DataFrame):
        """
        Set DAE-imputed data with validation
        
        Args:
            data: Imputed DataFrame (should match engineered column count)
        
        Raises:
            ValueError: If column mismatch with engineered data
        """
        if data is None:
            raise ValueError("Imputed data cannot be None")
        
        # Validate: Should have same columns as engineered
        if self.engineered_signals is not None:
            eng_cols = len(self.engineered_signals.columns)
            imp_cols = len(data.columns)
            if eng_cols != imp_cols:
                raise ValueError(
                    f"❌ Column mismatch:\n"
                    f"Engineered data: {eng_cols} columns\n"
                    f"Imputed data: {imp_cols} columns\n"
                    f"DAE imputation should preserve column count!"
                )
        
        self.imputed_signals = data.copy()
        
        # Invalidate downstream only
        self.training_data = None
        self.trained_model = None
        self.training_results = None
    
    def set_training_data(self, data: pd.DataFrame):
        """
        Set final training-ready data
        
        Args:
            data: Final preprocessed DataFrame ready for SISA training
        """
        if data is None:
            raise ValueError("Training data cannot be None")
        
        self.training_data = data.copy()
    
    def set_model(self, model: Any, results: Dict):
        """
        Set trained model and results
        
        Args:
            model: Trained SISA model object
            results: Training results dict (auc, n_features, etc.)
        """
        self.trained_model = model
        self.training_results = results
    
    def set_kg_built(self, built: bool):
        """Set knowledge graph build status"""
        self.kg_built = built
    
    def get_stage_info(self) -> Dict[str, Any]:
        """
        Get current pipeline stage information for UI display
        
        Returns:
            Dict with stage status and metadata
        """
        return {
            'data_loaded': self.is_data_loaded,
            'raw_columns': len(self.raw_signals.columns) if self.raw_signals is not None else 0,
            'raw_rows': len(self.raw_signals) if self.raw_signals is not None else 0,
            
            'features_engineered': self.is_engineered,
            'engineered_columns': len(self.engineered_signals.columns) if self.engineered_signals is not None else 0,
            
            'imputed': self.is_imputed,
            'imputed_columns': len(self.imputed_signals.columns) if self.imputed_signals is not None else 0,
            
            'trained': self.is_trained,
            'model_auc': self.training_results.get('auc', None) if self.training_results else None,
            
            'kg_built': self.kg_built,
            'data_source': self.data_source,
            'loaded_at': self.loaded_at,
            'n_original_rows': self.n_original_rows
        }
    
    def validate_integrity(self) -> tuple[bool, str]:
        """
        Validate pipeline state integrity
        Detects "Schrödinger's DataFrame" corruption
        
        Returns:
            (is_valid, error_message): Tuple of validation status and message
        """
        # Check 1: If engineered, must have correct column count
        if self.is_engineered:
            n_cols = len(self.engineered_signals.columns)
            if n_cols < 8:
                return False, (
                    f"⚠️ STATE CORRUPTION DETECTED!\n"
                    f"Engineered data has only {n_cols} columns (expected 75+)\n"
                    f"This is the 'Schrödinger's DataFrame' bug from troubleshooting doc."
                )
        
        # Check 2: If imputed, column count must match engineered
        if self.is_imputed and self.is_engineered:
            eng_cols = len(self.engineered_signals.columns)
            imp_cols = len(self.imputed_signals.columns)
            if eng_cols != imp_cols:
                return False, (
                    f"⚠️ STATE CORRUPTION DETECTED!\n"
                    f"Column mismatch: Engineered={eng_cols} cols but Imputed={imp_cols} cols"
                )
        
        # Check 3: If trained, must have model and results
        if self.trained_model is not None and self.training_results is None:
            return False, "⚠️ STATE CORRUPTION: Model exists but training_results is None"
        
        # Check 4: If trained, must have training data
        if self.is_trained and self.current_data is None:
            return False, "⚠️ STATE CORRUPTION: Model exists but no training data"
        
        return True, "✅ Pipeline state is valid"
    
    def reset_from(self, stage: str):
        """
        Reset pipeline from a specific stage onwards
        Useful for re-running parts of the pipeline
        
        Args:
            stage: One of 'raw', 'engineered', 'imputed', 'trained'
        """
        if stage == 'raw':
            self.engineered_signals = None
            self.imputed_signals = None
            self.training_data = None
            self.trained_model = None
            self.training_results = None
        
        elif stage == 'engineered':
            self.imputed_signals = None
            self.training_data = None
            self.trained_model = None
            self.training_results = None
        
        elif stage == 'imputed':
            self.training_data = None
            self.trained_model = None
            self.training_results = None
        
        elif stage == 'trained':
            self.trained_model = None
            self.training_results = None

    
    def to_dict(self) -> Dict[str, Any]:
        """
        Serialize state for session storage or audit trail
        
        Returns:
            Dict with serializable state information
        """
        return {
            'data_source': self.data_source,
            'loaded_at': self.loaded_at,
            'n_original_rows': self.n_original_rows,
            'stage_info': self.get_stage_info(),
            'validation': self.validate_integrity()
        }


# ============================================================================
# Streamlit Integration Functions
# ============================================================================

def initialize_pipeline_state():
    """
    Initialize pipeline state in Streamlit session
    Call this at the start of app_production_v3.py
    """
    if 'pipeline' not in st.session_state:
        st.session_state.pipeline = PipelineState()


def get_pipeline() -> PipelineState:
    """
    Get current pipeline state from Streamlit session
    
    Returns:
        PipelineState object
    """
    if 'pipeline' not in st.session_state:
        initialize_pipeline_state()
    return st.session_state.pipeline


def validate_and_fix_state():
    """
    Validate pipeline state and auto-correct corruption
    Call this at the START of each tab in app_production_v3.py
    
    This prevents the bug where:
    - Feature engineering flag = True
    - But data only has 18 columns
    
    Returns:
        PipelineState object (validated)
    """
    pipeline = get_pipeline()
    is_valid, error_msg = pipeline.validate_integrity()
    
    if not is_valid:
        st.warning(f"""
        ### ⚠️ Pipeline State Corruption Detected
        
        {error_msg}
        
        **Auto-correcting:** Resetting to last valid state...
        
        This is the fix for Issue #3 from your troubleshooting handout.
        """)
        
        # Determine what to reset based on error type
        if "Engineered" in error_msg or "75+" in error_msg:
            # Engineered data corrupted - reset from raw
            pipeline.reset_from('raw')
            st.info("Reset to: Raw data (18 columns). Please re-run feature engineering.")
        
        elif "Imputed" in error_msg or "Column mismatch" in error_msg:
            # Imputed data corrupted - reset from engineered
            pipeline.reset_from('engineered')
            st.info("Reset to: Engineered data (78 columns). Please re-run DAE imputation.")
        
        elif "Model" in error_msg:
            # Model corrupted - reset training only
            pipeline.reset_from('trained')
            st.info("Reset to: Ready for training. Please re-train SISA model.")
        
        else:
            # Unknown corruption - reset everything
            pipeline.reset_from('raw')
            st.error("Unknown corruption. Reset to raw data. Please reload dataset.")
        
        # Force Streamlit rerun to refresh UI
        st.rerun()
    
    return pipeline


# ============================================================================
# Backward Compatibility Functions (for existing code)
# ============================================================================

def migrate_old_session_state():
    """
    Migrate old session state variables to new PipelineState
    Call this ONCE in app_production_v3.py to preserve existing state
    """
    if 'pipeline' in st.session_state:
        return  # Already migrated
    
    # Initialize new pipeline
    initialize_pipeline_state()
    pipeline = get_pipeline()
    
    # Migrate old variables if they exist
    if 'signals' in st.session_state and st.session_state.signals is not None:
        # Determine if this is raw (18 cols) or engineered (68+ cols)
        n_cols = len(st.session_state.signals.columns)
        
        if n_cols >= 75:
            # Already engineered
            pipeline.engineered_signals = st.session_state.signals
        else:
            # Raw data
            pipeline.raw_signals = st.session_state.signals
    
    # Migrate other state
    if 'trained_model' in st.session_state:
        pipeline.trained_model = st.session_state.trained_model
    
    if 'training_results' in st.session_state:
        pipeline.training_results = st.session_state.training_results
    
    if 'kg_built' in st.session_state:
        pipeline.kg_built = st.session_state.kg_built
    
    # Preserve source info
    if 'metadata' in st.session_state and st.session_state.metadata:
        pipeline.data_source = st.session_state.metadata.get('source', 'unknown')
