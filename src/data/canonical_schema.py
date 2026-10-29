"""
Canonical Schema - HOTFIX for case-insensitive column matching
"""

import pandas as pd
import numpy as np
from typing import List, Optional, Dict, Any
import warnings

from src.utils.logger import setup_logger

logger = setup_logger(__name__)
warnings.filterwarnings("ignore", category=FutureWarning)


class CanonicalSchema:
    """
    Defines the canonical schema for PV signal data across all sources.
    
    TIER 1 (Required): Minimal features for basic signal detection
    TIER 2 (Optional): Extended features from Demo 1M dataset
    """

    # Tier 1: Core required columns
    TIER1_REQUIRED = [
        "drug_name",      # Drug product name
        "event_term",     # Adverse event term (MedDRA PT)
        "case_count",     # Number of co-occurring cases (A in 2x2 table)
        "prr",           # Proportional Reporting Ratio
        "ror",           # Reporting Odds Ratio
        "chi2",          # Chi-Square statistic
        "label",         # Binary label (1=signal, 0=non-signal)
    ]

    # Tier 2: Extended features are dynamically detected (numeric columns not in Tier 1)

    @staticmethod
    def validate(df: pd.DataFrame, allow_missing_label: bool = False) -> None:
        """
        Validate that dataframe matches canonical schema.
        
        Args:
            df: DataFrame to validate
            allow_missing_label: If True, label column is optional
            
        Raises:
            ValueError: If schema validation fails
        """
        required = CanonicalSchema.TIER1_REQUIRED.copy()
        if allow_missing_label:
            required.remove("label")

        missing = [col for col in required if col not in df.columns]
        if missing:
            raise ValueError(
                f"Missing required canonical columns: {missing}. "
                f"Available: {list(df.columns)}"
            )

        # Check data types
        for col in ["prr", "ror", "chi2", "case_count"]:
            if col in df.columns and not pd.api.types.is_numeric_dtype(df[col]):
                raise ValueError(f"Column '{col}' must be numeric")

        if "label" in df.columns and not pd.api.types.is_integer_dtype(df["label"]):
            raise ValueError("Column 'label' must be integer (0 or 1)")

        # Check for inf/nan in critical columns
        # STRICT post-imputation validation (NO SILENT PASS)
        for col in ["prr", "ror", "chi2"]:
            if col in df.columns:
                if np.isinf(df[col]).any():
                    raise ValueError(
                        f"❌ CRITICAL: Column '{col}' contains inf values "
                        "after canonicalization / imputation"
                    )
                if np.isnan(df[col]).any():
                    logger.critical(
                        f"⚠️⚠️⚠️ Column '{col}' contains NaN values after processing. "
                        "DAE imputation required before ML training."
                    )
        
        logger.info(f"✅ Schema validation passed for {len(df)} signals")

    @staticmethod
    def get_tier1_feature_columns() -> List[str]:
        """Get Tier 1 feature columns (excludes identifiers and label)"""
        return ["case_count", "prr", "ror", "chi2"]

    @staticmethod
    def get_tier2_columns(df: pd.DataFrame) -> List[str]:
        """
        Auto-detect Tier 2 feature columns.
        
        Tier 2 = numeric columns not in Tier 1 required set
        """
        tier1_set = set(CanonicalSchema.TIER1_REQUIRED)
        
        tier2_cols = [
            col for col in df.columns
            if col not in tier1_set
            and pd.api.types.is_numeric_dtype(df[col])
        ]
        
        return tier2_cols


def _find_column_case_insensitive(df: pd.DataFrame, patterns: List[str]) -> Optional[str]:
    """
    Find column matching any pattern (case-insensitive).
    
    Args:
        df: DataFrame to search
        patterns: List of patterns to match (e.g., ["prr", "PRR_6m"])
    
    Returns:
        Matched column name or None
    """
    df_cols_lower = {col.lower(): col for col in df.columns}
    
    for pattern in patterns:
        pattern_lower = pattern.lower()
        
        # Exact match (case-insensitive)
        if pattern_lower in df_cols_lower:
            return df_cols_lower[pattern_lower]
        
        # Partial match (e.g., "prr" matches "PRR_6m")
        for col_lower, col_original in df_cols_lower.items():
            if pattern_lower in col_lower:
                return col_original
    
    return None


def canonicalize_signals(
    df: pd.DataFrame,
    source: str,
    preserve_tier2: bool = True,
) -> pd.DataFrame:
    """
    Convert signals from any source to canonical schema.
    
    Args:
        df: Input dataframe (any source format)
        source: Data source identifier ("demo_hf", "faers_live", "local", etc.)
        preserve_tier2: If True, keep all extra numeric columns as Tier 2 features
    
    Returns:
        DataFrame in canonical schema format
    """
    logger.info(f"[{source}] Canonicalizing {len(df)} signals...")
    
    canonical = pd.DataFrame()
    
    # ===== STEP 1: Map drug_name and event_term =====
    
    drug_col = _find_column_case_insensitive(df, ["drug_name", "drug", "drugname"])
    event_col = _find_column_case_insensitive(df, ["event_term", "event", "pt", "reaction"])
    
    if not drug_col:
        raise ValueError(f"[{source}] Cannot find drug name column")
    if not event_col:
        raise ValueError(f"[{source}] Cannot find event term column")
    
    canonical["drug_name"] = df[drug_col].astype(str).str.upper().str.strip()
    canonical["event_term"] = df[event_col].astype(str).str.upper().str.strip()
    
    # ===== STEP 2: Map case_count =====
    
    case_col = _find_column_case_insensitive(
        df, ["case_count", "count", "a", "a_6m", "cases", "n"]
    )
    
    if not case_col:
        raise ValueError(f"[{source}] Cannot find case count column")
    
    canonical["case_count"] = pd.to_numeric(df[case_col], errors="coerce").fillna(0).astype(int)
    
    # ===== STEP 3: Map PRR, ROR, Chi2 =====
    
    prr_col = _find_column_case_insensitive(df, ["prr", "prr_6m", "proportional_reporting_ratio"])
    ror_col = _find_column_case_insensitive(df, ["ror", "ror_6m", "reporting_odds_ratio"])
    chi2_col = _find_column_case_insensitive(df, ["chi2", "chisquare", "chi_square", "chisquare_6m"])
    
    if not prr_col:
        raise ValueError(f"[{source}] Cannot find PRR column. Available: {list(df.columns)}")
    if not ror_col:
        raise ValueError(f"[{source}] Cannot find ROR column")
    if not chi2_col:
        raise ValueError(f"[{source}] Cannot find Chi2 column")
    
    canonical["prr"] = pd.to_numeric(df[prr_col], errors="coerce").fillna(0.0)
    canonical["ror"] = pd.to_numeric(df[ror_col], errors="coerce").fillna(0.0)
    canonical["chi2"] = pd.to_numeric(df[chi2_col], errors="coerce").fillna(0.0)
    
    # Replace inf values with 999.0
    canonical["prr"] = canonical["prr"].replace([np.inf, -np.inf], 999.0)
    canonical["ror"] = canonical["ror"].replace([np.inf, -np.inf], 999.0)
    canonical["chi2"] = canonical["chi2"].replace([np.inf, -np.inf], 999.0)
    
    # ===== STEP 4: Map label =====
    
    label_col = _find_column_case_insensitive(df, ["label", "is_signal", "ground_truth"])
    
    if label_col:
        canonical["label"] = pd.to_numeric(df[label_col], errors="coerce").fillna(0).astype(int)
    else:
        # Use PRR-based heuristic
        canonical["label"] = (
            (canonical["prr"] >= 2.0) & (canonical["chi2"] >= 4.0)
        ).astype(int)
        logger.info(f"[{source}] No label column found, using PRR-based heuristic")
    
    # ===== STEP 5: Preserve Tier 2 features =====
    
    if preserve_tier2:
        tier1_mapped_cols = {drug_col, event_col, case_col, prr_col, ror_col, chi2_col}
        if label_col:
            tier1_mapped_cols.add(label_col)
        
        tier2_cols = [
            col for col in df.columns
            if col not in tier1_mapped_cols
            and pd.api.types.is_numeric_dtype(df[col])
        ]
        
        if tier2_cols:
            for col in tier2_cols:
                canonical[col] = pd.to_numeric(df[col], errors="coerce").fillna(0.0)
            
            logger.info(f"[{source}] Preserved {len(tier2_cols)} Tier 2 features")
    
    # ===== STEP 6: Validate canonical schema =====
    
    try:
        CanonicalSchema.validate(canonical)
    except ValueError as e:
        logger.error(f"[{source}] Canonical schema validation failed: {e}")
        raise
    
    tier2_count = len(CanonicalSchema.get_tier2_columns(canonical))
    
    logger.info(
        f"[{source}] ✅ Canonicalized {len(canonical)} signals "
        f"(Tier 1: {len(CanonicalSchema.TIER1_REQUIRED)} + Tier 2: {tier2_count})"
    )
    
    return canonical


def summarize_schema(df: pd.DataFrame) -> Dict[str, Any]:
    """
    Generate summary of canonical schema structure.
    
    Args:
        df: Canonical dataframe
    
    Returns:
        Dictionary with schema metadata
    """
    tier1_cols = CanonicalSchema.TIER1_REQUIRED
    tier2_cols = CanonicalSchema.get_tier2_columns(df)
    
    summary = {
        "total_signals": len(df),
        "tier1_required": tier1_cols,
        "tier1_feature_count": len(CanonicalSchema.get_tier1_feature_columns()),
        "tier2_features": tier2_cols,
        "tier2_feature_count": len(tier2_cols),
        "has_tier2": len(tier2_cols) > 0,
        "positive_signals": int(df["label"].sum()) if "label" in df.columns else None,
        "signal_rate": float(df["label"].mean()) if "label" in df.columns else None,
    }
    
    return summary
