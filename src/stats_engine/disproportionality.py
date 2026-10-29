"""
Disproportionality Analysis - Statistical Signal Detection
NOW WITH: Canonical schema output + VECTORIZED computation + CORRECT chi2 + FIXED label logic
"""

import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency
from typing import Optional
import warnings

from src.data.canonical_schema import CanonicalSchema
from src.utils.logger import setup_logger

logger = setup_logger(__name__)
warnings.filterwarnings("ignore", category=RuntimeWarning)


class DisproportionalityAnalysis:
    """
    Compute PRR, ROR, and Chi-Square for drug-event pairs.

    OUTPUT: Always in canonical schema format (no column name ambiguity)
    """

    def __init__(self, min_cases: int = 3):
        self.min_cases = min_cases

    def compute_signals(
        self, icsr_data: pd.DataFrame, label_col: Optional[str] = None
    ) -> pd.DataFrame:
        """
        Compute disproportionality signals from raw ICSR data.

        Args:
            icsr_data: Raw ICSR data with columns [case_id, drug_name, event_term]
            label_col: Optional ground truth label column (for Demo 1M dataset)

        Returns:
            DataFrame in CANONICAL SCHEMA format with columns:
            - drug_name, event_term, case_count, prr, ror, chi2, label
            - Plus any Tier 2 features from input (preserved if present)
        """
        logger.info(f"Computing signals from {len(icsr_data):,} ICSR records")

        # Validate required columns
        required_cols = ["drug_name", "event_term"]
        missing = [c for c in required_cols if c not in icsr_data.columns]
        if missing:
            raise ValueError(f"Missing required columns: {missing}")

        # Check for Tier 2 features in input (Demo 1M case)
        tier2_features = [
            col for col in icsr_data.columns
            if col not in CanonicalSchema.TIER1_REQUIRED
            and col not in ["case_id", "drug_name", "event_term"]
            and pd.api.types.is_numeric_dtype(icsr_data[col])
        ]

        has_tier2 = len(tier2_features) > 0
        if has_tier2:
            logger.info(f"Detected {len(tier2_features)} Tier 2 features in input data")

        # Group by drug-event pairs
        if has_tier2:
            # Aggregate Tier 2 features using mean
            agg_dict = {"case_id": "count"}
            for col in tier2_features:
                agg_dict[col] = "mean"

            df = (
                icsr_data.groupby(["drug_name", "event_term"])
                .agg(agg_dict)
                .reset_index()
            )
        else:
            # Simple case count
            df = (
                icsr_data.groupby(["drug_name", "event_term"])
                .size()
                .reset_index(name="case_id")
            )

        # Rename case_id count to case_count (canonical name)
        df.rename(columns={"case_id": "case_count"}, inplace=True)

        # Filter by minimum cases
        df = df[df["case_count"] >= self.min_cases].copy()
        logger.info(f"After filtering: {len(df):,} drug-event pairs")

        if len(df) == 0:
            logger.warning("No signals found after filtering")
            return pd.DataFrame(columns=CanonicalSchema.TIER1_REQUIRED)

        # ========================================================================
        # VECTORIZED COMPUTATION - Pre-compute marginal counts ONCE
        # ========================================================================
        total_reports = icsr_data["case_id"].nunique()
        
        # Pre-compute drug marginals (b values)
        drug_counts = icsr_data.groupby("drug_name")["case_id"].nunique().to_dict()
        
        # Pre-compute event marginals (c values)
        event_counts = icsr_data.groupby("event_term")["case_id"].nunique().to_dict()
        
        # Vectorized computation
        df["drug_total"] = df["drug_name"].map(drug_counts)
        df["event_total"] = df["event_term"].map(event_counts)
        
        # Contingency table:
        #              Event    Not Event   Total
        # Drug          a         b         a+b = drug_total
        # Not Drug      c         d         c+d = not_drug_total
        # Total        a+c       b+d          N
        a = df["case_count"].values.astype(float)
        b = (df["drug_total"].values - a).astype(float)
        c = (df["event_total"].values - a).astype(float)
        d = (total_reports - a - b - c).astype(float)
        
        # PRR = [a / (a + b)] / [c / (c + d)] = [a * (c + d)] / [c * (a + b)]
        # ROR = (a / b) / (c / d) = (a * d) / (b * c)
        with np.errstate(divide='ignore', invalid='ignore'):
            # PRR: proportion of event in suspect drug vs other drugs
            prr_denom = c * (a + b)
            prr_num = a * (c + d)
            prr_values = np.where(
                (prr_denom > 0) & (c > 0) & ((a + b) > 0),
                prr_num / prr_denom,
                999.0
            )
            df["prr"] = np.nan_to_num(prr_values, nan=999.0, posinf=999.0, neginf=999.0)
            
            # ROR: reporting odds ratio
            ror_denom = b * c
            ror_num = a * d
            ror_values = np.where(
                (ror_denom > 0) & (b > 0) & (c > 0) & (d > 0),
                ror_num / ror_denom,
                999.0
            )
            df["ror"] = np.nan_to_num(ror_values, nan=999.0, posinf=999.0, neginf=999.0)

            # Vectorized Yates' corrected Chi-Square:
            # chi2 = N * (|a*d - b*c| - N/2)^2 / [(a+b)*(c+d)*(a+c)*(b+d)]
            N = a + b + c + d
            chi2_denom = (a + b) * (c + d) * (a + c) * (b + d)
            diff = np.abs(a * d - b * c) - (N / 2.0)
            diff = np.maximum(diff, 0.0)
            chi2_values = np.where(
                chi2_denom > 0,
                (N * (diff ** 2)) / chi2_denom,
                0.0
            )
            df["chi2"] = np.nan_to_num(chi2_values, nan=0.0, posinf=999.0, neginf=0.0)
        
        # Drop temporary columns
        df.drop(columns=["drug_total", "event_total"], inplace=True)
        
        # ========================================================================
        # Handle ground truth labels if available
        # ========================================================================
        if label_col and label_col in icsr_data.columns:
            logger.info(f"Using ground truth labels from '{label_col}'")

            # Aggregate labels (majority vote)
            label_agg = (
                icsr_data.groupby(["drug_name", "event_term"])[label_col]
                .apply(lambda x: (x.sum() / len(x)) > 0.5)
                .astype(int)
                .reset_index()
            )

            df = df.merge(
                label_agg, on=["drug_name", "event_term"], how="left"
            )
            df["label"] = df[label_col].fillna(0).astype(int)
            df.drop(columns=[label_col], inplace=True)
        else:
            # Use PRR-based heuristic (PRR >= 2.0 and Chi2 >= 4.0)
            # EXCLUDE PRR=999 (error sentinel value)
            df["label"] = (
                (df["prr"] >= 2.0) & (df["prr"] < 999.0) & (df["chi2"] >= 4.0)
            ).astype(int)
            logger.info("Using PRR-based heuristic for labels (PRR>=2 and <999, Chi²>=4)")

        # Ensure all canonical columns are present
        for col in CanonicalSchema.TIER1_REQUIRED:
            if col not in df.columns:
                if col == "label":
                    df["label"] = 0
                else:
                    raise ValueError(f"Failed to create canonical column: {col}")

        # Replace inf values with large finite values
        df["prr"] = df["prr"].replace([np.inf, -np.inf], 999.0)
        df["ror"] = df["ror"].replace([np.inf, -np.inf], 999.0)
        df["chi2"] = df["chi2"].replace([np.inf, -np.inf], 999.0)

        # Final validation
        CanonicalSchema.validate(df)

        logger.info(
            f"✅ Computed {len(df):,} signals in canonical schema "
            f"({'with' if has_tier2 else 'without'} Tier 2 features)"
        )

        return df
