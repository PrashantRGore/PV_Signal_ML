"""
ICSR to Signal Converter
Converts Individual Case Safety Reports (ICSRs) to aggregated drug-event signals
with rigorous vectorized disproportionality statistics (PRR, ROR, Chi²) and
zero-cost OMOP benchmark ground-truth labeling.

Legal Attribution:
- Disproportionality: Evans SJ, Waller PC, Davis S. Use of proportional reporting
  ratios (PRRs) for signal detection from spontaneous adverse drug reaction reports.
  Pharmacoepidemiol Drug Saf 2001; 10(6): 483-486.
- OMOP Benchmark: Ryan PB et al. Empirical performance evaluation of methods for risk
  identification in healthcare data. Stat Med 2013; 32(30): 5601-5613.
"""

from typing import Dict, List, Optional, Tuple
import pandas as pd
import numpy as np

from src.utils.logger import setup_logger
from src.stats_engine.disproportionality import DisproportionalityAnalysis
from src.data.standards.omop_benchmark import OMOPBenchmark

logger = setup_logger(__name__)


class ICSRToSignalConverter:
    """
    Convert ICSR data (individual case reports) to signal statistics.

    Input: ICSRs with columns like [case_id, suspect_drug, event_pt, age, sex, ...]
    Output: Aggregated signals with [drug_name, event_term, case_count, prr, ror, chi2, label, ground_truth_label]
    """

    def __init__(self, min_cases: int = 3, benchmark_only: bool = False):
        """
        Initialize converter.

        Args:
            min_cases: Minimum case count for a drug-event pair to be considered.
            benchmark_only: If True, keep only verified OMOP benchmark control pairs.
        """
        self.min_cases = min_cases
        self.benchmark_only = benchmark_only
        self.stats_engine = DisproportionalityAnalysis(min_cases=min_cases)

    def convert(self, icsr_df: pd.DataFrame) -> pd.DataFrame:
        """
        Convert ICSR dataframe to signal statistics dataframe.

        Args:
            icsr_df: DataFrame with ICSR data.

        Returns:
            DataFrame with signal statistics (canonical schema compatible).
        """
        logger.info(f"Converting {len(icsr_df):,} ICSRs to drug-event signals...")

        # Identify drug and event columns
        drug_col = self._find_drug_column(icsr_df)
        event_col = self._find_event_column(icsr_df)

        if not drug_col or not event_col:
            raise ValueError(
                f"Cannot find drug/event columns. "
                f"Available: {list(icsr_df.columns)}"
            )

        logger.info(f"Using: drug_col='{drug_col}', event_col='{event_col}'")

        # Standardize and clean input columns
        icsr_clean = icsr_df.copy()
        icsr_clean["drug_name"] = icsr_clean[drug_col].astype(str).str.upper().str.strip()
        icsr_clean["event_term"] = icsr_clean[event_col].astype(str).str.upper().str.strip()

        # Remove invalid or blank entries
        invalid_tokens = {"", "NAN", "NONE", "UNKNOWN", "NULL"}
        icsr_clean = icsr_clean[
            (~icsr_clean["drug_name"].isin(invalid_tokens)) &
            (~icsr_clean["event_term"].isin(invalid_tokens))
        ]

        if len(icsr_clean) == 0:
            logger.warning("No valid ICSR records after cleaning")
            return pd.DataFrame(columns=[
                "drug_name", "event_term", "case_count", "prr", "ror", "chi2",
                "label", "ground_truth_label", "control_status"
            ])

        # Compute vectorized disproportionality using the validated stats engine
        signals = self.stats_engine.compute_signals(icsr_clean)

        # Apply OMOP reference benchmark labeling to decouple labels from PRR/Chi²
        signals = OMOPBenchmark.label_dataframe(
            signals,
            drug_col="drug_name",
            event_col="event_term",
            benchmark_only=self.benchmark_only
        )

        # Statistical exploratory screening signal (Evans criteria: PRR >= 2, Chi² >= 4, N >= 3)
        signals["is_evans_signal"] = (
            (signals["prr"] >= 2.0) &
            (signals["chi2"] >= 4.0) &
            (signals["case_count"] >= self.min_cases)
        ).astype(int)

        # Standard canonical 'label' column:
        # Prioritize true ground truth label from OMOP benchmark controls.
        # If unconfirmed, label as 0 to avoid target leakage.
        signals["label"] = signals["ground_truth_label"].fillna(0).astype(int)

        logger.info(
            f"Generated {len(signals):,} drug-event signals "
            f"({(signals['control_status'] == 'POSITIVE_CONTROL').sum():,} positive controls, "
            f"{(signals['control_status'] == 'NEGATIVE_CONTROL').sum():,} negative controls, "
            f"{signals['is_evans_signal'].sum():,} Evans statistical candidates)"
        )

        return signals

    def _find_drug_column(self, df: pd.DataFrame) -> Optional[str]:
        """Find drug column in ICSR data."""
        candidates = [
            "drug_name", "suspect_drug", "drugname",
            "prod_ai", "medicinalproduct", "drug"
        ]

        for col in candidates:
            if col in df.columns:
                return col

        df_cols_lower = {c.lower(): c for c in df.columns}
        for col in candidates:
            if col.lower() in df_cols_lower:
                return df_cols_lower[col.lower()]

        return None

    def _find_event_column(self, df: pd.DataFrame) -> Optional[str]:
        """Find event column in ICSR data."""
        candidates = [
            "event_term", "event_pt", "reaction", "pt",
            "preferred_term", "event", "adverse_event"
        ]

        for col in candidates:
            if col in df.columns:
                return col

        df_cols_lower = {c.lower(): c for c in df.columns}
        for col in candidates:
            if col.lower() in df_cols_lower:
                return df_cols_lower[col.lower()]

        return None
