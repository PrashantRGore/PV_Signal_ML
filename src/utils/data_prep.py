"""
Data preparation utilities:
- prepare_data_for_signal_detection (existing behaviour placeholder)
- normalize_signals: canonical schema across all sources
"""

from typing import Literal

import numpy as np
import pandas as pd


def prepare_data_for_signal_detection(df: pd.DataFrame) -> pd.DataFrame:
    """
    Existing hook for any pre-signal-detection cleaning.
    Currently acts as a pass-through to avoid changing behaviour.
    """
    return df


def normalize_signals(
    signals: pd.DataFrame,
    source: Literal["local", "demo_hf", "faers_live", "faers_live_signals"] = "local",
) -> pd.DataFrame:
    """
    Normalize signal tables from different sources into a canonical schema:

    Required canonical columns:
        - drug_name      : str
        - event_term     : str
        - case_count     : int
        - prr            : float
        - chi2           : float
        - label          : int (0/1)

    For richer sources (e.g. 1M demo), all extra columns are preserved.
    For FAERS, source-specific names are mapped to canonical ones.
    """

    df = signals.copy()

    # --- 1. Map identifiers ---
    # Try common columns first; fall back to FAERS names if needed.
    if "drug_name" not in df.columns:
        if "DRUGNAMENORM" in df.columns:
            df["drug_name"] = df["DRUGNAMENORM"]
        elif "drug" in df.columns:
            df["drug_name"] = df["drug"]
        else:
            raise ValueError("normalize_signals: could not find a drug_name column.")

    if "event_term" not in df.columns:
        if "pt" in df.columns:
            df["event_term"] = df["pt"]
        elif "event" in df.columns:
            df["event_term"] = df["event"]
        else:
            raise ValueError("normalize_signals: could not find an event_term column.")

    # --- 2. Map quantitative stats (case_count, prr, chi2) ---
    # case_count
    if "case_count" not in df.columns:
        if "A_6m" in df.columns:
            df["case_count"] = df["A_6m"]
        elif "A" in df.columns:
            df["case_count"] = df["A"]
        else:
            raise ValueError(
                "normalize_signals: could not find case_count (expected case_count / A_6m / A)."
            )

    # PRR
    if "prr" not in df.columns:
        if "prr_6m" in df.columns:
            df["prr"] = df["prr_6m"]
        elif "PRR" in df.columns:
            df["prr"] = df["PRR"]
        else:
            raise ValueError(
                "normalize_signals: could not find PRR (expected prr / prr_6m / PRR)."
            )

    # Chi-square
    if "chi2" not in df.columns:
        if "chi2_6m" in df.columns:
            df["chi2"] = df["chi2_6m"]
        elif "chi_sq" in df.columns:
            df["chi2"] = df["chi_sq"]
        else:
            raise ValueError(
                "normalize_signals: could not find chi-square (expected chi2 / chi2_6m / chi_sq)."
            )

    # --- 3. Map label column ---
    # Priority: explicit label_prr -> is_signal_prr -> is_signal -> fallback 0
    if "label" not in df.columns:
        if "label_prr" in df.columns:
            df["label"] = df["label_prr"].astype(int)
        elif "is_signal_prr" in df.columns:
            df["label"] = df["is_signal_prr"].astype(int)
        elif "is_signal" in df.columns:
            df["label"] = df["is_signal"].astype(int)
        else:
            # If no label exists (e.g. some local exploratory datasets), default to 0
            df["label"] = 0

    # --- 4. Type enforcement ---
    df["drug_name"] = df["drug_name"].astype(str)
    df["event_term"] = df["event_term"].astype(str)
    df["case_count"] = df["case_count"].astype(int)
    df["prr"] = df["prr"].astype(float)
    df["chi2"] = df["chi2"].astype(float)
    df["label"] = df["label"].astype(int)

    # --- 5. Final safety: no inf values in canonical stats ---
    for col in ["prr", "chi2"]:
        df[col] = df[col].replace([np.inf, -np.inf], np.nan)

    # Use conservative caps in case caller did not already cap
    prr_finite = df["prr"].dropna()
    if len(prr_finite) == 0:
        prr_cap = 50.0
    else:
        prr_cap = float(min(prr_finite.quantile(0.99), 50.0))
    df["prr"] = df["prr"].fillna(prr_cap).clip(upper=prr_cap)

    chi2_finite = df["chi2"].dropna()
    if len(chi2_finite) == 0:
        chi2_cap = 1000.0
    else:
        chi2_cap = float(min(chi2_finite.quantile(0.99), 1000.0))
    df["chi2"] = df["chi2"].fillna(chi2_cap).clip(upper=chi2_cap)

    return df
