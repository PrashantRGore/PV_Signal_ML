import pandas as pd


def compute_enhanced_signals(pairs: pd.DataFrame, signals_basic: pd.DataFrame) -> pd.DataFrame:
    """
    Take raw drug-event pairs and the basic disproportionality table,
    and return an enhanced signals table with 50+ engineered features.
    This keeps the original app behaviour but adds extra columns for ML.
    """

    signals = signals_basic.copy()

    # Example aggregations; keep consistent with your earlier 54-feature script
    # Demographics (per drug-event)
    demo = (
        pairs.groupby(["drug_name", "event_term"])
        .agg(
            age_mean=("age", "mean"),
            age_std=("age", "std"),
            weight_mean=("weight", "mean"),
            weight_std=("weight", "std"),
        )
        .reset_index()
    )

    # Temporal
    tto = (
        pairs.groupby(["drug_name", "event_term"])
        .agg(
            ttp_min=("ttp_days", "min"),
            ttp_max=("ttp_days", "max"),
            ttp_mean=("ttp_days", "mean"),
            ttp_median=("ttp_days", "median"),
            ttp_std=("ttp_days", "std"),
        )
        .reset_index()
    )

    # Seriousness / outcome proportions
    seriousness = (
        pairs.groupby(["drug_name", "event_term", "serious"])["case_id"]
        .nunique()
        .groupby(level=[0, 1])
        .apply(lambda s: s / s.sum())
        .unstack(fill_value=0)
        .add_prefix("serious_prop_")
        .reset_index()
    )

    outcome = (
        pairs.groupby(["drug_name", "event_term", "outcome"])["case_id"]
        .nunique()
        .groupby(level=[0, 1])
        .apply(lambda s: s / s.sum())
        .unstack(fill_value=0)
        .add_prefix("outcome_prop_")
        .reset_index()
    )

    # Concomitant drugs count
    pairs["has_concomitant"] = pairs["concomitant_drugs"].notna()
    concom = (
        pairs.groupby(["drug_name", "event_term"])["has_concomitant"]
        .sum()
        .reset_index()
        .rename(columns={"has_concomitant": "n_concomitant_reports"})
    )

    # Lab means
    labs = (
        pairs.groupby(["drug_name", "event_term"])
        .agg(
            alt_mean=("alt", "mean"),
            ast_mean=("ast", "mean"),
            bilirubin_mean=("bilirubin", "mean"),
            creatinine_mean=("creatinine", "mean"),
            bun_mean=("bun", "mean"),
        )
        .reset_index()
    )

    # Merge all blocks onto basic signals
    key = ["drug_name", "event_term"]
    for block in [demo, tto, seriousness, outcome, concom, labs]:
        signals = signals.merge(block, on=key, how="left")

    return signals
