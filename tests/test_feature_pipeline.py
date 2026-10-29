import pytest
import pandas as pd
import numpy as np
from src.ml.feature_engineering.feature_pipeline import FeatureEngineeringPipeline
from src.data.enrichment.faers_enriched_extractor import FAERSEnrichedExtractor


def test_feature_pipeline_bmi_and_polypharmacy():
    """Verify BMI calculation with explicit height handling and clean polypharmacy parsing."""
    pipeline = FeatureEngineeringPipeline()

    data = pd.DataFrame({
        "age": [25, 70, np.nan],
        "weight_kg": [70.0, 80.0, np.nan],
        "height_cm": [175.0, np.nan, np.nan],
        "concomitant_medications": [
            "ASPIRIN, METFORMIN, LISINOPRIL, ATORVASTATIN, OMEPRAZOLE",  # 5 drugs -> polypharmacy risk 1
            "none",                                                      # 0 drugs
            np.nan                                                       # 0 drugs
        ],
        "country": ["USA", "DEU", "IND"]
    })

    transformed = pipeline.transform(data)

    # Patient 0: 70kg, 175cm -> BMI = 70 / (1.75^2) = 22.857
    assert pytest.approx(transformed.loc[0, "bmi_estimated"], 0.1) == 22.86
    assert transformed.loc[0, "_data_quality_height_imputed"] == 0
    assert transformed.loc[0, "concomitant_count"] == 5
    assert transformed.loc[0, "polypharmacy_risk"] == 1

    # Patient 1: 80kg, missing height -> BMI calculated with 1.7m proxy but flagged
    assert pytest.approx(transformed.loc[1, "bmi_estimated"], 0.1) == 80.0 / (1.7 ** 2)
    assert transformed.loc[1, "_data_quality_height_imputed"] == 1
    assert transformed.loc[1, "concomitant_count"] == 0
    assert transformed.loc[1, "polypharmacy_risk"] == 0

    # Patient 2: missing weight and concomitant
    assert np.isnan(transformed.loc[2, "bmi_estimated"])
    assert transformed.loc[2, "concomitant_count"] == 0


def test_faers_merge_all_reactions():
    """Verify that FAERSEnrichedExtractor preserves all reactions per case instead of dropping them."""
    extractor = FAERSEnrichedExtractor()

    drug_df = pd.DataFrame({
        "primaryid": [101, 102],
        "drugname": ["WARFARIN", "ASPIRIN"],
        "role_cod": ["PS", "PS"],
        "dose_amt": [5, 100]
    })

    # Case 101 has 3 distinct adverse reactions
    reac_df = pd.DataFrame({
        "primaryid": [101, 101, 101, 102],
        "pt": ["HAEMORRHAGE", "EPISTAXIS", "HAEMATURIA", "GASTRIC ULCER"]
    })

    merged = extractor._merge_tables(
        demo_df=None,
        drug_df=drug_df,
        reac_df=reac_df,
        outc_df=None,
        rpsr_df=None,
        ther_df=None,
        indi_df=None
    )

    # Must preserve all 3 reactions for case 101 + 1 for case 102 = 4 total rows
    assert len(merged) == 4
    case_101_events = set(merged[merged["primaryid"] == 101]["pt"])
    assert case_101_events == {"HAEMORRHAGE", "EPISTAXIS", "HAEMATURIA"}

