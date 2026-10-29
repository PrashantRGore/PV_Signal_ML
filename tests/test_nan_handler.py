import pytest
import pandas as pd
import numpy as np
from src.ml.nan_handler.denoising_autoencoder import DenoisingAutoencoder
from src.causal.causal_features import CausalFeatureEngineer


def test_mice_imputer():
    """Verify that MICE IterativeImputer imputes missing values with multivariate logic."""
    df = pd.DataFrame({
        "age": [20.0, 30.0, np.nan, 50.0, 60.0, 70.0],
        "weight": [55.0, 65.0, 75.0, np.nan, 95.0, 105.0],
        "sbp": [110.0, 120.0, 130.0, 140.0, np.nan, 160.0]
    })

    imputer = DenoisingAutoencoder(strategy="mice")
    imputed = imputer.fit_transform(df)

    # Must have 0 NaNs after imputation
    assert imputed.isna().sum().sum() == 0
    # Values should be within plausible ranges
    assert 20.0 <= imputed.loc[2, "age"] <= 60.0
    assert 60.0 <= imputed.loc[3, "weight"] <= 100.0


def test_pytorch_dae_imputer():
    """Verify that PyTorch Tabular DAE trains and reconstructs corrupted tabular data."""
    try:
        import torch
    except ImportError:
        pytest.skip("PyTorch not installed in this environment; skipping PyTorch DAE test")

    np.random.seed(42)
    # Generate synthetic correlated clinical data
    x1 = np.random.normal(50, 10, 100)
    x2 = x1 * 1.5 + np.random.normal(0, 2, 100)
    df = pd.DataFrame({"feat_a": x1, "feat_b": x2})

    imputer = DenoisingAutoencoder(strategy="dae")
    imputer.fit(df, epochs=5)

    assert imputer.strategy == "dae"
    assert imputer.dae_net is not None

    # Test with injected NaNs
    test_df = df.copy()
    test_df.loc[0:5, "feat_a"] = np.nan
    reconstructed = imputer.transform(test_df)

    assert reconstructed.isna().sum().sum() == 0
    assert len(reconstructed) == len(test_df)


def test_deterministic_pubmed_cooccurrence():
    """Verify that _get_pubmed_cooccurrence gives identical deterministic results across runs."""
    engine = CausalFeatureEngineer()
    count_1 = engine._get_pubmed_cooccurrence("IBUPROFEN", "GASTRIC ULCER")
    count_2 = engine._get_pubmed_cooccurrence("IBUPROFEN", "GASTRIC ULCER")

    assert count_1 == count_2
    assert isinstance(count_1, int)
    assert 0 <= count_1 <= 100
