import pytest
import pandas as pd
import numpy as np
from pathlib import Path
from src.ml.sisa_trainer import SISATrainer


def test_sisa_training_and_unlearning(tmp_path):
    """Test SISA ensemble training, prediction, and unlearning."""
    np.random.seed(42)
    n = 200
    
    data = pd.DataFrame({
        "case_id": [f"CASE_{i}" for i in range(n)],
        "drug_name": np.random.choice(["DrugX", "DrugY"], n),
        "event_term": np.random.choice(["Event1", "Event2"], n),
        "case_count": np.random.randint(5, 50, n),
        "prr": np.random.uniform(1.0, 5.0, n),
        "ror": np.random.uniform(1.0, 5.0, n),
        "chi2": np.random.uniform(2.0, 30.0, n),
        "label": np.random.choice([0, 1], n, p=[0.7, 0.3])
    })
    
    trainer = SISATrainer(model_dir=tmp_path / "sisa_models", n_shards=4)
    meta = trainer.train(data, n_shards=4)
    
    assert "auc" in meta
    assert len(trainer.shard_models) == 4
    
    # Test ensemble predict and predict_proba
    probas = trainer.predict_proba(data)
    assert len(probas) == n
    assert np.all((probas >= 0.0) & (probas <= 1.0))
    
    preds = trainer.predict(data)
    assert len(preds) == n
    assert set(np.unique(preds)).issubset({0, 1})
    
    # Test unlearning
    target_case = "CASE_10"
    res = trainer.unlearn(target_case)
    
    assert res["status"] == "unlearned"
    assert "affected_shard" in res
    assert "new_model_fingerprint" in res

def test_sisa_training_single_class_error(tmp_path):
    """Test that SISA training with a single-class target raises a clear, descriptive ValueError."""
    data = pd.DataFrame({
        "case_id": [f"CASE_{i}" for i in range(50)],
        "drug_name": ["DrugX"] * 50,
        "event_term": ["Event1"] * 50,
        "case_count": [10] * 50,
        "prr": [2.5] * 50,
        "ror": [2.6] * 50,
        "chi2": [5.0] * 50,
        "label": [0] * 50  # all zeroes (unconfirmed reference controls)
    })
    trainer = SISATrainer(model_dir=tmp_path / "sisa_single_class", n_shards=2)
    with pytest.raises(ValueError, match="target 'label' contains only one class"):
        trainer.train(data, n_shards=2)
