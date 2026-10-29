import pytest
from src.regulatory.causality_assessor import (
    CausalityAssessor,
    WHOUMCCausality,
    NaranjoCausality,
    CausalityCategory
)


def test_naranjo_causality_dynamic_scoring():
    """Verify that Naranjo causality evaluates real case fields dynamically."""
    # Case with clear positive rechallenge, dechallenge, known adverse event, lab confirmation, and no other drugs
    strong_case = {
        "ground_truth_label": 1,
        "time_to_onset_days": 5,
        "dechallenge_positive": 1,
        "rechallenge_positive": 1,
        "concomitant_count": 0,
        "concomitant_drugs": [],
        "alt_elevated": 1
    }

    result = NaranjoCausality.assess(strong_case)
    # Score should be >= 5 (Probable or Definite)
    assert result["score"] >= 5
    assert result["category"] in ("Probable", "Definite")
    assert result["item_scores"][1] == 2  # Q2: Appeared after drug
    assert result["item_scores"][2] == 1  # Q3: Dechallenge positive
    assert result["item_scores"][3] == 2  # Q4: Rechallenge positive
    assert result["item_scores"][9] == 1  # Q10: Lab confirmed


def test_who_umc_causality_categories():
    """Verify WHO-UMC causality assignment based on clinical criteria."""
    # Definite / Certain case
    certain_case = {
        "time_to_onset_days": 3,
        "dechallenge_positive": 1,
        "rechallenge_positive": 1,
        "concomitant_count": 0
    }
    assessment = WHOUMCCausality.assess(certain_case)
    assert assessment["category"] == CausalityCategory.CERTAIN.value

    # Temporally impossible case (event before drug)
    unlikely_case = {
        "time_to_onset_days": -10,
        "dechallenge_positive": 0,
        "rechallenge_positive": 0
    }
    assessment_unlikely = WHOUMCCausality.assess(unlikely_case)
    assert assessment_unlikely["category"] == CausalityCategory.UNLIKELY.value


def test_causality_assessor_consensus():
    """Verify combined consensus assessment."""
    assessor = CausalityAssessor()
    signal = {
        "ground_truth_label": 1,
        "time_to_onset_days": 4,
        "dechallenge_positive": 1,
        "concomitant_count": 1
    }
    consensus = assessor.assess_signal(signal, method="both")

    assert "who_umc" in consensus
    assert "naranjo" in consensus
    assert "consensus" in consensus
    assert consensus["consensus"]["level"] in ("High Causality", "Moderate Causality")

