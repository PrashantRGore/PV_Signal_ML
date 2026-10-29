import pytest
import pandas as pd
import numpy as np
from src.stats_engine.disproportionality import DisproportionalityAnalysis


def test_disproportionality_formulas():
    """Verify PRR, ROR and Chi2 calculations against known analytical results."""
    # Build synthetic ICSR dataset
    # Total reports = 100
    # Drug A with Event E = 10 (a=10)
    # Drug A with Event Other = 10 (b=10) -> Drug A total = 20
    # Other Drug with Event E = 10 (c=10)
    # Other Drug with Event Other = 70 (d=70) -> Other Drug total = 80
    # Total = 100
    
    records = []
    # a=10
    for i in range(10):
        records.append({"case_id": f"c_{len(records)}", "drug_name": "DrugA", "event_term": "EventE"})
    # b=10
    for i in range(10):
        records.append({"case_id": f"c_{len(records)}", "drug_name": "DrugA", "event_term": "EventOther"})
    # c=10
    for i in range(10):
        records.append({"case_id": f"c_{len(records)}", "drug_name": "DrugB", "event_term": "EventE"})
    # d=70
    for i in range(70):
        records.append({"case_id": f"c_{len(records)}", "drug_name": "DrugB", "event_term": "EventOther"})
        
    df = pd.DataFrame(records)
    
    engine = DisproportionalityAnalysis(min_cases=3)
    signals = engine.compute_signals(df)
    
    # Filter for DrugA + EventE
    signal_a = signals[(signals["drug_name"] == "DrugA") & (signals["event_term"] == "EventE")].iloc[0]
    
    # Analytical expected:
    # a = 10, b = 10, c = 10, d = 70
    # PRR = [a / (a + b)] / [c / (c + d)] = (10 / 20) / (10 / 80) = 0.5 / 0.125 = 4.0
    # ROR = (a * d) / (b * c) = (10 * 70) / (10 * 10) = 700 / 100 = 7.0
    assert pytest.approx(signal_a["prr"], 0.01) == 4.0
    assert pytest.approx(signal_a["ror"], 0.01) == 7.0
    assert signal_a["prr"] != signal_a["ror"], "PRR and ROR must have distinct formulas"
    
    # Chi2 with Yates correction:
    # N = 100, |ad - bc| = |700 - 100| = 600, diff = 600 - 50 = 550
    # denom = (20) * (80) * (20) * (80) = 2,560,000
    # chi2 = 100 * (550^2) / 2,560,000 = 100 * 302,500 / 2,560,000 = 11.816
    assert pytest.approx(signal_a["chi2"], 0.1) == 11.816

