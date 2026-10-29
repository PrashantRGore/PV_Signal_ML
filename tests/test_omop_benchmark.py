import pytest
import pandas as pd
import numpy as np
from src.data.standards.omop_benchmark import OMOPBenchmark
from src.data.icsr_to_signals import ICSRToSignalConverter


def test_omop_benchmark_lookup():
    """Verify positive and negative control lookup in OMOP benchmark."""
    # Positive control: Celecoxib - AMI
    label, hoi, status = OMOPBenchmark.get_ground_truth("CELECOXIB", "MYOCARDIAL INFARCTION")
    assert label == 1
    assert hoi == "AMI"
    assert status == "POSITIVE_CONTROL"

    # Negative control: Amoxicillin - AMI
    label_neg, hoi_neg, status_neg = OMOPBenchmark.get_ground_truth("AMOXICILLIN", "MYOCARDIAL INFARCTION")
    assert label_neg == 0
    assert hoi_neg == "AMI"
    assert status_neg == "NEGATIVE_CONTROL"

    # Positive control: Acetaminophen - Acute Hepatic Failure
    label_dili, hoi_dili, status_dili = OMOPBenchmark.get_ground_truth("ACETAMINOPHEN", "HEPATIC FAILURE")
    assert label_dili == 1
    assert hoi_dili == "DILI"
    assert status_dili == "POSITIVE_CONTROL"

    # Unlabeled pair
    label_unk, hoi_unk, status_unk = OMOPBenchmark.get_ground_truth("UNKNOWN_DRUG", "UNKNOWN_EVENT")
    assert label_unk is None
    assert status_unk == "UNLABELED"


def test_icsr_to_signals_converter():
    """Verify ICSR conversion without mock data and with OMOP labeling."""
    records = []
    # 5 reports of CELECOXIB with MYOCARDIAL INFARCTION
    for i in range(5):
        records.append({"case_id": f"c_{len(records)}", "drug_name": "CELECOXIB", "event_term": "MYOCARDIAL INFARCTION"})
    # 5 reports of CELECOXIB with HEADACHE
    for i in range(5):
        records.append({"case_id": f"c_{len(records)}", "drug_name": "CELECOXIB", "event_term": "HEADACHE"})
    # 5 reports of AMOXICILLIN with MYOCARDIAL INFARCTION
    for i in range(5):
        records.append({"case_id": f"c_{len(records)}", "drug_name": "AMOXICILLIN", "event_term": "MYOCARDIAL INFARCTION"})
    # 35 reports of AMOXICILLIN with HEADACHE
    for i in range(35):
        records.append({"case_id": f"c_{len(records)}", "drug_name": "AMOXICILLIN", "event_term": "HEADACHE"})

    df = pd.DataFrame(records)
    converter = ICSRToSignalConverter(min_cases=3)
    signals = converter.convert(df)

    assert len(signals) > 0
    assert "drug_name" in signals.columns
    assert "event_term" in signals.columns
    assert "prr" in signals.columns
    assert "ror" in signals.columns
    assert "chi2" in signals.columns
    assert "ground_truth_label" in signals.columns
    assert "is_evans_signal" in signals.columns

    # Verify Celecoxib - AMI got positive control label (1)
    celecoxib_ami = signals[(signals["drug_name"] == "CELECOXIB") & (signals["event_term"] == "MYOCARDIAL INFARCTION")]
    assert len(celecoxib_ami) == 1
    assert celecoxib_ami.iloc[0]["ground_truth_label"] == 1.0

    # Verify Amoxicillin - AMI got negative control label (0)
    amoxicillin_ami = signals[(signals["drug_name"] == "AMOXICILLIN") & (signals["event_term"] == "MYOCARDIAL INFARCTION")]
    assert len(amoxicillin_ami) == 1
    assert amoxicillin_ami.iloc[0]["ground_truth_label"] == 0.0

