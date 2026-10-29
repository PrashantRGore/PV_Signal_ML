import pytest
from src.rag.sar_generator import SARGenerator


def test_sar_generator_deterministic_fallback():
    """Verify that SARGenerator produces a complete, compliant SAR even when Ollama is offline."""
    # Point to an unreachable port to test offline fallback
    generator = SARGenerator(model_name="llama3.2", base_url="http://localhost:9999")

    signal_data = {
        "drug_name": "IBUPROFEN",
        "event_name": "GASTROINTESTINAL HAEMORRHAGE",
        "prr": 3.45,
        "ror": 4.12,
        "chi_square": 24.8,
        "count": 42
    }

    papers = [
        {
            "authors": ["Smith JA", "Doe R"],
            "year": "2023",
            "title": "NSAID-induced gastrointestinal bleeding risk",
            "journal": "Pharmacoepidemiol Drug Saf",
            "pmid": "12345678"
        }
    ]

    result = generator.generate_sar(signal_data=signal_data, literature_papers=papers)

    assert result is not None
    assert "report" in result
    report = result["report"]

    # Verify standard regulatory headings
    assert "# SAFETY ASSESSMENT REPORT (SAR)" in report
    assert "## 1. Executive Summary" in report
    assert "## 3. Disproportionality & Statistical Evidence" in report
    assert "## 7. Pharmacovigilance Recommendations" in report
    assert "## 9. References" in report

    # Verify quantitative numbers are accurately reflected
    assert "3.45" in report
    assert "4.12" in report
    assert "24.80" in report or "24.8" in report
    assert "42" in report

    # Verify literature citation was included
    assert "Smith JA" in report or "12345678" in report

