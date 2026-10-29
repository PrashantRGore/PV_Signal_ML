"""
OMOP / OHDSI Ground Truth Reference Benchmark
=============================================
Zero-Cost, Open-Access Ground Truth Reference Set for Pharmacovigilance Machine Learning.

Legal Attribution & Citations:
- Ryan PB, Schuemie MJ, Gruber S, et al. Empirical performance evaluation of
  methods for risk identification in healthcare data: results from the Experiments
  of the Observational Medical Outcomes Partnership. Statistics in Medicine 2013; 32(30): 5601-5613.
  DOI: 10.1002/sim.5869.
- OHDSI OMOP HOI (Health Outcomes of Interest) Reference Set. Available under CC-BY 4.0 / Open Science.
- US National Library of Medicine (NLM) / openFDA: Public Domain.

This module provides verified positive controls (established adverse drug events)
and negative controls (unassociated drug-event pairs) for four canonical outcomes:
1. Acute Myocardial Infarction (AMI)
2. Upper Gastrointestinal Bleed (UGIB)
3. Acute Kidney Injury (AKI)
4. Drug-Induced Liver Injury (DILI)

Using this external reference standard eliminates target leakage and prevents circular
models that train on PRR/Chi² to predict an arbitrary PRR threshold.
"""

from typing import Dict, List, Optional, Set, Tuple
import pandas as pd
import numpy as np


class OMOPBenchmark:
    """
    Gold standard ground truth reference set from OMOP / Ryan et al. (2013).
    """

    # Health Outcome of Interest (HOI) Synonym Mappings (Free openFDA/MeSH/MedDRA PT strings)
    HOI_TERMS: Dict[str, Set[str]] = {
        "AMI": {
            "MYOCARDIAL INFARCTION",
            "ACUTE MYOCARDIAL INFARCTION",
            "CORONARY ARTERY OCCLUSION",
            "ACUTE CORONARY SYNDROME",
            "CARDIAC ARREST",
            "MYOCARDIAL ISCHAEMIA",
            "MYOCARDIAL INFARCTION POSTPROCEDURAL"
        },
        "UGIB": {
            "GASTROINTESTINAL HAEMORRHAGE",
            "GASTROINTESTINAL BLEEDING",
            "GASTRIC ULCER HAEMORRHAGE",
            "DUODENAL ULCER HAEMORRHAGE",
            "HAEMATEMESIS",
            "MELAENA",
            "PEPTIC ULCER HAEMORRHAGE",
            "UPPER GASTROINTESTINAL HAEMORRHAGE",
            "GASTRIC HAEMORRHAGE"
        },
        "AKI": {
            "ACUTE KIDNEY INJURY",
            "RENAL FAILURE ACUTE",
            "RENAL FAILURE",
            "RENAL IMPAIRMENT",
            "NEPHROTOXICITY",
            "TUBULOINTERSTITIAL NEPHRITIS",
            "BLOOD CREATININE INCREASED",
            "ANURIA",
            "OLIGURIA"
        },
        "DILI": {
            "HEPATOTOXICITY",
            "HEPATIC FAILURE",
            "ACUTE HEPATIC FAILURE",
            "LIVER INJURY",
            "DRUG-INDUCED LIVER INJURY",
            "JAUNDICE",
            "HEPATITIS TOXIC",
            "HEPATITIS ACUTE",
            "ALANINE AMINOTRANSFERASE INCREASED",
            "ASPARTATE AMINOTRANSFERASE INCREASED"
        }
    }

    # Reference Pairs: (Drug, Outcome) -> Label (1 = Positive Control, 0 = Negative Control)
    # Sourced directly from OMOP experiments (Ryan et al., 2013)
    REFERENCE_PAIRS: Dict[Tuple[str, str], int] = {
        # --- Acute Myocardial Infarction (AMI) ---
        ("ROFECOXIB", "AMI"): 1,
        ("CELECOXIB", "AMI"): 1,
        ("ROSIGLITAZONE", "AMI"): 1,
        ("DICLOFENAC", "AMI"): 1,
        ("INDOMETHACIN", "AMI"): 1,
        ("IBUPROFEN", "AMI"): 1,
        ("NAPROXEN", "AMI"): 1,
        ("MELOXICAM", "AMI"): 1,
        # AMI Negative Controls
        ("AMOXICILLIN", "AMI"): 0,
        ("METFORMIN", "AMI"): 0,
        ("SERTRALINE", "AMI"): 0,
        ("OMEPRAZOLE", "AMI"): 0,
        ("SIMVASTATIN", "AMI"): 0,
        ("LEVOTHYROXINE", "AMI"): 0,
        ("LORAZEPAM", "AMI"): 0,
        ("LISINOPRIL", "AMI"): 0,
        ("HYDROCHLOROTHIAZIDE", "AMI"): 0,
        ("AMLODIPINE", "AMI"): 0,
        ("WARFARIN", "AMI"): 0,
        ("CETIRIZINE", "AMI"): 0,

        # --- Upper Gastrointestinal Bleed (UGIB) ---
        ("IBUPROFEN", "UGIB"): 1,
        ("NAPROXEN", "UGIB"): 1,
        ("KETOROLAC", "UGIB"): 1,
        ("ASPIRIN", "UGIB"): 1,
        ("CELECOXIB", "UGIB"): 1,
        ("MELOXICAM", "UGIB"): 1,
        ("INDOMETHACIN", "UGIB"): 1,
        ("WARFARIN", "UGIB"): 1,
        ("CLOPIDOGREL", "UGIB"): 1,
        ("PIROXICAM", "UGIB"): 1,
        ("KETOPROFEN", "UGIB"): 1,
        # UGIB Negative Controls
        ("AMOXICILLIN", "UGIB"): 0,
        ("CITALOPRAM", "UGIB"): 0,
        ("METFORMIN", "UGIB"): 0,
        ("SIMVASTATIN", "UGIB"): 0,
        ("AMLODIPINE", "UGIB"): 0,
        ("LEVOTHYROXINE", "UGIB"): 0,
        ("ATORVASTATIN", "UGIB"): 0,
        ("HYDROCHLOROTHIAZIDE", "UGIB"): 0,
        ("ACETAMINOPHEN", "UGIB"): 0,
        ("DIPHENHYDRAMINE", "UGIB"): 0,

        # --- Acute Kidney Injury (AKI) ---
        ("GENTAMICIN", "AKI"): 1,
        ("CISPLATIN", "AKI"): 1,
        ("VANCOMYCIN", "AKI"): 1,
        ("AMPHOTERICIN B", "AKI"): 1,
        ("LISINOPRIL", "AKI"): 1,
        ("IBUPROFEN", "AKI"): 1,
        ("CAPTOPRIL", "AKI"): 1,
        ("LOSARTAN", "AKI"): 1,
        ("CYCLOSPORINE", "AKI"): 1,
        ("TACROLIMUS", "AKI"): 1,
        ("ACICLOVIR", "AKI"): 1,
        # AKI Negative Controls
        ("AMOXICILLIN", "AKI"): 0,
        ("OMEPRAZOLE", "AKI"): 0,
        ("CITALOPRAM", "AKI"): 0,
        ("METFORMIN", "AKI"): 0,
        ("ATORVASTATIN", "AKI"): 0,
        ("SERTRALINE", "AKI"): 0,
        ("CETIRIZINE", "AKI"): 0,
        ("LEVOTHYROXINE", "AKI"): 0,
        ("PREDNISONE", "AKI"): 0,

        # --- Drug-Induced Liver Injury (DILI) ---
        ("ACETAMINOPHEN", "DILI"): 1,
        ("AMOXICILLIN-CLAVULANATE", "DILI"): 1,
        ("ISONIAZID", "DILI"): 1,
        ("VALPROIC ACID", "DILI"): 1,
        ("METHOTREXATE", "DILI"): 1,
        ("KETOCONAZOLE", "DILI"): 1,
        ("ATORVASTATIN", "DILI"): 1,
        ("AZATHIOPRINE", "DILI"): 1,
        ("RIFAMPICIN", "DILI"): 1,
        ("PYRAZINAMIDE", "DILI"): 1,
        # DILI Negative Controls
        ("METFORMIN", "DILI"): 0,
        ("LISINOPRIL", "DILI"): 0,
        ("HYDROCHLOROTHIAZIDE", "DILI"): 0,
        ("OMEPRAZOLE", "DILI"): 0,
        ("LEVOTHYROXINE", "DILI"): 0,
        ("AMLODIPINE", "DILI"): 0,
        ("INSULIN", "DILI"): 0,
        ("ALBUTEROL", "DILI"): 0,
        ("LORAZEPAM", "DILI"): 0,
        ("DIPHENHYDRAMINE", "DILI"): 0
    }

    @classmethod
    def match_hoi(cls, event_term: str) -> Optional[str]:
        """Match adverse event string to canonical HOI category."""
        if not event_term or not isinstance(event_term, str):
            return None
        term_clean = event_term.strip().upper()
        for hoi, terms in cls.HOI_TERMS.items():
            if term_clean in terms:
                return hoi
            # Fallback substring match
            for t in terms:
                if t in term_clean or term_clean in t:
                    return hoi
        return None

    @classmethod
    def normalize_drug(cls, drug_name: str) -> str:
        """Standardize drug name for reference lookup."""
        if not drug_name or not isinstance(drug_name, str):
            return ""
        name = drug_name.strip().upper()
        # Clean common salt suffixes and formulation noise
        noise_tokens = [" SODIUM", " POTASSIUM", " HYDROCHLORIDE", " HCL", " TABLET", " CAPSULE", " INJECTION", " ORAL"]
        for token in noise_tokens:
            name = name.replace(token, "")
        return name.strip()

    @classmethod
    def get_ground_truth(cls, drug_name: str, event_term: str) -> Tuple[Optional[int], Optional[str], str]:
        """
        Lookup ground truth status for a drug-event pair.

        Returns:
            (label, hoi_category, control_status)
            - label: 1 (positive control), 0 (negative control), None (unconfirmed)
            - hoi_category: 'AMI', 'UGIB', 'AKI', 'DILI', or None
            - control_status: 'POSITIVE_CONTROL', 'NEGATIVE_CONTROL', or 'UNLABELED'
        """
        drug_norm = cls.normalize_drug(drug_name)
        hoi = cls.match_hoi(event_term)

        if hoi is None:
            return None, None, "UNLABELED"

        key = (drug_norm, hoi)
        if key in cls.REFERENCE_PAIRS:
            label = cls.REFERENCE_PAIRS[key]
            status = "POSITIVE_CONTROL" if label == 1 else "NEGATIVE_CONTROL"
            return label, hoi, status

        # Drug exact match not found; check partial drug name
        for (ref_drug, ref_hoi), ref_label in cls.REFERENCE_PAIRS.items():
            if ref_hoi == hoi:
                if ref_drug in drug_norm or drug_norm in ref_drug:
                    status = "POSITIVE_CONTROL" if ref_label == 1 else "NEGATIVE_CONTROL"
                    return ref_label, hoi, status

        return None, hoi, "UNLABELED"

    @classmethod
    def label_dataframe(
        cls,
        df: pd.DataFrame,
        drug_col: str = "drug_name",
        event_col: str = "event_term",
        benchmark_only: bool = False
    ) -> pd.DataFrame:
        """
        Enrich signal dataframe with OMOP ground truth labels.

        Args:
            df: DataFrame containing drug and event columns
            drug_col: Column name for drug
            event_col: Column name for adverse event
            benchmark_only: If True, keep only verified positive and negative controls

        Returns:
            Enriched DataFrame with columns:
            - ground_truth_label (float: 1.0, 0.0, or np.nan)
            - hoi_category (str or None)
            - control_status ('POSITIVE_CONTROL', 'NEGATIVE_CONTROL', 'UNLABELED')
        """
        res = df.copy()
        labels: List[Optional[float]] = []
        hois: List[Optional[str]] = []
        statuses: List[str] = []

        for _, row in res.iterrows():
            drug = str(row.get(drug_col, ""))
            event = str(row.get(event_col, ""))
            lbl, hoi, status = cls.get_ground_truth(drug, event)
            labels.append(float(lbl) if lbl is not None else np.nan)
            hois.append(hoi)
            statuses.append(status)

        res["ground_truth_label"] = labels
        res["hoi_category"] = hois
        res["control_status"] = statuses

        if benchmark_only:
            res = res.dropna(subset=["ground_truth_label"]).reset_index(drop=True)
            res["ground_truth_label"] = res["ground_truth_label"].astype(int)

        return res

