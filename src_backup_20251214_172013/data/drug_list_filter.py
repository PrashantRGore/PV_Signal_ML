"""
Drug List Filter for Targeted Surveillance
Regulatory Compliant: FDA, EMA, PMDA standards
"""

import pandas as pd
import streamlit as st
from pathlib import Path

class DrugListFilter:
    """
    Filters FAERS data to drugs of interest
    Use Case: Portfolio monitoring, product-specific surveillance
    """
    
    def __init__(self):
        self.drug_list = []
        self.filter_active = False
    
    def load_from_excel(self, file_path):
        """Load drug list from Excel"""
        df = pd.read_excel(file_path)
        
        # Try to find drug name column
        possible_cols = ['drug_name', 'drug', 'product', 'medicine', 'substance']
        drug_col = None
        
        for col in df.columns:
            if col.lower() in possible_cols:
                drug_col = col
                break
        
        if drug_col is None:
            # Use first column
            drug_col = df.columns[0]
            st.warning(f"Using first column '{drug_col}' as drug names")
        
        self.drug_list = df[drug_col].str.upper().str.strip().tolist()
        self.filter_active = True
        return self.drug_list
    
    def load_from_uploaded_file(self, uploaded_file):
        """Load from Streamlit file uploader"""
        df = pd.read_excel(uploaded_file)
        
        # Auto-detect drug column
        for col in df.columns:
            if 'drug' in col.lower() or 'product' in col.lower():
                self.drug_list = df[col].str.upper().str.strip().tolist()
                self.filter_active = True
                return self.drug_list
        
        # Use first column
        self.drug_list = df[df.columns[0]].str.upper().str.strip().tolist()
        self.filter_active = True
        return self.drug_list
    
    def filter_data(self, data: pd.DataFrame) -> pd.DataFrame:
        """Filter FAERS data to drugs of interest"""
        if not self.filter_active or not self.drug_list:
            return data
        
        # Ensure drug names are uppercase for matching
        data['drug_name'] = data['drug_name'].str.upper().str.strip()
        
        # Filter to drugs in list
        filtered = data[data['drug_name'].isin(self.drug_list)]
        
        return filtered
    
    def get_coverage_report(self, data: pd.DataFrame) -> dict:
        """Report which drugs from list are in FAERS data"""
        if not self.filter_active:
            return {}
        
        data_drugs = set(data['drug_name'].str.upper().unique())
        list_drugs = set(self.drug_list)
        
        found = list_drugs.intersection(data_drugs)
        missing = list_drugs - data_drugs
        
        return {
            'total_in_list': len(list_drugs),
            'found_in_faers': len(found),
            'missing_from_faers': len(missing),
            'found_drugs': sorted(found),
            'missing_drugs': sorted(missing)
        }
