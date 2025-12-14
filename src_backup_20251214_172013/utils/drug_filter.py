"""
Drug Portfolio Filter for Targeted Pharmacovigilance
"""

import pandas as pd
import streamlit as st

class DrugFilter:
    def __init__(self):
        self.drug_list = []
        self.active = False
    
    def load_from_excel(self, uploaded_file):
        """Load drug list from uploaded Excel file"""
        try:
            df = pd.read_excel(uploaded_file)
            
            # Auto-detect drug column
            drug_col = None
            possible_names = ['drug', 'product', 'medicine', 'substance', 'name']
            
            for col in df.columns:
                if any(name in col.lower() for name in possible_names):
                    drug_col = col
                    break
            
            if drug_col is None:
                drug_col = df.columns[0]
                st.info(f"Using first column '{drug_col}' as drug names")
            
            # Extract and clean drug names
            self.drug_list = df[drug_col].astype(str).str.upper().str.strip().tolist()
            self.drug_list = [d for d in self.drug_list if d and d != 'NAN']
            self.active = True
            
            return len(self.drug_list)
            
        except Exception as e:
            st.error(f"Error reading Excel: {e}")
            return 0
    
    def filter_signals(self, signals_df):
        """Filter signals to only show drugs in portfolio"""
        if not self.active or not self.drug_list:
            return signals_df
        
        # Ensure drug names are uppercase for matching
        signals_df['drug_name_upper'] = signals_df['drug_name'].str.upper().str.strip()
        
        # Filter to portfolio drugs
        filtered = signals_df[signals_df['drug_name_upper'].isin(self.drug_list)].copy()
        filtered = filtered.drop('drug_name_upper', axis=1)
        
        return filtered
    
    def get_coverage_report(self, signals_df):
        """Report which portfolio drugs have signals"""
        if not self.active:
            return {}
        
        signal_drugs = set(signals_df['drug_name'].str.upper().unique())
        portfolio_set = set(self.drug_list)
        
        found = portfolio_set.intersection(signal_drugs)
        missing = portfolio_set - signal_drugs
        
        return {
            'total_portfolio': len(portfolio_set),
            'drugs_with_signals': len(found),
            'drugs_without_signals': len(missing),
            'found_list': sorted(found),
            'missing_list': sorted(missing)
        }
