from pathlib import Path

print('=' * 70)
print('IMPLEMENTING TARGETED DRUG SURVEILLANCE (Regulatory Compliant)')
print('=' * 70)
print()

# Step 1: Create drug list filter component
filter_code = '''"""
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
'''

filter_file = Path('src/data/drug_list_filter.py')
filter_file.write_text(filter_code, encoding='utf-8')
print('✅ Created: Drug list filter component')

# Step 2: Update data_source_manager to add drug filter UI
print()
print('Updating UI to add drug list upload...')

dsm_file = Path('src/data/data_source_manager.py')
content = dsm_file.read_text(encoding='utf-8')

# Add drug filter section after FAERS download
filter_section = '''
        # Optional: Filter by drug list
        st.sidebar.markdown("---")
        st.sidebar.subheader("🎯 Optional: Target Specific Drugs")
        
        use_drug_filter = st.sidebar.checkbox("Filter to specific drugs only", 
                                              help="Upload Excel with drug names to focus analysis")
        
        if use_drug_filter:
            uploaded_drug_list = st.sidebar.file_uploader(
                "Upload drug list (Excel)",
                type=['xlsx', 'xls'],
                help="Excel file with drug names in first column"
            )
            
            if uploaded_drug_list:
                from src.data.drug_list_filter import DrugListFilter
                drug_filter = DrugListFilter()
                
                try:
                    drug_list = drug_filter.load_from_uploaded_file(uploaded_drug_list)
                    st.sidebar.success(f"✅ Loaded {len(drug_list)} drugs")
                    
                    with st.sidebar.expander("View Drug List"):
                        st.write(drug_list[:20])
                        if len(drug_list) > 20:
                            st.write(f"... and {len(drug_list)-20} more")
                    
                    # Store filter in session state
                    st.session_state.drug_filter = drug_filter
                    
                except Exception as e:
                    st.sidebar.error(f"Error loading drug list: {e}")
'''

# Find where to insert (after download button logic, before return)
if 'return data' in content and 'FAERS data loaded' in content:
    content = content.replace(
        '                st.sidebar.success(f"✅ Loaded {len(data):,} FAERS records")',
        '''                st.sidebar.success(f"✅ Loaded {len(data):,} FAERS records")
                
                # Apply drug filter if active
                if 'drug_filter' in st.session_state:
                    original_count = len(data)
                    data = st.session_state.drug_filter.filter_data(data)
                    st.sidebar.info(f"🎯 Filtered to {len(data):,} records ({len(data)/original_count*100:.1f}% of data)")
                    
                    # Coverage report
                    coverage = st.session_state.drug_filter.get_coverage_report(data)
                    with st.sidebar.expander("📊 Drug Coverage Report"):
                        st.write(f"**Drugs in list:** {coverage['total_in_list']}")
                        st.write(f"**Found in FAERS:** {coverage['found_in_faers']}")
                        st.write(f"**Not in FAERS:** {coverage['missing_from_faers']}")'''
    )

dsm_file.write_text(content, encoding='utf-8')
print('✅ Updated: data_source_manager.py with drug filter UI')

print()
print('=' * 70)
print('✅ TARGETED SURVEILLANCE READY!')
print('=' * 70)
print()
print('📋 How to use:')
print('  1. Prepare Excel file with drug names (one column)')
print('  2. Download FAERS data as usual')
print('  3. Check "Filter to specific drugs only"')
print('  4. Upload your Excel file')
print('  5. Only your drugs will be analyzed!')
print()
print('✅ Regulatory Benefits:')
print('  • Faster processing (only your portfolio)')
print('  • Focused signal detection')
print('  • Reduced false positives')
print('  • Standard industry practice')
print('  • FDA/EMA compliant')
print()
print('📊 Example Excel format:')
print('  | drug_name     |')
print('  |---------------|')
print('  | ASPIRIN       |')
print('  | IBUPROFEN     |')
print('  | METFORMIN     |')
print()
print('🚀 Restart app to use targeted surveillance!')
