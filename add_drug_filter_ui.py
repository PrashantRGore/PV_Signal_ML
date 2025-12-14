from pathlib import Path

print('=' * 70)
print('IMPLEMENTING DRUG PORTFOLIO FILTER (Targeted Surveillance)')
print('=' * 70)
print()

# Step 1: Ensure drug filter component exists
print('1️⃣ Creating drug filter component...')
filter_code = '''"""
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
'''

filter_file = Path('src/utils/drug_filter.py')
filter_file.parent.mkdir(parents=True, exist_ok=True)
filter_file.write_text(filter_code, encoding='utf-8')
print(f'✅ Created: {filter_file}')

# Step 2: Add filter UI to app
print()
print('2️⃣ Adding drug filter UI to app...')

app_file = Path('app_enhanced.py')
content = app_file.read_text(encoding='utf-8')

# Add after signal detection results
filter_ui_code = '''
    # Drug Portfolio Filter
    with st.sidebar.expander("🎯 Filter by Drug Portfolio"):
        st.write("**Targeted Surveillance**")
        st.write("Upload Excel with your portfolio drugs to focus analysis")
        
        uploaded_drug_list = st.file_uploader(
            "Upload drug list (Excel)",
            type=['xlsx', 'xls'],
            help="Excel file with drug names (any column)",
            key="drug_filter_upload"
        )
        
        if uploaded_drug_list:
            from src.utils.drug_filter import DrugFilter
            
            if 'drug_filter' not in st.session_state:
                st.session_state.drug_filter = DrugFilter()
            
            drug_count = st.session_state.drug_filter.load_from_excel(uploaded_drug_list)
            
            if drug_count > 0:
                st.success(f"✅ Loaded {drug_count} portfolio drugs")
                
                with st.expander("View Portfolio"):
                    drugs = st.session_state.drug_filter.drug_list
                    st.write(drugs[:20])
                    if len(drugs) > 20:
                        st.write(f"... and {len(drugs)-20} more")
    
    # Apply drug filter to signals if active
    if st.session_state.signals is not None and 'drug_filter' in st.session_state:
        if st.session_state.drug_filter.active:
            original_count = len(st.session_state.signals)
            filtered_signals = st.session_state.drug_filter.filter_signals(st.session_state.signals)
            
            st.info(f"🎯 **Portfolio Filter Active:** {len(filtered_signals):,} signals from {original_count:,} total")
            
            # Coverage report
            coverage = st.session_state.drug_filter.get_coverage_report(st.session_state.signals)
            
            col1, col2, col3 = st.columns(3)
            col1.metric("Portfolio Drugs", coverage['total_portfolio'])
            col2.metric("Drugs with Signals", coverage['drugs_with_signals'])
            col3.metric("Drugs without Signals", coverage['drugs_without_signals'])
            
            if coverage['drugs_without_signals'] > 0:
                with st.expander(f"View {coverage['drugs_without_signals']} drugs without signals"):
                    st.write(coverage['missing_list'])
            
            # Use filtered signals for display
            st.session_state.signals = filtered_signals
'''

# Find where to insert (after signal detection but before Top Signals display)
marker = 'st.success(f"✅ Detected {len(signals)} signals")'

if marker in content:
    content = content.replace(marker, marker + filter_ui_code)
    app_file.write_text(content, encoding='utf-8')
    print('✅ Added filter UI to app')
else:
    print('⚠️  Marker not found, trying alternative location...')
    # Alternative: add before "Top Signals" section
    marker2 = 'st.subheader("Top Signals")'
    if marker2 in content:
        content = content.replace(marker2, filter_ui_code + '\n    ' + marker2)
        app_file.write_text(content, encoding='utf-8')
        print('✅ Added filter UI (alternative location)')

print()
print('=' * 70)
print('✅ DRUG PORTFOLIO FILTER READY!')
print('=' * 70)
print()
print('📋 How to use:')
print('  1. Run signal detection as usual')
print('  2. In sidebar, expand "🎯 Filter by Drug Portfolio"')
print('  3. Upload Excel file with your drugs')
print('  4. Signals instantly filtered to your portfolio!')
print()
print('📊 Excel format (any column name works):')
print('  | drug_name     | or | Product      | or | Medicine     |')
print('  |---------------|    |--------------|    |--------------|')
print('  | ASPIRIN       |    | IBUPROFEN    |    | METFORMIN    |')
print('  | WARFARIN      |    | ATORVASTATIN |    | LISINOPRIL   |')
print()
print('🚀 Restart app to see the filter!')
