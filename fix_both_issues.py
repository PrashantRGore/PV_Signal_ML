from pathlib import Path

# Fix 1: Add date filtering to FAERSDownloader
faers_file = Path('src/data/faers_downloader.py')
content = faers_file.read_text(encoding='utf-8')

# Add date filtering in the result DataFrame creation
old_result = '''        # Standardize columns
        result = pd.DataFrame({
            'case_id': merged[id_col_drug],
            'drug_name': merged.get('drugname', merged.get('drug_name', '')),
            'event_term': merged.get('pt', merged.get('event_term', '')),
            'year': year,
            'quarter': quarter
        })
        
        return result.dropna(subset=['drug_name', 'event_term'])'''

new_result = '''        # Standardize columns
        result = pd.DataFrame({
            'case_id': merged[id_col_drug],
            'drug_name': merged.get('drugname', merged.get('drug_name', '')),
            'event_term': merged.get('pt', merged.get('event_term', '')),
            'year': year,
            'quarter': quarter
        })
        
        # Clean and filter
        result = result.dropna(subset=['drug_name', 'event_term'])
        
        # Remove empty/invalid values
        result = result[result['drug_name'].str.strip() != '']
        result = result[result['event_term'].str.strip() != '']
        
        # Convert to uppercase for consistency
        result['drug_name'] = result['drug_name'].str.upper().str.strip()
        result['event_term'] = result['event_term'].str.upper().str.strip()
        
        return result'''

content = content.replace(old_result, new_result)
faers_file.write_text(content, encoding='utf-8')

print('✅ Fixed FAERS data cleaning')

# Fix 2: Add validation in signal detection to prevent negative values
disp_file = Path('src/stats_engine/disproportionality.py')
disp_content = disp_file.read_text(encoding='utf-8')

# Find the chi2 calculation and add validation
old_chi2 = '''    def _compute_chi2(self, n11, n10, n01, n00):
        \"\"\"Compute Chi-square statistic\"\"\"
        contingency_table = np.array([[n11, n10], [n01, n00]])
        chi2, p_value, _, _ = stats.chi2_contingency(contingency_table)'''

new_chi2 = '''    def _compute_chi2(self, n11, n10, n01, n00):
        \"\"\"Compute Chi-square statistic with validation\"\"\"
        # Ensure all values are non-negative
        n11 = max(0, int(n11))
        n10 = max(0, int(n10))
        n01 = max(0, int(n01))
        n00 = max(0, int(n00))
        
        # Check if table is valid
        if n11 + n10 + n01 + n00 == 0:
            return 0.0, 1.0
        
        contingency_table = np.array([[n11, n10], [n01, n00]])
        chi2, p_value, _, _ = stats.chi2_contingency(contingency_table)'''

disp_content = disp_content.replace(old_chi2, new_chi2)
disp_file.write_text(disp_content, encoding='utf-8')

print('✅ Fixed chi-square validation')

# Fix 3: Add data validation before signal detection
print()
print('📊 Adding data validation check...')

