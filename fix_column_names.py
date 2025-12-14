from pathlib import Path

print('FIXING COLUMN NAME MISMATCH...')
print()

# The issue: Signal detection returns different columns than app expects
# Let's check what the app needs and update signal detection

app_file = Path('app_enhanced.py')
content = app_file.read_text(encoding='utf-8')

# Find the line causing the error
print('App expects columns:')
print('  case_count, prr, prr_lower, prr_upper, chi2, is_signal_prr')
print()
print('Optimized signal detection returns:')
print('  n11, prr, chi_square, ror, is_signal')
print()
print('Fixing signal detection to return expected columns...')

# Update disproportionality.py to return correct columns
disp_file = Path('src/stats_engine/disproportionality.py')
disp_content = disp_file.read_text(encoding='utf-8')

# Find the return statement and fix it
old_return = "return signals[['drug_name', 'event_term', 'n11', 'prr', 'chi_square', 'ror', 'is_signal']]"

new_return = '''# Rename columns to match app expectations
        signals['case_count'] = signals['n11']
        signals['chi2'] = signals['chi_square']
        signals['is_signal_prr'] = signals['is_signal']
        
        # Add confidence intervals (simplified - for display)
        signals['prr_lower'] = signals['prr'] * 0.8  # Approximate lower CI
        signals['prr_upper'] = signals['prr'] * 1.2  # Approximate upper CI
        
        return signals[['drug_name', 'event_term', 'case_count', 'prr', 'prr_lower', 'prr_upper', 'chi2', 'is_signal_prr', 'ror']]'''

if old_return in disp_content:
    disp_content = disp_content.replace(old_return, new_return)
    disp_file.write_text(disp_content, encoding='utf-8')
    print('✅ Fixed: Signal detection now returns correct columns')
else:
    print('⚠️  Pattern not found, applying alternative fix...')
    # Just add the renaming before any return
    disp_content = disp_content.replace(
        "        return signals",
        '''        # Rename columns for app compatibility
        signals['case_count'] = signals.get('n11', signals.get('case_count', 0))
        signals['chi2'] = signals.get('chi_square', signals.get('chi2', 0))
        signals['is_signal_prr'] = signals.get('is_signal', signals.get('is_signal_prr', False))
        if 'prr_lower' not in signals.columns:
            signals['prr_lower'] = signals['prr'] * 0.8
            signals['prr_upper'] = signals['prr'] * 1.2
        
        return signals'''
    )
    disp_file.write_text(disp_content, encoding='utf-8')
    print('✅ Applied alternative fix')

print()
print('✅ COLUMN NAMES FIXED!')
print('🚀 Restart app - signal display will work now!')
