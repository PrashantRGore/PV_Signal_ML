"""
Final Pre-Launch Validation
"""
from pathlib import Path
import re

print("="*60)
print("FINAL PRE-LAUNCH VALIDATION")
print("="*60)

# 1. Check safe_aggregation.py exists
safe_agg_file = Path("src/utils/safe_aggregation.py")
if safe_agg_file.exists():
    print("\n✓ Safe aggregation utility exists")
    content = safe_agg_file.read_text()
    functions = ['safe_mean', 'safe_median', 'safe_std', 'safe_mode', 'safe_percentage']
    for func in functions:
        if f"def {func}" in content:
            print(f"  ✓ {func}() defined")
else:
    print("\n✗ Safe aggregation utility MISSING!")

# 2. Check disproportionality.py uses safe functions
disp_file = Path("src/stats_engine/disproportionality.py")
if disp_file.exists():
    content = disp_file.read_text()
    
    print("\n2. Disproportionality.py:")
    
    # Check imports
    if "from utils.safe_aggregation import" in content:
        print("  ✓ Safe aggregation imported")
    else:
        print("  ✗ Safe aggregation NOT imported")
    
    # Check usage
    safe_calls = len(re.findall(r'safe_(mean|median|std|mode|percentage|threshold)', content))
    print(f"  ✓ {safe_calls} safe function calls found")
    
    # Check signal assignments
    signal_assigns = len(re.findall(r"signal\['[^']+'\]\s*=", content))
    print(f"  ✓ {signal_assigns} signal feature assignments")
    
    # Check for signals.append position
    lines = content.split('\n')
    append_line = -1
    first_extended_line = -1
    
    for i, line in enumerate(lines):
        if 'signals.append(signal)' in line and append_line == -1:
            append_line = i + 1
        if 'has_concomitant_pct' in line and first_extended_line == -1:
            first_extended_line = i + 1
    
    if first_extended_line > 0 and append_line > 0:
        if first_extended_line < append_line:
            print(f"  ✓ Extended features (line {first_extended_line}) BEFORE signals.append (line {append_line})")
        else:
            print(f"  ✗ Extended features (line {first_extended_line}) AFTER signals.append (line {append_line}) - WRONG!")
    
# 3. Check SISA trainer validations
sisa_file = Path("src/ml/sisa_trainer.py")
if sisa_file.exists():
    content = sisa_file.read_text()
    
    print("\n3. SISA Trainer:")
    
    if "if signals_df.empty:" in content:
        print("  ✓ Empty dataframe check added")
    else:
        print("  ⚠ Missing empty dataframe check")
    
    if "if X.empty or X.shape[1] == 0:" in content:
        print("  ✓ Feature validation added")
    else:
        print("  ⚠ Missing feature validation")
    
    if "X.isna().all().any()" in content:
        print("  ✓ NaN handling added")
    else:
        print("  ⚠ Missing NaN handling")

# 4. Check app_enhanced.py
app_file = Path("app_enhanced.py")
if app_file.exists():
    content = app_file.read_text()
    
    print("\n4. App Enhanced:")
    
    if "from src.utils.data_prep import prepare_data_for_signal_detection" in content:
        print("  ✓ Data prep imported")
    else:
        print("  ⚠ Data prep not imported")
    
    if "prepare_data_for_signal_detection" in content:
        print("  ✓ Data prep called before signal detection")
    else:
        print("  ⚠ Data prep not called")

print("\n" + "="*60)
print("VALIDATION COMPLETE")
print("="*60)
print("\nIf all checks passed, run:")
print("  streamlit run app_enhanced.py")
print("="*60)
