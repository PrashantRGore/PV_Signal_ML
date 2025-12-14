from pathlib import Path

print('FORCING CONTINGENCY FIX...')
disp_file = Path('src/stats_engine/disproportionality.py')
content = disp_file.read_text(encoding='utf-8')

print('Current groupby line:')
import re
match = re.search(r"drug_events = data\.groupby.*", content)
if match:
    print(f'  {match.group(0)}')

# FORCE replacement
old_line = "data.groupby(['drug_name', 'event_term']).size().reset_index(name='n11')"
new_line = "data.groupby(['drug_name', 'event_term'])['case_id'].nunique().reset_index(name='n11')"

if old_line in content:
    content = content.replace(old_line, new_line)
    print('✅ Replaced .size() with .nunique()')
else:
    print('⚠️  Old pattern not found, checking if already fixed...')
    if 'nunique()' in content:
        print('✅ Already uses nunique()')
    else:
        print('❌ PROBLEM: Cannot find pattern')

# Replace all len(unique()) with nunique()
content = content.replace(
    "len(data['case_id'].unique())",
    "data['case_id'].nunique()"
)
content = content.replace(
    "len(data[data['drug_name'] == drug]['case_id'].unique())",
    "data[data['drug_name'] == drug]['case_id'].nunique()"
)
content = content.replace(
    "len(data[data['event_term'] == event]['case_id'].unique())",
    "data[data['event_term'] == event]['case_id'].nunique()"
)

disp_file.write_text(content, encoding='utf-8')
print('✅ File written')

# VERIFY
content_check = disp_file.read_text(encoding='utf-8')
if '.nunique()' in content_check and 'drug_events = data.groupby' in content_check:
    line = [l for l in content_check.split('\n') if 'drug_events = data.groupby' in l][0]
    print(f'✅ VERIFIED: {line.strip()}')
else:
    print('❌ VERIFICATION FAILED')
