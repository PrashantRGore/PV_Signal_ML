from pathlib import Path

disp_file = Path('src/stats_engine/disproportionality.py')
content = disp_file.read_text(encoding='utf-8')

# Show the current calculation method
print('=== Current contingency table calculation ===')
import re
match = re.search(r'n11 = row\[.n11.\].*?n00 =.*', content, re.DOTALL)
if match:
    print(match.group(0)[:500])
    print()

# Fix: Add robust contingency table calculation
old_calc = '''              n11 = row['n11']  # Drug + Event'''

# Check if we need to see more context first
lines = content.split('\n')
for i, line in enumerate(lines):
    if 'n11 = row' in line and 'n11' in line:
        print(f'Found at line {i}:')
        for j in range(max(0, i-2), min(len(lines), i+15)):
            print(f'{j}: {lines[j]}')
        break

print()
print('Please share the output above so I can create the exact fix!')
