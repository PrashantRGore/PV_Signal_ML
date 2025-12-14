import re

# Read the app
with open('app_enhanced.py', 'r', encoding='utf-8') as f:
    content = f.read()

print('Searching for all SISATrainer initialization issues...\n')

# Fix ALL variations
fixes_made = 0

# Pattern 1: SISATrainer(num_shards=...)
pattern1 = r'SISATrainer\s*\(\s*num_shards\s*=\s*\w+\s*\)'
matches1 = re.findall(pattern1, content)
if matches1:
    print(f'Found {len(matches1)} instances of SISATrainer(num_shards=...)')
    for match in matches1:
        print(f'  - {match}')
    content = re.sub(pattern1, 'SISATrainer(model_dir)', content)
    fixes_made += len(matches1)

# Pattern 2: Make sure model_dir is defined before use
# Check if 'from pathlib import Path' and 'model_dir = Path' exist nearby
if 'SISATrainer(model_dir)' in content:
    # Find all occurrences and ensure model_dir is defined
    lines = content.split('\n')
    new_lines = []
    
    for i, line in enumerate(lines):
        if 'SISATrainer(model_dir)' in line:
            # Check if model_dir is defined in previous lines (within 10 lines)
            has_model_dir = False
            for j in range(max(0, i-10), i):
                if 'model_dir = Path(' in lines[j]:
                    has_model_dir = True
                    break
            
            if not has_model_dir:
                # Add model_dir definition before this line
                indent = len(line) - len(line.lstrip())
                new_lines.append(' ' * indent + 'from pathlib import Path')
                new_lines.append(' ' * indent + \"model_dir = Path('models/sisa')\")
                print(f'\\nAdded model_dir definition at line {i+1}')
                fixes_made += 1
        
        new_lines.append(line)
    
    content = '\\n'.join(new_lines)

# Write back
with open('app_enhanced.py', 'w', encoding='utf-8') as f:
    f.write(content)

print(f'\\n✅ Fixed {fixes_made} issues in app_enhanced.py')
print('All SISATrainer initializations now use: SISATrainer(model_dir)')
