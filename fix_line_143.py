# Read the file
with open('app_enhanced.py', 'r', encoding='utf-8') as f:
    lines = f.readlines()

# Find line 143 and add model_dir definition before it
new_lines = []
for i, line in enumerate(lines):
    # When we hit line 142 (if enable_sisa:), add model_dir on next line
    if i == 142 and 'trainer = SISATrainer(model_dir)' in line:
        # Add model_dir definition with proper indentation
        indent = ' ' * 24  # Match the indentation of line 143
        new_lines.append(indent + 'from pathlib import Path\n')
        new_lines.append(indent + \"model_dir = Path('models/sisa')\n\")
        print(f'✅ Added model_dir definition before line {i+1}')
    
    new_lines.append(line)

# Write back
with open('app_enhanced.py', 'w', encoding='utf-8') as f:
    f.writelines(new_lines)

print('✅ Fixed line 143 - model_dir now defined')
