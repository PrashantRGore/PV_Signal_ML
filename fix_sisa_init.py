import re

# Read the app
with open('app_enhanced.py', 'r', encoding='utf-8') as f:
    content = f.read()

# Fix 1: Replace incorrect SISATrainer initialization
content = re.sub(
    r'trainer = SISATrainer\(num_shards=n_shards\)',
    'trainer = SISATrainer(model_dir)',
    content
)

# Fix 2: Also look for any variant
content = re.sub(
    r'SISATrainer\(\s*num_shards\s*=\s*[^)]+\)',
    'SISATrainer(model_dir)',
    content
)

# Write back
with open('app_enhanced.py', 'w', encoding='utf-8') as f:
    f.write(content)

print('✅ Fixed SISATrainer initialization in app_enhanced.py')
