from huggingface_hub import HfApi
from pathlib import Path

api = HfApi()
repo_id = 'PrashantRGore/pv-signal-ml-data'
upload_dir = Path(r'C:\Users\koreo\PV_Signal_ML\data\huggingface_upload')

print('📤 Adding standardized FAERS-format files to HuggingFace dataset...\n')

# Upload new files to root (not subdirectories)
files_to_upload = {
    'demo.parquet': 'Patient demographics (standardized)',
    'drug.parquet': 'Drug exposures (standardized)',
    'reac.parquet': 'Adverse events (standardized)',
    'outc.parquet': 'Case outcomes (standardized)',
    'README.md': 'Updated documentation'
}

for filename, description in files_to_upload.items():
    file_path = upload_dir / filename
    if file_path.exists():
        print(f'📤 Uploading {filename}...')
        print(f'   {description}')
        
        api.upload_file(
            path_or_fileobj=str(file_path),
            path_in_repo=filename,  # Upload to root, not subdirectory
            repo_id=repo_id,
            repo_type='dataset',
            commit_message=f'Add {filename} - standardized FAERS format for easier loading'
        )
        
        size_mb = file_path.stat().st_size / 1024 / 1024
        print(f'   ✅ Uploaded ({size_mb:.1f} MB)\n')

print('='*70)
print('🎉 SUCCESS! Dataset updated with standardized files')
print('='*70)
print(f'\nDataset URL: https://huggingface.co/datasets/{repo_id}')
print(f'\nNew files available at root level:')
print('  ✅ demo.parquet - Easy to load')
print('  ✅ drug.parquet - Standard naming')
print('  ✅ reac.parquet - FAERS compatible')
print('  ✅ outc.parquet - Clean structure')
print('\nOriginal files preserved in data/clean/ subdirectory')
