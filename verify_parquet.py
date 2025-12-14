import pyarrow.parquet as pq

files = ['data/huggingface_upload/demo.parquet', 
         'data/huggingface_upload/drug.parquet',
         'data/huggingface_upload/outc.parquet', 
         'data/huggingface_upload/reac.parquet']

for f in files:
    try:
        table = pq.read_table(f)
        print(f'✓ {f}: {len(table)} rows - VALID')
    except Exception as e:
        print(f'✗ {f}: CORRUPTED - {e}')
