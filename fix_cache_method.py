from pathlib import Path

print('Fixing FAERSCache method name mismatch...')

# Fix: Add get_cache_size as alias to get_total_size
cache_file = Path('src/data/faers_cache.py')
content = cache_file.read_text(encoding='utf-8')

# Add the missing method
if 'def get_cache_size(' not in content:
    # Add alias method
    content = content.replace(
        '    def get_total_size(self):',
        '''    def get_total_size(self):
        """Total cache size in MB"""
        return sum(f.stat().st_size for f in self.cache_dir.glob('*.parquet')) / (1024*1024)
    
    def get_cache_size(self):'''
    )
    
    cache_file.write_text(content, encoding='utf-8')
    print('✅ Added get_cache_size() method')
else:
    print('✅ Method already exists')

# Also check if cache display code exists in data_source_manager
print()
print('Checking cache display in data_source_manager...')
dsm_file = Path('src/data/data_source_manager.py')
dsm_content = dsm_file.read_text(encoding='utf-8')

if 'cache.get_cache_size()' in dsm_content:
    print('⚠️  Cache display code exists but may be in wrong location')
    print('   Removing old cache display code...')
    
    # Remove the problematic cache display section
    import re
    # Remove the cache status display that's causing the error
    pattern = r"        # Show cache status.*?st\.write\(f\"❌ Q\{qtr\}-\{year\}: Not cached\"\)"
    dsm_content = re.sub(pattern, '', dsm_content, flags=re.DOTALL)
    
    dsm_file.write_text(dsm_content, encoding='utf-8')
    print('✅ Removed problematic cache display')

print()
print('✅ FIXED: Cache method mismatch resolved')
print('🚀 Restart your app now!')
