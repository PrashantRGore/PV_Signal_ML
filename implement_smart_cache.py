from pathlib import Path

print()
print('IMPLEMENTING SMART DATE-AWARE CACHE...')

cache_code = '''"""
FAERS Cache with Date-Range Awareness
"""

import pandas as pd
from pathlib import Path
from datetime import datetime
import json
import hashlib

class FAERSCache:
    def __init__(self, cache_dir='data/cache/faers'):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.metadata_file = self.cache_dir / 'metadata.json'
        self.metadata = self._load_metadata()
    
    def _load_metadata(self):
        if self.metadata_file.exists():
            with open(self.metadata_file, 'r') as f:
                return json.load(f)
        return {}
    
    def _save_metadata(self):
        with open(self.metadata_file, 'w') as f:
            json.dump(self.metadata, f, indent=2)
    
    def get_cache_key(self, year, quarter):
        """Simple key for full quarter data"""
        return f'{year}Q{quarter}'
    
    def get_cache_path(self, year, quarter):
        """Path for cached quarter file"""
        return self.cache_dir / f'{year}Q{quarter}.parquet'
    
    def is_cached(self, year, quarter):
        """Check if full quarter is cached"""
        return self.get_cache_path(year, quarter).exists()
    
    def save(self, year, quarter, data, start_date=None, end_date=None):
        """Save quarter data with optional date range metadata"""
        path = self.get_cache_path(year, quarter)
        data.to_parquet(path, index=False, compression='snappy')
        
        key = self.get_cache_key(year, quarter)
        self.metadata[key] = {
            'cached_at': datetime.now().isoformat(),
            'records': len(data),
            'size_mb': path.stat().st_size / (1024*1024),
            'start_date': str(start_date) if start_date else None,
            'end_date': str(end_date) if end_date else None,
            'full_quarter': True  # Full quarter data, can be filtered
        }
        self._save_metadata()
        return path
    
    def load(self, year, quarter):
        """Load full quarter data (can be filtered by date after loading)"""
        path = self.get_cache_path(year, quarter)
        if path.exists():
            return pd.read_parquet(path)
        return None
    
    def get_info(self, year, quarter):
        """Get cache metadata"""
        key = self.get_cache_key(year, quarter)
        return self.metadata.get(key, {})
    
    def list_cached(self):
        """List all cached quarters"""
        return [(int(k.split('Q')[0]), int(k.split('Q')[1])) 
                for k in self.metadata.keys() if 'Q' in k]
    
    def get_total_size(self):
        """Total cache size in MB"""
        return sum(f.stat().st_size for f in self.cache_dir.glob('*.parquet')) / (1024*1024)
    
    def clear_quarter(self, year, quarter):
        """Remove quarter from cache"""
        path = self.get_cache_path(year, quarter)
        if path.exists():
            path.unlink()
        
        key = self.get_cache_key(year, quarter)
        if key in self.metadata:
            del self.metadata[key]
            self._save_metadata()
'''

cache_file = Path('src/data/faers_cache.py')
cache_file.write_text(cache_code, encoding='utf-8')
print('✅ Updated: Smart cache (caches full quarters, filters on load)')

print()
print('📋 Cache Strategy:')
print('  • Downloads full quarter from FDA (e.g., all of Q1-2025)')
print('  • Caches to: data/cache/faers/2025Q1.parquet')
print('  • NOT in Windows Downloads folder')
print('  • Date filtering applied AFTER loading from cache')
print('  • Same quarter = instant load, no duplicate download')
