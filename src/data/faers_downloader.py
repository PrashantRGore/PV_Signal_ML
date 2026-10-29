"""
FAERS Downloader - Download quarterly FAERS data from FDA
"""

import requests
from pathlib import Path
from typing import Optional
import zipfile
from src.utils.logger import setup_logger

logger = setup_logger(__name__)


class FAERSDownloader:
    """Download FAERS quarterly data from FDA."""
    
    BASE_URL = "https://fis.fda.gov/content/Exports/"
    
    def __init__(self, cache_dir: Path = Path("data/faers_cache")):
        self.cache_dir = cache_dir
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        
    def download_quarter(self, quarter: str) -> Optional[Path]:
        """
        Download FAERS quarterly data (wrapper for compatibility).
        
        Args:
            quarter: Quarter string like "24Q3" or "2024Q3"
            
        Returns:
            Path to downloaded zip file
        """
        # Parse quarter string
        if len(quarter) == 4:  # Format: "24Q3"
            year = f"20{quarter[:2]}"
            q = quarter[2:]
        elif len(quarter) == 6:  # Format: "2024Q3"
            year = quarter[:4]
            q = quarter[4:]
        else:
            logger.error(f"Invalid quarter format: {quarter}. Use '24Q3' or '2024Q3'")
            return None
        
        return self.download(year, q)
    
    def download(self, year: str, quarter: str) -> Optional[Path]:
        """
        Download FAERS data for specific year and quarter.
        
        Args:
            year: Year (e.g., "2024")
            quarter: Quarter (e.g., "Q3")
            
        Returns:
            Path to downloaded zip file
        """
        # Normalize quarter format
        if not quarter.startswith('Q'):
            quarter = f"Q{quarter}"
        
        # Check cache first and verify zip integrity
        cache_file = self.cache_dir / f"faers_{year}{quarter}.zip"
        if cache_file.exists():
            if zipfile.is_zipfile(cache_file):
                logger.info(f"Using valid cached FAERS data: {cache_file}")
                return cache_file
            else:
                logger.warning(f"Cached FAERS file {cache_file} is corrupted. Removing and re-downloading...")
                try:
                    cache_file.unlink()
                except Exception as e:
                    logger.error(f"Failed to remove corrupt cache file: {e}")
        
        # Construct download URL
        # Example: https://fis.fda.gov/content/Exports/faers_ascii_2024q3.zip
        filename = f"faers_ascii_{year.lower()}{quarter.lower()}.zip"
        url = f"{self.BASE_URL}{filename}"
        
        logger.info(f"Downloading FAERS {year}{quarter} from {url}")
        
        try:
            response = requests.get(url, stream=True, timeout=300)
            
            if response.status_code != 200:
                logger.error(f"Download failed: HTTP {response.status_code}")
                return None
            
            # Download with progress
            total_size = int(response.headers.get('content-length', 0))
            downloaded = 0
            
            tmp_file = cache_file.with_suffix('.tmp')
            with open(tmp_file, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
                        downloaded += len(chunk)
                        
                        if total_size > 0:
                            progress = (downloaded / total_size) * 100
                            if downloaded % (1024 * 1024 * 10) == 0:  # Log every 10MB
                                logger.info(f"Progress: {progress:.1f}% ({downloaded / (1024*1024):.1f} MB)")
            
            if zipfile.is_zipfile(tmp_file):
                tmp_file.replace(cache_file)
                logger.info(f"✅ Downloaded and verified FAERS {year}{quarter} to {cache_file}")
                return cache_file
            else:
                logger.error(f"Downloaded file for {year}{quarter} failed zip validation.")
                tmp_file.unlink(missing_ok=True)
                return None
            
        except requests.exceptions.RequestException as e:
            logger.error(f"Download failed: {e}")
            return None
        except Exception as e:
            logger.error(f"Unexpected error: {e}")
            return None
    
    def list_available_quarters(self) -> list:
        """List quarters available in cache."""
        cached = list(self.cache_dir.glob("faers_*.zip"))
        return [f.stem.replace("faers_", "") for f in cached]
