from pathlib import Path
from datetime import datetime

file = Path('src/data/data_source_manager.py')
content = file.read_text(encoding='utf-8')

# Find and replace the method
old_method = '''    def _compute_quarters(self, start_date: datetime, end_date: datetime) -> List[Tuple[int, int]]:
        """
        Compute all quarters overlapping the date range
        Returns: List of (year, quarter) tuples
        """
        quarters = []
        current = start_date
        
        while current <= end_date:
            year = current.year
            month = current.month
            quarter = (month - 1) // 3 + 1
            
            if (year, quarter) not in quarters:
                quarters.append((year, quarter))
            
            # Move to next quarter
            if quarter == 4:
                current = datetime(year + 1, 1, 1)
            else:
                current = datetime(year, (quarter * 3) + 1, 1)
        
        return quarters'''

new_method = '''    def _compute_quarters(self, start_date, end_date) -> List[Tuple[int, int]]:
        """
        Compute all quarters overlapping the date range
        Returns: List of (year, quarter) tuples
        """
        # Convert date objects to datetime for consistent comparison
        if not isinstance(start_date, datetime):
            start_date = datetime.combine(start_date, datetime.min.time())
        if not isinstance(end_date, datetime):
            end_date = datetime.combine(end_date, datetime.min.time())
        
        quarters = []
        current = start_date
        
        while current <= end_date:
            year = current.year
            month = current.month
            quarter = (month - 1) // 3 + 1
            
            if (year, quarter) not in quarters:
                quarters.append((year, quarter))
            
            # Move to next quarter
            if quarter == 4:
                current = datetime(year + 1, 1, 1)
            else:
                current = datetime(year, (quarter * 3) + 1, 1)
        
        return quarters'''

content = content.replace(old_method, new_method)
file.write_text(content, encoding='utf-8')
print('✅ File updated successfully!')
