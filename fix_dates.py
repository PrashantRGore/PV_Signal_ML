from pathlib import Path
from datetime import datetime

file = Path('src/data/data_source_manager.py')
content = file.read_text(encoding='utf-8')

# Fix the date defaults - make it flexible and future-proof
old_dates = '''        with col1:
            start_date = st.date_input(
                "Start Date",
                value=datetime(2024, 1, 1),
                min_value=datetime(2012, 1, 1),
                max_value=datetime.now()
            )
        with col2:
            end_date = st.date_input(
                "End Date",
                value=datetime(2024, 12, 31),
                min_value=start_date,
                max_value=datetime.now()
            )'''

new_dates = '''        # Get today's date - works for any year (2025, 2026, 2027...)
        today = datetime.now()
        
        with col1:
            start_date = st.date_input(
                "Start Date",
                value=datetime(2024, 1, 1),
                min_value=datetime(2004, 1, 1),  # FAERS started ~2004
                max_value=today,  # Cannot select future dates
                help="FAERS data available from 2004 onwards"
            )
        with col2:
            # End date: must be >= start_date AND <= today (no future dates allowed)
            default_end = datetime(2024, 12, 31) if datetime(2024, 12, 31) <= today else today
            end_date = st.date_input(
                "End Date",
                value=default_end,
                min_value=start_date,
                max_value=today,  # Always today's date, no future
                help="Cannot be before start date or after today"
            )'''

content = content.replace(old_dates, new_dates)
file.write_text(content, encoding='utf-8')
print('✅ Date inputs updated - future-proof and today-limited!')
