from pathlib import Path

file = Path('src/data/data_source_manager.py')
content = file.read_text(encoding='utf-8')

# Find and replace the entire load_faers_dataset method
old_section = '''    def load_faers_dataset(self) -> Optional[pd.DataFrame]:
        """Load live FAERS data for user-selected date range"""
        st.sidebar.subheader("Live FAERS Quarterly Data")
        
        # Get today's date - works for any year (2025, 2026, 2027...)
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
            )
        
        # Compute quarters
        quarters = self._compute_quarters(start_date, end_date)'''

new_section = '''    def load_faers_dataset(self) -> Optional[pd.DataFrame]:
        """Load live FAERS data for user-selected date range"""
        st.sidebar.subheader("Live FAERS Quarterly Data")
        
        # Get today's date - works for any year (2025, 2026, 2027...)
        today = datetime.now()
        
        # Date inputs - independent selection
        col1, col2 = st.sidebar.columns(2)
        with col1:
            start_date = st.date_input(
                "Start Date",
                value=datetime(2024, 1, 1),
                min_value=datetime(2004, 1, 1),  # FAERS started ~2004
                max_value=today,
                help="FAERS data available from 2004 onwards"
            )
        with col2:
            end_date = st.date_input(
                "End Date",
                value=today,  # Default to today
                min_value=datetime(2004, 1, 1),
                max_value=today,
                help="Select end date"
            )
        
        # Validate dates before proceeding
        if start_date > end_date:
            st.sidebar.error("⚠️ End date must be >= Start date")
            return None
        
        # Configure button
        configure_dates = st.sidebar.button("📅 Configure Date Range", help="Click after selecting dates")
        
        if not configure_dates:
            st.sidebar.info("👆 Click 'Configure Date Range' to proceed")
            return None
        
        # Compute quarters after validation
        quarters = self._compute_quarters(start_date, end_date)'''

content = content.replace(old_section, new_section)
file.write_text(content, encoding='utf-8')
print('✅ Configure Date button added successfully!')
