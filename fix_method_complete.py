from pathlib import Path
import re

file = Path('src/data/data_source_manager.py')
content = file.read_text(encoding='utf-8')

# Use regex to find and replace the entire load_faers_dataset method
pattern = r'(    def load_faers_dataset\(self\) -> Optional\[pd\.DataFrame\]:.*?)(    def \w+|class \w+|\Z)'

new_method = '''    def load_faers_dataset(self) -> Optional[pd.DataFrame]:
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
        quarters = self._compute_quarters(start_date, end_date)
        st.sidebar.write(f"**Detected Quarters:** {len(quarters)}")
        
        selected_quarters = st.sidebar.multiselect(
            "Confirm quarters to download:",
            options=quarters,
            default=quarters,
            format_func=lambda q: f"Q{q[1]}-{q[0]}"
        )
        
        load_faers = st.sidebar.button("📡 Download & Process FAERS Data")
        
        if load_faers and selected_quarters:
            try:
                downloader = FAERSDownloader()
                with st.spinner(f"Downloading {len(selected_quarters)} quarters..."):
                    data = downloader.download_and_process(
                        quarters=selected_quarters,
                        start_date=start_date,
                        end_date=end_date
                    )
                
                self.metadata = {
                    "source": "faers_live",
                    "quarters": selected_quarters,
                    "date_range": f"{start_date} to {end_date}",
                    "load_time": datetime.now().isoformat(),
                    "record_count": len(data)
                }
                
                st.sidebar.success(f"✅ Loaded {len(data):,} FAERS records")
                logger.info(f"FAERS data loaded: {len(data)} records from {len(selected_quarters)} quarters")
                return data
                
            except Exception as e:
                st.sidebar.error(f"❌ FAERS download failed: {str(e)}")
                logger.error(f"FAERS load failed: {e}")
                return None

'''

# Replace using regex with DOTALL flag
content = re.sub(pattern, new_method + r'\2', content, flags=re.DOTALL)

file.write_text(content, encoding='utf-8')
print('✅ Method completely replaced with Configure Date button!')
