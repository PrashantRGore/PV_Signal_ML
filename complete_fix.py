from pathlib import Path

dsm_file = Path('src/data/data_source_manager.py')
content = dsm_file.read_text(encoding='utf-8')

# Replace the ENTIRE load_faers_dataset method with working session state version
import re

# Find and replace the method
pattern = r'    def load_faers_dataset\(self\) -> Optional\[pd\.DataFrame\]:.*?(?=\n    def _compute_quarters)'

replacement = '''    def load_faers_dataset(self) -> Optional[pd.DataFrame]:
        """Load live FAERS data for user-selected date range"""
        st.sidebar.subheader("Live FAERS Quarterly Data")
        
        # Initialize session state
        if 'faers_configured' not in st.session_state:
            st.session_state.faers_configured = False
            st.session_state.faers_quarters = []
            st.session_state.faers_start = None
            st.session_state.faers_end = None
        
        # Get today's date
        today = datetime.now()
        
        # Date inputs
        col1, col2 = st.sidebar.columns(2)
        with col1:
            start_date = st.date_input(
                "Start Date",
                value=datetime(2024, 1, 1),
                min_value=datetime(2004, 1, 1),
                max_value=today,
                help="FAERS data available from 2004 onwards"
            )
        with col2:
            end_date = st.date_input(
                "End Date",
                value=today,
                min_value=datetime(2004, 1, 1),
                max_value=today,
                help="Select end date"
            )
        
        # Validate dates
        if start_date > end_date:
            st.sidebar.error("⚠️ End date must be >= Start date")
            return None
        
        # Configure button - SAVE TO SESSION STATE
        if st.sidebar.button("📅 Configure Date Range"):
            st.session_state.faers_configured = True
            st.session_state.faers_start = start_date
            st.session_state.faers_end = end_date
            st.session_state.faers_quarters = self._compute_quarters(start_date, end_date)
            st.rerun()  # Force rerun to show quarters
        
        # Check if configured
        if not st.session_state.faers_configured:
            st.sidebar.info("👆 Click 'Configure Date Range' to proceed")
            return None
        
        # Show quarters (from session state)
        quarters = st.session_state.faers_quarters
        st.sidebar.success(f"✅ Detected {len(quarters)} quarters: {', '.join([f'Q{q[1]}-{q[0]}' for q in quarters])}")
        
        # Quarter selection
        selected_quarters = st.sidebar.multiselect(
            "Select quarters to download:",
            options=quarters,
            default=quarters,
            format_func=lambda q: f"Q{q[1]}-{q[0]}"
        )
        
        # Download button - ACTUALLY DOWNLOAD WHEN CLICKED
        if st.sidebar.button("📡 Download & Process FAERS Data", type="primary"):
            if not selected_quarters:
                st.sidebar.error("⚠️ Select at least one quarter")
                return None
            
            st.sidebar.warning(f"🚀 DOWNLOADING {len(selected_quarters)} QUARTERS FROM FDA...")
            
            try:
                from src.data.faers_downloader import FAERSDownloader
                downloader = FAERSDownloader()
                
                # This will show progress bars in sidebar
                data = downloader.download_and_process(
                    quarters=selected_quarters,
                    start_date=st.session_state.faers_start,
                    end_date=st.session_state.faers_end
                )
                
                self.metadata = {
                    "source": "faers_live",
                    "quarters": selected_quarters,
                    "date_range": f"{st.session_state.faers_start} to {st.session_state.faers_end}",
                    "load_time": datetime.now().isoformat(),
                    "record_count": len(data)
                }
                
                st.sidebar.success(f"✅ Loaded {len(data):,} FAERS records")
                logger.info(f"FAERS loaded: {len(data)} records from {len(selected_quarters)} quarters")
                
                # Reset configuration for next download
                st.session_state.faers_configured = False
                
                return data
                
            except Exception as e:
                st.sidebar.error(f"❌ Download failed: {str(e)}")
                logger.error(f"FAERS download error: {e}", exc_info=True)
                return None
        
        return None

'''

content = re.sub(pattern, replacement, content, flags=re.DOTALL)
dsm_file.write_text(content, encoding='utf-8')
print('✅ COMPLETE FIX APPLIED - Session state will persist quarters!')
print('✅ Download button will now trigger actual download!')
