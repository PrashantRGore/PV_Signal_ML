from pathlib import Path

app_file = Path('app_enhanced.py')
content = app_file.read_text(encoding='utf-8')

# Add data validation before compute_signals
old_detection = '''        st.header("Statistical Signal Detection")

        if st.button("🚀 Run Signal Detection", type="primary"):
            with st.spinner("Computing disproportionality signals..."):
                stats_engine = DisproportionalityAnalysis()
                signals = stats_engine.compute_signals(st.session_state.raw_data)
                st.session_state.signals = signals'''

new_detection = '''        st.header("Statistical Signal Detection")

        if st.button("🚀 Run Signal Detection", type="primary"):
            # Validate data structure
            data = st.session_state.raw_data
            
            st.info(f"📊 Data shape: {data.shape[0]:,} rows × {data.shape[1]} columns")
            st.write("**Columns:**", list(data.columns))
            
            # Check required columns
            required_cols = ['drug_name', 'event_term']
            missing_cols = [col for col in required_cols if col not in data.columns]
            
            if missing_cols:
                st.error(f"❌ Missing required columns: {missing_cols}")
                st.stop()
            
            # Check for empty values
            empty_drugs = data['drug_name'].isna().sum()
            empty_events = data['event_term'].isna().sum()
            
            if empty_drugs > 0 or empty_events > 0:
                st.warning(f"⚠️ Found {empty_drugs} empty drugs, {empty_events} empty events. Cleaning...")
                data = data.dropna(subset=['drug_name', 'event_term'])
                st.session_state.raw_data = data
            
            # Show sample
            with st.expander("🔍 View Data Sample"):
                st.dataframe(data.head(100))
            
            with st.spinner("Computing disproportionality signals..."):
                stats_engine = DisproportionalityAnalysis()
                signals = stats_engine.compute_signals(data)
                st.session_state.signals = signals'''

content = content.replace(old_detection, new_detection)
app_file.write_text(content, encoding='utf-8')

print('✅ Added data validation before signal detection')
print('✅ Will show data structure and clean before processing')
