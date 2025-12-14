from pathlib import Path

dsm_file = Path('src/data/data_source_manager.py')
content = dsm_file.read_text(encoding='utf-8')

# Add cache display after Configure button
old_section = '''        # Show quarters (from session state)
        quarters = st.session_state.faers_quarters
        st.sidebar.success(f"✅ Detected {len(quarters)} quarters: {', '.join([f'Q{q[1]}-{q[0]}' for q in quarters])}")'''

new_section = '''        # Show quarters (from session state)
        quarters = st.session_state.faers_quarters
        st.sidebar.success(f"✅ Detected {len(quarters)} quarters: {', '.join([f'Q{q[1]}-{q[0]}' for q in quarters])}")
        
        # Show cache status
        from src.data.faers_cache import FAERSCache
        cache = FAERSCache()
        cached_quarters = [q for q in quarters if cache.is_cached(q[0], q[1])]
        
        if cached_quarters:
            st.sidebar.info(f"💾 {len(cached_quarters)}/{len(quarters)} quarters available in cache")
            with st.sidebar.expander(f"📦 View Cache ({cache.get_cache_size():.1f} MB total)"):
                for year, qtr in quarters:
                    if cache.is_cached(year, qtr):
                        info = cache.get_cache_info(year, qtr)
                        st.write(f"✅ Q{qtr}-{year}: {info.get('record_count', 0):,} records ({info.get('file_size_mb', 0):.1f} MB)")
                    else:
                        st.write(f"❌ Q{qtr}-{year}: Not cached")'''

content = content.replace(old_section, new_section)
dsm_file.write_text(content, encoding='utf-8')

print('✅ Added cache status display to UI')
print('✅ Users will see which quarters are cached')
