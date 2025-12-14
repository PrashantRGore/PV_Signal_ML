# Search for the ML Validation section and replace SISA training code
import re

# Read current app
with open('app_enhanced.py', 'r', encoding='utf-8') as f:
    content = f.read()

# Find and replace SISA training section
# Look for the checkbox Enable SISA Training pattern
sisa_pattern = r"if st\.checkbox\('Enable SISA Training'.*?\n(?:.*?\n)*?.*?st\.error.*?Training failed.*?\)"

replacement = '''if st.checkbox('Enable SISA Training', key='sisa_training_enable'):
    st.info('**ML-based validation for signal prioritization**')
    
    n_shards = st.number_input('Number of Shards', min_value=5, max_value=20, value=10, 
                                help='More shards = better unlearning capability')
    
    if st.button('🧠 Train ML Model', key='train_sisa'):
        with st.spinner('Training SISA model...'):
            try:
                from src.ml.sisa_trainer import SISATrainer
                from pathlib import Path
                
                # Initialize trainer
                model_dir = Path('models/sisa')
                trainer = SISATrainer(model_dir)
                
                # Check if we have signals to train on
                if 'signals_df' in st.session_state and st.session_state.signals_df is not None:
                    signals_df = st.session_state.signals_df
                    
                    st.info(f'Training on {len(signals_df)} signals...')
                    
                    # Train model
                    results = trainer.train(signals_df, n_shards=n_shards)
                    
                    st.success('✅ SISA Model Trained Successfully!')
                    
                    # Display results
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        st.metric('AUC Score', f"{results['auc']:.4f}")
                    with col2:
                        st.metric('Training Samples', results['n_train'])
                    with col3:
                        st.metric('Features', results['n_features'])
                    
                    # Show classification report
                    with st.expander('📊 Detailed Classification Report'):
                        st.code(results['report'])
                    
                    # Predict on all signals
                    validated_signals = trainer.predict(signals_df)
                    st.session_state['validated_signals'] = validated_signals
                    
                    # Show validated signals
                    st.subheader('🎯 ML-Validated Signals')
                    high_confidence = validated_signals[validated_signals['ml_confidence'] >= 0.7]
                    st.write(f'High confidence signals: **{len(high_confidence)}**')
                    st.dataframe(high_confidence.head(20))
                    
                else:
                    st.error('⚠️ No signals found. Run signal detection first.')
                    
            except Exception as e:
                st.error(f'❌ Training failed: {str(e)}')
                import traceback
                st.code(traceback.format_exc())'''

# Try pattern matching replacement
if re.search(sisa_pattern, content, re.DOTALL):
    content = re.sub(sisa_pattern, replacement, content, flags=re.DOTALL)
    print('✅ Found and replaced SISA section using pattern matching')
else:
    print('⚠️ Pattern not found. Trying alternative approach...')
    
    # Alternative: Look for specific marker and insert after it
    if 'Enable SISA Training' in content:
        # Find the section
        lines = content.split('\n')
        new_lines = []
        skip_until = None
        replaced = False
        
        for i, line in enumerate(lines):
            if "st.checkbox('Enable SISA Training'" in line and not replaced:
                # Found start, now find the end of this block
                indent_level = len(line) - len(line.lstrip())
                
                # Add new code
                new_lines.append(replacement)
                
                # Skip old code until we find a line at same or lower indent
                j = i + 1
                while j < len(lines):
                    next_line = lines[j]
                    if next_line.strip() and not next_line.strip().startswith('#'):
                        next_indent = len(next_line) - len(next_line.lstrip())
                        if next_indent <= indent_level:
                            break
                    j += 1
                
                skip_until = j
                replaced = True
                continue
            
            if skip_until is not None and i < skip_until:
                continue
            
            new_lines.append(line)
        
        content = '\n'.join(new_lines)
        print('✅ Replaced SISA section using line-by-line approach')

# Write updated content
with open('app_enhanced.py', 'w', encoding='utf-8') as f:
    f.write(content)

print('\n✅ app_enhanced.py updated successfully!')
print('Backup saved as: app_enhanced.py.backup')
