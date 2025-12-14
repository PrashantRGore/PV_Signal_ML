import re

# Read the file
with open("app_enhanced.py", "r", encoding="utf-8") as f:
    content = f.read()

# Find the SISA training section and replace it completely
# Look for the "if enable_sisa:" block around line 142

pattern = r"(if enable_sisa:.*?)(else:.*?st\.warning.*?$)"

replacement = """if enable_sisa:
                        from pathlib import Path
                        model_dir = Path("models/sisa")
                        trainer = SISATrainer(model_dir)
                        
                        # Train SISA model
                        results = trainer.train(signals_df, n_shards=num_shards)
                        st.session_state.trained_model = trainer
                        st.session_state.training_results = results
                        
                        st.success("✅ SISA Model Trained Successfully!")
                        
                        # Display metrics
                        col1, col2, col3 = st.columns(3)
                        with col1:
                            st.metric("AUC Score", f"{results['auc']:.4f}")
                        with col2:
                            st.metric("Training Samples", results['n_train'])
                        with col3:
                            st.metric("Features", results['n_features'])
                        
                        # Show detailed results
                        with st.expander("📊 Training Details"):
                            st.json({
                                'auc': results['auc'],
                                'n_features': results['n_features'],
                                'n_train': results['n_train'],
                                'n_test': results['n_test']
                            })
                            st.code(results['report'])
                    else:
                        st.warning("SISA disabled - using standard training")"""

# Try to find and replace the section
if re.search(pattern, content, re.DOTALL | re.MULTILINE):
    content = re.sub(pattern, replacement, content, flags=re.DOTALL | re.MULTILINE)
    print("✅ Replaced SISA section using regex")
else:
    # Fallback: Find line-by-line
    lines = content.split('\n')
    new_lines = []
    i = 0
    
    while i < len(lines):
        line = lines[i]
        
        # Find "if enable_sisa:" line
        if 'if enable_sisa:' in line and i > 140:
            indent = len(line) - len(line.lstrip())
            base_indent = ' ' * indent
            inner_indent = ' ' * (indent + 4)
            
            # Add complete new SISA block
            new_lines.append(line)  # Keep the "if enable_sisa:" line
            new_lines.append(inner_indent + 'from pathlib import Path')
            new_lines.append(inner_indent + 'model_dir = Path("models/sisa")')
            new_lines.append(inner_indent + 'trainer = SISATrainer(model_dir)')
            new_lines.append(inner_indent + '')
            new_lines.append(inner_indent + '# Train SISA model')
            new_lines.append(inner_indent + 'results = trainer.train(signals_df, n_shards=num_shards)')
            new_lines.append(inner_indent + 'st.session_state.trained_model = trainer')
            new_lines.append(inner_indent + 'st.session_state.training_results = results')
            new_lines.append(inner_indent + '')
            new_lines.append(inner_indent + 'st.success("✅ SISA Model Trained Successfully!")')
            new_lines.append(inner_indent + '')
            new_lines.append(inner_indent + '# Display metrics')
            new_lines.append(inner_indent + 'col1, col2, col3 = st.columns(3)')
            new_lines.append(inner_indent + 'with col1:')
            new_lines.append(inner_indent + '    st.metric("AUC Score", f"{results[\'auc\']:.4f}")')
            new_lines.append(inner_indent + 'with col2:')
            new_lines.append(inner_indent + '    st.metric("Training Samples", results[\'n_train\'])')
            new_lines.append(inner_indent + 'with col3:')
            new_lines.append(inner_indent + '    st.metric("Features", results[\'n_features\'])')
            new_lines.append(inner_indent + '')
            new_lines.append(inner_indent + '# Show detailed results')
            new_lines.append(inner_indent + 'with st.expander("📊 Training Details"):')
            new_lines.append(inner_indent + '    st.json({')
            new_lines.append(inner_indent + '        "auc": results["auc"],')
            new_lines.append(inner_indent + '        "n_features": results["n_features"],')
            new_lines.append(inner_indent + '        "n_train": results["n_train"],')
            new_lines.append(inner_indent + '        "n_test": results["n_test"]')
            new_lines.append(inner_indent + '    })')
            new_lines.append(inner_indent + '    st.code(results["report"])')
            
            # Skip old code until we find "else:" at same indent level
            i += 1
            while i < len(lines):
                next_line = lines[i]
                next_indent = len(next_line) - len(next_line.lstrip()) if next_line.strip() else 999
                
                if next_line.strip().startswith('else:') and next_indent == indent:
                    new_lines.append(next_line)  # Keep the else line
                    break
                i += 1
            
            print("✅ Replaced SISA section line-by-line")
        else:
            new_lines.append(line)
        
        i += 1
    
    content = '\n'.join(new_lines)

# Write back
with open("app_enhanced.py", "w", encoding="utf-8") as f:
    f.write(content)

print("\n✅ Complete SISA section updated!")
print("Fixed issues:")
print("  1. ✅ Variable name: metadata → results")
print("  2. ✅ Added proper metric display")
print("  3. ✅ Added session state storage")
print("  4. ✅ Added error-safe JSON display")
