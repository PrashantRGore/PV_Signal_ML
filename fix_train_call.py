# Read the file
with open("app_enhanced.py", "r", encoding="utf-8") as f:
    lines = f.readlines()

# Find and fix line 146 (index 145)
new_lines = []
for i, line in enumerate(lines):
    if i == 145 and "metadata = trainer.train(X, y, case_ids)" in line:
        # Replace with correct call using signals_df
        indent = " " * 24  # Match indentation
        new_lines.append(indent + "# Get n_shards from the UI input\n")
        new_lines.append(indent + "results = trainer.train(signals_df, n_shards=num_shards)\n")
        print(f"✅ Fixed line {i+1} - Updated train() call")
    else:
        new_lines.append(line)

# Write back
with open("app_enhanced.py", "w", encoding="utf-8") as f:
    f.writelines(new_lines)

print("✅ Fixed train() method call")
