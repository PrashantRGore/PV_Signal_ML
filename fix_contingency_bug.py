from pathlib import Path

disp_file = Path('src/stats_engine/disproportionality.py')
content = disp_file.read_text(encoding='utf-8')

# Fix the contingency table calculation
old_calc = '''        drug_events = data.groupby(['drug_name', 'event_term']).size().reset_index(name='n11')

        total_cases = len(data['case_id'].unique())

        for _, row in drug_events.iterrows():
            drug = row['drug_name']
            event = row['event_term']
            n11 = row['n11']  # Drug + Event

            # Skip if below minimum threshold
            if n11 < self.min_case_count:
                continue

            # Compute contingency table
            n1_ = len(data[data['drug_name'] == drug]['case_id'].unique())  # Drug total
            n_1 = len(data[data['event_term'] == event]['case_id'].unique())  # Event total
            n10 = n1_ - n11  # Drug but not Event
            n01 = n_1 - n11  # Event but not Drug
            n00 = total_cases - n11 - n10 - n01  # Neither'''

new_calc = '''        # Count UNIQUE case_ids for each drug-event pair (not row count!)
        drug_events = data.groupby(['drug_name', 'event_term'])['case_id'].nunique().reset_index(name='n11')

        total_cases = data['case_id'].nunique()

        for _, row in drug_events.iterrows():
            drug = row['drug_name']
            event = row['event_term']
            n11 = row['n11']  # Drug + Event (unique cases)

            # Skip if below minimum threshold
            if n11 < self.min_case_count:
                continue

            # Compute contingency table (all using unique case counts)
            n1_ = data[data['drug_name'] == drug]['case_id'].nunique()  # Drug total
            n_1 = data[data['event_term'] == event]['case_id'].nunique()  # Event total
            n10 = n1_ - n11  # Drug but not Event
            n01 = n_1 - n11  # Event but not Drug
            n00 = total_cases - n11 - n10 - n01  # Neither
            
            # Safety check: ensure no negative values
            if n10 < 0 or n01 < 0 or n00 < 0:
                logger.warning(f"Negative contingency value for {drug}-{event}: n10={n10}, n01={n01}, n00={n00}")
                continue'''

content = content.replace(old_calc, new_calc)
disp_file.write_text(content, encoding='utf-8')

print('✅ FIXED: Changed groupby().size() to groupby().nunique()')
print('✅ Now counts UNIQUE case_ids, preventing negative values')
print('✅ Added safety check for negative values')
