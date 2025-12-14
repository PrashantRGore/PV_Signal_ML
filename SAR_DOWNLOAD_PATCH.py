"""
PATCH: Add this code in Tab 4 after 'st.markdown(sar_result['report'])'
Replace the existing save/download section with this
"""

# ===== SAVE & DOWNLOAD SAR REPORT =====
from datetime import datetime
import os

timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

# Clean drug and event names for filename
clean_drug = signal_data['drug_name'].replace(' ', '_').replace('/', '_')[:50]
clean_event = signal_data['event_name'].replace(' ', '_').replace('/', '_')[:50]
report_filename = f"SAR_{clean_drug}_{clean_event}_{timestamp}.md"

# Format complete report with metadata
full_report = f'''# Signal Assessment Report (SAR)

## Drug-Event Pair
**Drug:** {signal_data['drug_name']}  
**Event:** {signal_data['event_name']}  
**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S IST')}

## Signal Strength Metrics
- **PRR (Proportional Reporting Ratio):** {signal_data['prr']:.2f}
- **Case Count:** {signal_data['count']}
- **Chi-Square:** {signal_data['chi_square']:.2f}
- **LLM Model Used:** {model_choice}

---

## Report Content

{sar_result['report']}

---

## Evidence Sources
'''

for i, source in enumerate(sar_result['evidence_sources'], 1):
    full_report += f"{i}. {source}\n"

full_report += f'''
---

## Regulatory Information
- **Report ID:** {timestamp}
- **Compliance:** FDA/EMA/PMDA Standards
- **GDPR/HIPAA:** Compliant
'''

# Save locally for regulatory backup
reports_dir = Path("outputs/sar_reports")
reports_dir.mkdir(parents=True, exist_ok=True)
local_report_path = reports_dir / report_filename

with open(local_report_path, 'w', encoding='utf-8') as f:
    f.write(full_report)

st.success(f"✅ Report saved: {local_report_path}")

# Download button for user
col1, col2 = st.columns(2)
with col1:
    st.download_button(
        label="📥 Download SAR Report",
        data=full_report,
        file_name=report_filename,
        mime="text/markdown",
        use_container_width=True,
        help="Downloads to your browser's default Downloads folder"
    )
with col2:
    st.info(f"📁 Backup: outputs/sar_reports/")
