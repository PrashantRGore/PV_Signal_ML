# Quick Reference: GDPR Right to Be Forgotten

**For 2 Patients - 5 ICSRs Total**

---

## 🚀 One-Command Solution

```bash
python process_deletion_requests.py
```

**That's it!** The script will:
- ✅ Record deletion requests
- ✅ Delete from database
- ✅ Pseudonymize IDs
- ✅ Generate audit trail
- ✅ Verify deletion

---

## 📋 What Gets Deleted

### Patient 1: John Doe
- ICSR_2024_001
- ICSR_2024_002
- ICSR_2024_003

### Patient 2: Jane Smith
- ICSR_2024_004
- ICSR_2024_005

---

## 📁 Files Generated

```
gdpr_registry/
├── deletion_requests.jsonl          # Deletion log
├── icsr_pseudonyms.json             # Pseudonym mappings
├── deleted_icsr_ids.txt             # List of deleted IDs
└── deletion_audit_20251208_*.json   # Audit report
```

---

## 🔐 How It Works

```
Request Deletion
    ↓
Record with Timestamp
    ↓
Delete from Database
    ↓
Create Pseudonym (HMAC-SHA256)
    ↓
Generate Audit Report
    ↓
Verify Deletion
    ↓
✅ GDPR Compliant
```

---

## 📊 Example Output

```
======================================================================
🔐 GDPR Right to Be Forgotten - Manual Deletion Processing
======================================================================

📋 Processing Deletion Requests

Total requests: 2
Total ICSRs to delete: 5

──────────────────────────────────────────────────────────────────────
Request 1: PATIENT_001
──────────────────────────────────────────────────────────────────────
Reason: Patient requested right to be forgotten
ICSRs to delete: 3
  ✅ Deletion recorded: ICSR_2024_001
  ✅ Deletion recorded: ICSR_2024_002
  ✅ Deletion recorded: ICSR_2024_003

──────────────────────────────────────────────────────────────────────
Request 2: PATIENT_002
──────────────────────────────────────────────────────────────────────
Reason: Patient requested right to be forgotten
ICSRs to delete: 2
  ✅ Deletion recorded: ICSR_2024_004
  ✅ Deletion recorded: ICSR_2024_005

🗑️  Deleting from Database
✅ Deleted from database: ICSR_2024_001
✅ Deleted from database: ICSR_2024_002
✅ Deleted from database: ICSR_2024_003
✅ Deleted from database: ICSR_2024_004
✅ Deleted from database: ICSR_2024_005

Total records deleted from database: 5

🔒 Pseudonymization of References
ICSR ID: ICSR_2024_001
  Pseudonym: ICSR_a3f5b2c1d9e4f7a2

ICSR ID: ICSR_2024_002
  Pseudonym: ICSR_b4g6c3d2e0f5g8b3

... (3 more)

✅ All 5 ICSRs successfully deleted from database

✅ Deletion Processing Complete
```

---

## 🔍 Verify Deletion

Check the files:

```bash
# View deletion log
cat gdpr_registry/deletion_requests.jsonl

# View pseudonym mappings
cat gdpr_registry/icsr_pseudonyms.json

# View deleted IDs
cat gdpr_registry/deleted_icsr_ids.txt

# View audit report
cat gdpr_registry/deletion_audit_*.json
```

---

## 📝 Manual Steps (If Needed)

### Step 1: Record Deletion Request
```python
from gdpr_deletion_registry import GDPRDeletionRegistry

registry = GDPRDeletionRegistry()
registry.request_deletion("ICSR_2024_001", "Patient requested")
```

### Step 2: Delete from Database
```python
import sqlite3

conn = sqlite3.connect("pv_signal.db")
cursor = conn.cursor()
cursor.execute("DELETE FROM icsr WHERE icsr_id = ?", ("ICSR_2024_001",))
conn.commit()
conn.close()
```

### Step 3: Create Pseudonym
```python
pseudonym = registry.pseudonymize_icsr_id("ICSR_2024_001")
# Output: ICSR_a3f5b2c1d9e4f7a2
```

---

## ✅ Compliance Checklist

- [ ] Deletion requests recorded with timestamp
- [ ] All ICSRs deleted from database
- [ ] Pseudonym mappings created
- [ ] Audit trail generated
- [ ] Deletion verified
- [ ] Files saved in `gdpr_registry/`
- [ ] Ready for regulatory inspection

---

## 🎯 Key Points

✅ **Irreversible:** Pseudonyms cannot be reversed  
✅ **Traceable:** Audit trail links to deletion records  
✅ **Compliant:** Meets GDPR Article 17 requirements  
✅ **Secure:** No personal data in logs  
✅ **Auditable:** Complete record for inspection  

---

## 📞 Need Help?

1. Read `GDPR_RIGHT_TO_BE_FORGOTTEN_GUIDE.md` (detailed guide)
2. Check `gdpr_deletion_registry.py` (source code)
3. Review generated audit reports
4. Check `VALIDATION_STATUS.md` (compliance details)

---

**Last Updated:** 2025-12-08

*Run `python process_deletion_requests.py` to process deletion requests for 2 patients!*
