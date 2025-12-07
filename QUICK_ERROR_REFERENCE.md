# Quick Error Reference Guide

**Quick lookup for common errors and solutions**

---

## 🔴 Critical Errors

### Error: `404 Not Found` on `/signal-report`
**Cause:** API endpoint doesn't exist  
**Fix:** Restart API with latest code
```bash
python -m uvicorn api:app --host 127.0.0.1 --port 8000 --reload
```

### Error: `AttributeError: 'PVSignalRAGLangChain' object has no attribute 'generate_sar'`
**Cause:** Wrong method name  
**Fix:** Use `generate_sar_report()` instead of `generate_sar()`  
**File:** `api.py` line 83

### Error: `FileNotFoundError: [Errno 2] No such file or directory: 'sar_reports/full_signals_1M.csv'`
**Cause:** Required CSV file missing  
**Fix:** Create dummy data or download FAERS data
```bash
python faers_build_signals.py 2024-01-01 2024-03-31
```

### Error: `ConnectionRefusedError: [Errno 111] Connection refused` (Ollama)
**Cause:** Ollama service not running  
**Fix:** Start Ollama in separate terminal
```bash
ollama serve
# In another terminal:
ollama pull llama3.2:3b
```

---

## 🟡 Common Warnings

### Warning: `DeprecationWarning: datetime.utcnow() is deprecated`
**Cause:** Python 3.13 deprecation  
**Impact:** None - code still works  
**Fix:** Planned for Phase 2 (not critical)

### Warning: `Markdown file not found on disk`
**Cause:** Report generation didn't create markdown file  
**Fix:** Check if JSON was created, verify template exists
```bash
ls -la sar_reports/reports/
ls -la templates/signal_report_template.md
```

---

## 🟢 Quick Fixes

### Fix 1: Restart Everything
```bash
# Terminal 1: Start API
python -m uvicorn api:app --host 127.0.0.1 --port 8000 --reload

# Terminal 2: Start Streamlit
streamlit run pv_fullstack.py
```

### Fix 2: Clear Caches
```bash
# Clear Streamlit cache
rm -rf ~/.streamlit/cache

# Clear Python cache
find . -type d -name __pycache__ -exec rm -rf {} +
find . -type f -name "*.pyc" -delete
```

### Fix 3: Verify Dependencies
```bash
# Check all packages installed
pip list | grep -E "pandas|numpy|streamlit|fastapi|langchain|ollama"

# Reinstall if needed
pip install -r requirements.txt --upgrade
```

### Fix 4: Check Ports
```bash
# Check if ports are in use
netstat -ano | findstr :8000
netstat -ano | findstr :8501

# Kill process if needed
taskkill /PID <PID> /F
```

---

## 📊 Error Diagnosis Flowchart

```
Error occurs
    ↓
Is it a 404 error?
    ├─ YES → Restart API with --reload
    └─ NO → Continue
    ↓
Is it a FileNotFoundError?
    ├─ YES → Check file exists, create if needed
    └─ NO → Continue
    ↓
Is it a ConnectionError?
    ├─ YES → Check Ollama/API/Database running
    └─ NO → Continue
    ↓
Is it an AttributeError?
    ├─ YES → Check method names and imports
    └─ NO → Continue
    ↓
Check logs and error message
    ↓
Review CODE_REVIEW_AND_FIXES.md for similar issues
    ↓
If not found, check:
    - File paths are correct
    - CSV columns match expected names
    - Data types are correct
    - Required directories exist
```

---

## 🛠️ Maintenance Commands

### Daily Checks
```bash
# Check API is running
curl http://127.0.0.1:8000/

# Check Streamlit is running
curl http://127.0.0.1:8501/

# Check audit logs
tail -f audit_logs/audit.log

# Check GDPR registry
ls -la gdpr_registry/
```

### Weekly Checks
```bash
# Verify all CSV files exist
ls -la sar_reports/*.csv

# Check disk space
df -h

# Verify backups
ls -la backups/

# Check for errors in logs
grep -i error audit_logs/*.log
```

### Monthly Checks
```bash
# Update dependencies
pip install -r requirements.txt --upgrade

# Run full test suite
pytest tests/

# Generate compliance report
python -c "from audit_logging import AuditLogger; logger = AuditLogger(); print(logger.get_audit_report())"

# Backup data
tar -czf backups/backup_$(date +%Y%m%d).tar.gz sar_reports/ gdpr_registry/ audit_logs/
```

---

## 📞 Support Resources

### Documentation Files
- `README.md` – Project overview
- `ANALYSIS_AND_COMPLIANCE_REPORT.md` – Compliance details
- `CODE_REVIEW_AND_FIXES.md` – Detailed error analysis
- `IMPLEMENTATION_SUMMARY.md` – Implementation details
- `DEPLOYMENT_CHECKLIST.md` – Deployment verification

### Key Files to Check
- `api.py` – API endpoints
- `rag_langchain.py` – RAG pipeline
- `pv_fullstack.py` – Streamlit UI
- `audit_logging.py` – Audit trail
- `gdpr_deletion_registry.py` – GDPR compliance

### External Resources
- **Streamlit Docs:** https://docs.streamlit.io/
- **FastAPI Docs:** https://fastapi.tiangolo.com/
- **LangChain Docs:** https://python.langchain.com/
- **Ollama Docs:** https://ollama.ai/
- **Pandas Docs:** https://pandas.pydata.org/

---

## 🚨 Emergency Procedures

### If API Crashes
```bash
# 1. Kill all Python processes
taskkill /F /IM python.exe

# 2. Check for port conflicts
netstat -ano | findstr :8000

# 3. Restart API
python -m uvicorn api:app --host 127.0.0.1 --port 8000 --reload
```

### If Streamlit Freezes
```bash
# 1. Kill Streamlit process
taskkill /F /IM streamlit.exe

# 2. Clear cache
rm -rf ~/.streamlit/cache

# 3. Restart
streamlit run pv_fullstack.py
```

### If Database Locked
```bash
# 1. Check for locked files
lsof | grep mlflow.db

# 2. Kill process holding lock
kill -9 <PID>

# 3. Restart MLflow
mlflow ui --backend-store-uri file:///C:/Users/koreo/mlruns
```

### If Out of Disk Space
```bash
# 1. Check disk usage
du -sh *

# 2. Clean old logs
rm -rf audit_logs/*.log.old

# 3. Archive old reports
tar -czf backups/old_reports.tar.gz sar_reports/reports/*.json

# 4. Remove archived files
rm -rf sar_reports/reports/*.json.old
```

---

## ✅ Verification Checklist

After fixing an error, verify:

- [ ] API is running: `curl http://127.0.0.1:8000/`
- [ ] Streamlit is running: Open http://127.0.0.1:8501
- [ ] No errors in logs: `tail -f audit_logs/audit.log`
- [ ] Required files exist: `ls -la sar_reports/`
- [ ] Ollama is running: `ollama list`
- [ ] Database is accessible: `sqlite3 mlflow.db ".tables"`
- [ ] Ports are available: `netstat -ano | findstr :8000`
- [ ] Permissions are correct: `ls -la gdpr_registry/`

---

**Last Updated:** 2025-12-07  
**Status:** QUICK REFERENCE READY
