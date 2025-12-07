# Comprehensive Testing Guide

**Date:** 2025-12-08  
**Status:** Ready for Testing  
**Scope:** All features including auto-start API, PSMF Annex Generator, and MLflow audit logging

---

## 🎯 Testing Objectives

1. ✅ Verify API auto-starts when Streamlit app loads
2. ✅ Verify SAR generation works without manual API startup
3. ✅ Verify PSMF Annex D generator creates proper documentation
4. ✅ Verify MLflow runs are captured in audit trail
5. ✅ Verify no errors occur during report generation
6. ✅ Verify markdown files are created successfully
7. ✅ Verify all file paths are Windows-compatible

---

## 📋 Pre-Testing Checklist

### Environment Setup
- [ ] Python 3.13 installed
- [ ] All dependencies installed: `pip install -r requirements.txt`
- [ ] Ollama running: `ollama serve` (in separate terminal)
- [ ] Model pulled: `ollama pull llama3.2:3b`
- [ ] Port 8000 available (check: `netstat -ano | findstr :8000`)
- [ ] Port 8501 available (check: `netstat -ano | findstr :8501`)

### File Checks
- [ ] `pv_fullstack.py` exists and is updated
- [ ] `psmf_annex_generator.py` exists
- [ ] `api.py` exists and is updated
- [ ] `audit_logging.py` exists and is updated
- [ ] `sar_reports/full_signals_1M.csv` exists
- [ ] `templates/signal_report_template.md` exists

### Directory Structure
```
pv-signal-ml/
├── pv_fullstack.py ✅
├── psmf_annex_generator.py ✅
├── api.py ✅
├── audit_logging.py ✅
├── sar_reports/
│   ├── full_signals_1M.csv ✅
│   ├── reports/ (will be created)
│   └── model_run_metadata.json
├── templates/
│   └── signal_report_template.md ✅
├── audit_logs/ (will be created)
└── psmf_annexes/ (will be created)
```

---

## 🧪 Test Cases

### Test 1: Streamlit App Startup
**Objective:** Verify app loads without errors

**Steps:**
1. Open terminal in project directory
2. Run: `streamlit run pv_fullstack.py`
3. Wait for "You can now view your Streamlit app in your browser"
4. Open http://127.0.0.1:8501

**Expected Results:**
- ✅ App loads without errors
- ✅ "Starting API server..." message appears
- ✅ "✅ API server started successfully!" message appears
- ✅ Dashboard tab is visible
- ✅ PSMF Annex Generator tab is visible
- ✅ Methods & Governance tab is visible

**Pass/Fail:** ___________

---

### Test 2: API Auto-Start
**Objective:** Verify API starts automatically without manual intervention

**Steps:**
1. Check if API is already running: `netstat -ano | findstr :8000`
2. If running, kill it: `taskkill /PID <PID> /F`
3. Start Streamlit app: `streamlit run pv_fullstack.py`
4. Wait for API startup message
5. Verify API is running: `curl http://127.0.0.1:8000/`

**Expected Results:**
- ✅ API starts automatically
- ✅ No manual terminal needed
- ✅ API responds to health check
- ✅ Success message appears in Streamlit

**Pass/Fail:** ___________

---

### Test 3: SAR Generation (Dashboard Tab)
**Objective:** Verify signal assessment report generation works

**Steps:**
1. Go to Dashboard tab
2. Select signal: "OncoKill (Cisplatin)" → "Acute Kidney Injury"
3. Set date range: 2024-01-01 to 2024-03-31
4. Click "Generate signal assessment report"
5. Wait for report generation
6. Check for success message

**Expected Results:**
- ✅ No WinError 123 (filename syntax error)
- ✅ Success message: "Report generated and stored under..."
- ✅ JSON evidence displayed
- ✅ Markdown report displayed
- ✅ Download button available
- ✅ Files created in `sar_reports/reports/`

**Files to Check:**
```bash
ls -la sar_reports/reports/OncoKill*
# Should see:
# OncoKill__Cisplatin____Acute_Kidney_Injury__2024-01-01_2024-03-31.json
# OncoKill__Cisplatin____Acute_Kidney_Injury__2024-01-01_2024-03-31.md
```

**Pass/Fail:** ___________

---

### Test 4: PSMF Annex D Generation
**Objective:** Verify PSMF Annex D generator creates proper documentation

**Steps:**
1. Go to "PSMF Annex Generator" tab
2. Read the description
3. Click "🔄 Generate PSMF Annex D"
4. Wait for generation
5. Check for success message

**Expected Results:**
- ✅ Success message appears
- ✅ JSON metadata displayed (File, Products, Generated At, Size)
- ✅ Document preview shown (first 2000 characters)
- ✅ Download button available
- ✅ File created in `psmf_annexes/`

**Files to Check:**
```bash
ls -la psmf_annexes/
# Should see:
# PSMF_Annex_D_Signal_Management_*.md
```

**Content to Verify:**
- [ ] Title: "PSMF Annex D: Signal Management Procedures"
- [ ] Sections: Introduction, Methodology, System Description, Products, Workflow, Quality Assurance
- [ ] Contains: PRR formula, thresholds, regulatory references
- [ ] Contains: System architecture diagram
- [ ] Contains: Products list
- [ ] Contains: Compliance mapping

**Pass/Fail:** ___________

---

### Test 5: Markdown File Creation
**Objective:** Verify markdown files are created without errors

**Steps:**
1. Generate SAR (Test 3)
2. Check if markdown file exists
3. Open file in text editor
4. Verify content is readable

**Expected Results:**
- ✅ Markdown file created
- ✅ File is readable (no encoding errors)
- ✅ Contains proper markdown formatting
- ✅ No "Markdown file not found on disk" error

**Files to Check:**
```bash
cat sar_reports/reports/OncoKill__Cisplatin____Acute_Kidney_Injury__2024-01-01_2024-03-31.md
```

**Pass/Fail:** ___________

---

### Test 6: Audit Logging
**Objective:** Verify all operations are logged to audit trail

**Steps:**
1. Generate SAR (Test 3)
2. Generate PSMF Annex D (Test 4)
3. Check audit logs

**Expected Results:**
- ✅ `audit_logs/access_log.jsonl` created
- ✅ `audit_logs/report_generation.jsonl` created
- ✅ `audit_logs/mlflow_runs.jsonl` created
- ✅ `audit_logs/audit.log` created
- ✅ All operations logged with timestamp

**Files to Check:**
```bash
# Check audit logs exist
ls -la audit_logs/

# View recent entries
tail -10 audit_logs/audit.log
tail -5 audit_logs/report_generation.jsonl
tail -5 audit_logs/mlflow_runs.jsonl
```

**Pass/Fail:** ___________

---

### Test 7: MLflow Integration
**Objective:** Verify MLflow runs are captured

**Steps:**
1. Check MLflow directory: `C:/Users/koreo/mlruns`
2. View MLflow UI: `mlflow ui`
3. Check for runs in experiment

**Expected Results:**
- ✅ MLflow directory exists
- ✅ Runs are logged
- ✅ Experiments are visible in UI
- ✅ Metrics and parameters are recorded

**Commands:**
```bash
# View MLflow runs
ls -la C:/Users/koreo/mlruns/

# Start MLflow UI
mlflow ui --backend-store-uri file:///C:/Users/koreo/mlruns

# Then open http://127.0.0.1:5000
```

**Pass/Fail:** ___________

---

### Test 8: Special Characters Handling
**Objective:** Verify special characters in drug/event names are handled correctly

**Steps:**
1. Try generating report for drug with special characters
2. Example: "Drug/Name" or "Event:Name"
3. Verify no filename errors

**Expected Results:**
- ✅ No WinError 123
- ✅ Filenames use safe characters (underscores)
- ✅ Report generated successfully

**Pass/Fail:** ___________

---

### Test 9: Error Handling
**Objective:** Verify errors are handled gracefully

**Steps:**
1. Try generating report with missing CSV file
2. Try generating PSMF with missing signals database
3. Check error messages

**Expected Results:**
- ✅ Graceful error messages
- ✅ No crashes
- ✅ Helpful troubleshooting tips
- ✅ Errors logged to audit trail

**Pass/Fail:** ___________

---

### Test 10: Download Functionality
**Objective:** Verify download buttons work

**Steps:**
1. Generate SAR
2. Click "Download report (.md)"
3. Check downloaded file

**Expected Results:**
- ✅ File downloads successfully
- ✅ Filename is correct
- ✅ Content is readable
- ✅ No encoding issues

**Pass/Fail:** ___________

---

## 🔍 Verification Checklist

### Code Quality
- [ ] No syntax errors in pv_fullstack.py
- [ ] No syntax errors in psmf_annex_generator.py
- [ ] No syntax errors in api.py
- [ ] No syntax errors in audit_logging.py
- [ ] All imports are available
- [ ] No undefined variables

### Functionality
- [ ] API auto-starts
- [ ] SAR generation works
- [ ] PSMF generation works
- [ ] Markdown files created
- [ ] JSON files created
- [ ] Download buttons work
- [ ] Error messages are helpful

### Compliance
- [ ] Audit logs created
- [ ] MLflow runs logged
- [ ] No PII in logs
- [ ] Timestamps correct
- [ ] File paths Windows-compatible

### Performance
- [ ] App loads in <5 seconds
- [ ] API starts in <30 seconds
- [ ] SAR generation in <60 seconds
- [ ] PSMF generation in <30 seconds
- [ ] No memory leaks

---

## 📊 Test Results Summary

| Test # | Test Name | Status | Notes |
|---|---|---|---|
| 1 | Streamlit Startup | _____ | |
| 2 | API Auto-Start | _____ | |
| 3 | SAR Generation | _____ | |
| 4 | PSMF Generation | _____ | |
| 5 | Markdown Creation | _____ | |
| 6 | Audit Logging | _____ | |
| 7 | MLflow Integration | _____ | |
| 8 | Special Characters | _____ | |
| 9 | Error Handling | _____ | |
| 10 | Download Functionality | _____ | |

**Overall Status:** _______________

---

## 🐛 Bug Report Template

If you find issues, please document them:

```
**Bug Title:** [Brief description]

**Severity:** Critical / High / Medium / Low

**Steps to Reproduce:**
1. [Step 1]
2. [Step 2]
3. [Step 3]

**Expected Result:**
[What should happen]

**Actual Result:**
[What actually happened]

**Error Message:**
[Full error message if applicable]

**Screenshots:**
[Attach if helpful]

**Environment:**
- OS: Windows 10/11
- Python: 3.13
- Streamlit: [version]
- FastAPI: [version]
```

---

## ✅ Sign-Off

**Tested By:** ___________________  
**Date:** ___________________  
**Status:** ___________________  
**Ready for Deployment:** Yes / No

---

**Next Steps:**
1. Complete all tests
2. Document any issues
3. Fix critical bugs
4. Re-test if needed
5. Proceed to GitHub deployment

---

*This testing guide ensures comprehensive validation before production deployment.*
