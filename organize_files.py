"""
File Organization Script for PV-Signal-ML

Moves experimental and iteration files to Experimental/ folder
to keep production code clean for Git commits and Streamlit deployment.
"""

import shutil
from pathlib import Path
import json
from datetime import datetime


def organize_project_files():
    """
    Organize project files into production and experimental categories.
    """
    base_dir = Path.cwd()
    experimental_dir = base_dir / "Experimental"
    experimental_dir.mkdir(exist_ok=True)
    
    # Files to move to Experimental/
    experimental_files = [
        # UI iterations
        "pv_ui_complete_ORIGINAL.py",
        "pv_ui_complete.py",
        "pv_ui_complete_enhanced.py",
        "pv_ui_minimal_fixed.py",
        "pv_ui_production.py",
        
        # RAG iterations
        "rag_langchain_fixed.py",
        
        # SHAP iterations
        "shap_analysis.py",
        "shap_analysis_fixed.py",
        
        # FAERS period-specific
        "faers_build_signals_q1_2025.py",
        
        # Notebooks
        "faers_ingestion_exploration.ipynb",
        "faers_q1_2024_explore.ipynb",
        
        # PowerShell scripts (legacy)
        "add_compliance_tab.ps1",
        "cioms_compliance_update.ps1",
        "launch_production_ui.ps1",
        "pipeline_langchain_update.ps1",
        "run_complete_regulatory_system.ps1",
        "run_complete_regulatory_system_ORIGINAL.ps1",
        "run_full_pipeline.ps1",
        "run_full_regulatory_pipeline.ps1",
        "ui_sar_factory_tab.ps1",
        
        # Utility scripts
        "debug.py",
        "check_mlflow.py",
        "load_1M_sqlite.py",
        "add_rag_to_ui.py",
        "api_signal_reports.py",
    ]
    
    # Files to keep in root (production)
    production_files = [
        # Main apps
        "pv_fullstack.py",
        "pv_ui.py",
        "api.py",
        
        # Core pipeline
        "faers_build_signals.py",
        "stats_engine.py",
        "prepare_ml_features.py",
        
        # ML
        "pv_signal_ml_pipeline.py",
        "shap_analysis_simple.py",
        "train_with_mlflow.py",
        
        # RAG
        "rag_langchain.py",
        "rag_signal_evidence.py",
        "rag_engine.py",
        "rag_pv_signals.py",
        
        # Reports
        "signal_report_builder.py",
        "generate_psmf.py",
        "export_assessment_bundle.py",
        
        # Governance & compliance
        "data_lineage.py",
        "gdpr_deletion_registry.py",
        "audit_logging.py",
        "governance_dpia.md",
        "change_control.py",
        "drift_monitor.py",
        "fairness_analyzer.py",
        "run_metadata.py",
        
        # Documentation
        "README.md",
        "ANALYSIS_AND_COMPLIANCE_REPORT.md",
        "PSMF_v1.0.md",
        
        # Config
        "requirements.txt",
    ]
    
    moved_files = []
    skipped_files = []
    errors = []
    
    print("🔄 Organizing project files...\n")
    
    for filename in experimental_files:
        file_path = base_dir / filename
        
        if file_path.exists():
            try:
                dest_path = experimental_dir / filename
                shutil.move(str(file_path), str(dest_path))
                moved_files.append(filename)
                print(f"✅ Moved: {filename} → Experimental/")
            except Exception as e:
                errors.append((filename, str(e)))
                print(f"❌ Error moving {filename}: {e}")
        else:
            skipped_files.append(filename)
    
    print(f"\n📊 Summary:")
    print(f"  ✅ Moved: {len(moved_files)} files")
    print(f"  ⏭️  Skipped (not found): {len(skipped_files)} files")
    print(f"  ❌ Errors: {len(errors)} files")
    
    if errors:
        print(f"\n⚠️  Errors encountered:")
        for filename, error in errors:
            print(f"  - {filename}: {error}")
    
    # Create organization report
    report = {
        "timestamp": datetime.utcnow().isoformat() + "Z",
        "action": "file_organization",
        "moved_files": moved_files,
        "skipped_files": skipped_files,
        "errors": [{"file": f, "error": e} for f, e in errors],
        "production_files": production_files,
        "summary": {
            "total_moved": len(moved_files),
            "total_skipped": len(skipped_files),
            "total_errors": len(errors)
        }
    }
    
    # Save report
    report_path = base_dir / "file_organization_report.json"
    with open(report_path, 'w', encoding='utf-8') as f:
        json.dump(report, f, indent=2, default=str)
    
    print(f"\n📋 Report saved: {report_path}")
    
    # Create .gitignore entry for Experimental folder
    gitignore_path = base_dir / ".gitignore"
    gitignore_entry = "\n# Experimental files (iterations and prototypes)\nExperimental/\n"
    
    if gitignore_path.exists():
        with open(gitignore_path, 'r', encoding='utf-8') as f:
            content = f.read()
        if "Experimental/" not in content:
            with open(gitignore_path, 'a', encoding='utf-8') as f:
                f.write(gitignore_entry)
            print(f"✅ Updated .gitignore to exclude Experimental/")
    else:
        with open(gitignore_path, 'w', encoding='utf-8') as f:
            f.write(gitignore_entry)
        print(f"✅ Created .gitignore with Experimental/ exclusion")
    
    print(f"\n✅ File organization complete!")
    print(f"\n📁 Production files ready for Git commit:")
    for f in sorted(production_files):
        if (base_dir / f).exists():
            print(f"  - {f}")
    
    return report


def create_experimental_readme():
    """
    Create README for Experimental folder explaining its purpose.
    """
    experimental_dir = Path.cwd() / "Experimental"
    experimental_dir.mkdir(exist_ok=True)
    
    readme_content = """# Experimental Files

This folder contains iteration files, prototypes, and experimental code that are not part of the production system.

## Contents

- **UI Iterations:** Earlier versions of Streamlit apps (pv_ui_complete.py, pv_ui_complete_enhanced.py, etc.)
- **RAG Iterations:** Earlier versions of RAG pipeline (rag_langchain_fixed.py)
- **SHAP Iterations:** Earlier versions of explainability code
- **Notebooks:** Jupyter notebooks for exploration and analysis
- **Legacy Scripts:** PowerShell scripts and utility scripts from earlier development phases

## Why These Files Are Here

These files represent the development process and are kept for reference, but they are not used in the production system. They are excluded from Git commits to keep the repository clean.

## Using Production Files Instead

For production use, refer to the main files in the root directory:

- **UI:** `pv_fullstack.py` or `pv_ui.py`
- **RAG:** `rag_langchain.py`
- **SHAP:** `shap_analysis_simple.py`
- **FAERS:** `faers_build_signals.py`

## Cleanup

If you want to permanently remove these files, you can safely delete the entire `Experimental/` folder.

---

**Last Updated:** 2025-12-07
"""
    
    readme_path = experimental_dir / "README.md"
    with open(readme_path, 'w', encoding='utf-8') as f:
        f.write(readme_content)
    
    print(f"✅ Created Experimental/README.md")


if __name__ == "__main__":
    print("=" * 70)
    print("PV-Signal-ML File Organization Tool")
    print("=" * 70)
    print()
    
    # Organize files
    report = organize_project_files()
    
    # Create Experimental README
    print()
    create_experimental_readme()
    
    print("\n" + "=" * 70)
    print("✅ File organization complete!")
    print("=" * 70)
    print("\n📝 Next steps:")
    print("  1. Review the file_organization_report.json for details")
    print("  2. Commit changes to Git: git add -A && git commit -m 'Organize files'")
    print("  3. Deploy to Streamlit: streamlit run pv_fullstack.py")
    print()
