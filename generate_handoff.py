"""
PV Signal ML - Project Handoff Generator
Generates complete project state for thread transitions
"""

from pathlib import Path
import json
from datetime import datetime

class ProjectHandoff:
    def __init__(self, project_root="."):
        self.project_root = Path(project_root)
        self.state = {}
        
    def capture_state(self):
        """Capture complete project state"""
        self.state = {
            "metadata": self._get_metadata(),
            "architecture": self._get_architecture(),
            "implementation_status": self._get_implementation_status(),
            "file_structure": self._get_file_structure(),
            "technical_stack": self._get_technical_stack(),
            "recent_changes": self._get_recent_changes(),
            "data_sources": self._get_data_sources(),
            "ml_models": self._get_ml_models(),
            "deployment": self._get_deployment_info(),
            "next_steps": self._get_next_steps(),
            "known_issues": self._get_known_issues()
        }
        return self.state
    
    def _get_metadata(self):
        return {
            "project_name": "PV Signal ML - Pharmacovigilance Signal Detection",
            "version": "2.0 (Enhanced with Causality ML)",
            "last_updated": datetime.now().isoformat(),
            "environment": "Windows 10/11, Python 3.13, Conda base",
            "working_directory": "C:\\Users\\koreo\\PV_Signal_ML"
        }
    
    def _get_architecture(self):
        return {
            "pipeline": [
                {
                    "stage": 1,
                    "name": "Signal Detection",
                    "method": "Statistical (PRR, Chi-square, ROR)",
                    "input": "FAERS adverse event data",
                    "output": "Drug-event signal pairs",
                    "status": "✅ Implemented and working"
                },
                {
                    "stage": 2,
                    "name": "Causality Assessment",
                    "method": "ML (SISA-trained model)",
                    "input": "Detected signals + features",
                    "output": "Causality probability scores",
                    "status": "✅ Integrated, needs model training"
                },
                {
                    "stage": 3,
                    "name": "Explainability",
                    "method": "SHAP values",
                    "input": "ML predictions",
                    "output": "Feature importance",
                    "status": "⚠️ Partially implemented"
                }
            ],
            "key_components": {
                "data_ingestion": "DataSourceManager (Live FAERS + Demo dataset)",
                "signal_detection": "DisproportionalityAnalysis (PRR, Chi², ROR)",
                "causality_ml": "CausalityScorer + SISATrainer (GDPR-compliant)",
                "drug_filtering": "Portfolio-based targeted surveillance",
                "ui": "Streamlit multi-tab dashboard"
            }
        }
    
    def _get_implementation_status(self):
        return {
            "completed": [
                "✅ Live FAERS quarterly data download (2025 Q1)",
                "✅ Statistical signal detection (PRR, Chi², ROR)",
                "✅ Drug portfolio filtering (4 drugs loaded)",
                "✅ SISA model architecture for machine unlearning",
                "✅ CausalityScorer module for ML integration",
                "✅ Multi-tab Streamlit interface",
                "✅ Data source manager (Live/Demo switch)",
                "✅ CSV export functionality"
            ],
            "in_progress": [
                "🔄 SHAP explainability tab (needs model fix)",
                "🔄 Causality ML integration (code ready, needs testing)",
                "🔄 Model training on live FAERS data"
            ],
            "pending": [
                "⏳ RAG integration for evidence generation",
                "⏳ Literature mining (PubMed/Embase)",
                "⏳ GDPR right-to-be-forgotten workflow",
                "⏳ Regulatory report generation (PSMF format)",
                "⏳ MLflow experiment tracking"
            ]
        }
    
    def _get_file_structure(self):
        return {
            "core_files": {
                "app_enhanced.py": "Main Streamlit application (enhanced with ML)",
                "src/data/data_source_manager.py": "Handles Live FAERS and Demo datasets",
                "src/analysis/disproportionality_analysis.py": "Statistical signal detection",
                "src/ml/sisa_trainer.py": "SISA sharding for GDPR compliance",
                "src/ml/causality_scorer.py": "ML causality assessment (NEW)",
                "src/utils/logger.py": "Logging configuration"
            },
            "data_directories": {
                "data/faers_quarterly/": "Downloaded FAERS Q1 2025 data",
                "data/processed/": "Processed datasets",
                "models/": "Trained ML models (sisa_model.pkl)",
                "logs/": "Application and MLflow logs"
            },
            "config_files": {
                "requirements.txt": "Python dependencies",
                "README.md": "Project documentation",
                ".gitignore": "Git exclusions"
            }
        }
    
    def _get_technical_stack(self):
        return {
            "core": {
                "python": "3.13",
                "streamlit": "Latest",
                "pandas": "Latest",
                "numpy": "Latest"
            },
            "ml_frameworks": {
                "scikit-learn": "Model training",
                "xgboost": "Ensemble methods",
                "transformers": "For BERT models (future)",
                "shap": "Model explainability"
            },
            "data_processing": {
                "requests": "API calls and downloads",
                "zipfile": "FAERS data extraction",
                "pickle": "Model serialization"
            },
            "development": {
                "conda": "Environment management",
                "git": "Version control",
                "powershell": "Windows scripting"
            }
        }
    
    def _get_recent_changes(self):
        return {
            "session_summary": [
                "Added drug portfolio filtering (Excel upload)",
                "Integrated SISA machine unlearning architecture",
                "Created CausalityScorer for Stage 2 ML assessment",
                "Fixed column display with conditional causality columns",
                "Added High/Moderate/Low causality metrics",
                "Prepared SHAP explainability (needs debugging)"
            ],
            "key_decisions": [
                "Two-stage pipeline: Statistical → ML Causality",
                "SISA sharding for GDPR right-to-be-forgotten",
                "Conditional column display for gradual feature rollout",
                "Portfolio-based filtering for targeted surveillance"
            ]
        }
    
    def _get_data_sources(self):
        return {
            "live_faers": {
                "source": "FDA FAERS Quarterly Data Files",
                "current_period": "2025Q1 (Jan-Mar 2025)",
                "download_url": "https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html",
                "status": "Downloaded and processed",
                "records": "~50M total, filtered by date range",
                "date_range_used": "2025/01/15 to 2025/04/16"
            },
            "demo_dataset": {
                "source": "HuggingFace (1M records)",
                "purpose": "Quick testing and ML training",
                "status": "Available as fallback"
            },
            "drug_portfolio": {
                "file": "Dummy company drug list.xlsx",
                "drugs": ["ASPIRIN", "IBUPROFEN", "ALLOPURINOL", "METFORMIN"],
                "count": 4
            }
        }
    
    def _get_ml_models(self):
        return {
            "sisa_model": {
                "path": "models/sisa_model.pkl",
                "architecture": "SISA (Sharded, Isolated, Sliced, Aggregated)",
                "purpose": "Causality prediction with machine unlearning",
                "performance": "~0.85 accuracy (reported)",
                "status": "Trained on demo data, needs FAERS retraining",
                "features": ["prr", "chi2", "case_count", "ror"]
            },
            "future_models": {
                "drug_causality_bert": "BERT-based causality from text narratives",
                "ensemble_model": "Combined statistical + NLP features"
            }
        }
    
    def _get_deployment_info(self):
        return {
            "local_dev": {
                "status": "Active",
                "url": "http://localhost:8501",
                "command": "streamlit run app_enhanced.py"
            },
            "planned_deployment": {
                "streamlit_cloud": "For demo/presentation",
                "huggingface_spaces": "Model hosting",
                "github": "Source code repository"
            }
        }
    
    def _get_next_steps(self):
        return {
            "immediate": [
                "1. Test causality ML integration with live FAERS data",
                "2. Fix SHAP explainability tab (model attribute error)",
                "3. Verify High/Moderate/Low causality metrics display",
                "4. Train SISA model on actual FAERS dataset"
            ],
            "short_term": [
                "5. Implement RAG for evidence generation",
                "6. Add literature mining (PubMed API)",
                "7. Create GDPR data deletion workflow",
                "8. Generate regulatory report templates"
            ],
            "long_term": [
                "9. Integrate Drug-Causality-BERT v2",
                "10. Build automated monitoring dashboard",
                "11. Deploy to production environment",
                "12. Add user authentication and audit trails"
            ]
        }
    
    def _get_known_issues(self):
        return {
            "critical": [],
            "high_priority": [
                "SHAP tab error: 'SISATrainer' object has no attribute 'shard_models'",
                "Causality columns not showing (needs model training first)"
            ],
            "medium_priority": [
                "Model training slow on large FAERS dataset",
                "Memory usage high during signal computation"
            ],
            "low_priority": [
                "UI styling improvements needed",
                "Add more comprehensive logging"
            ]
        }
    
    def generate_handoff_prompt(self):
        """Generate complete handoff prompt for new thread"""
        state = self.capture_state()
        
        prompt = f"""# PV Signal ML Project Handoff

## Project Context
I'm continuing development of **{state['metadata']['project_name']}**, a production-grade pharmacovigilance signal detection system with ML-based causality assessment.

**Last Update:** {state['metadata']['last_updated']}
**Working Directory:** `{state['metadata']['working_directory']}`
**Environment:** {state['metadata']['environment']}

---

## Architecture Overview

### Two-Stage Pipeline
"""
        
        for stage in state['architecture']['pipeline']:
            prompt += f"""
**Stage {stage['stage']}: {stage['name']}** {stage['status']}
- Method: {stage['method']}
- Input: {stage['input']}
- Output: {stage['output']}
"""
        
        prompt += f"""
### Key Components
"""
        for comp, desc in state['architecture']['key_components'].items():
            prompt += f"- **{comp}**: {desc}\n"
        
        prompt += f"""
---

## Current Implementation Status

### ✅ Completed
"""
        for item in state['implementation_status']['completed']:
            prompt += f"{item}\n"
        
        prompt += f"""
### 🔄 In Progress
"""
        for item in state['implementation_status']['in_progress']:
            prompt += f"{item}\n"
        
        prompt += f"""
### ⏳ Pending
"""
        for item in state['implementation_status']['pending']:
            prompt += f"{item}\n"
        
        prompt += f"""
---

## File Structure

### Core Application Files
"""
        for file, desc in state['file_structure']['core_files'].items():
            prompt += f"- `{file}`: {desc}\n"
        
        prompt += f"""
### Data & Models
"""
        for dir, desc in state['file_structure']['data_directories'].items():
            prompt += f"- `{dir}`: {desc}\n"
        
        prompt += f"""
---

## Technical Stack

**Core:** Python {state['technical_stack']['core']['python']}, Streamlit, Pandas, NumPy
**ML:** scikit-learn, XGBoost, SHAP, (Transformers for future BERT integration)
**Development:** Conda, Git, PowerShell

---

## Data Sources

### Live FAERS Data
- **Period:** {state['data_sources']['live_faers']['current_period']}
- **Status:** {state['data_sources']['live_faers']['status']}
- **Date Range:** {state['data_sources']['live_faers']['date_range_used']}
- **Records:** {state['data_sources']['live_faers']['records']}

### Drug Portfolio
- **Drugs:** {', '.join(state['data_sources']['drug_portfolio']['drugs'])}
- **Count:** {state['data_sources']['drug_portfolio']['count']} drugs loaded

---

## ML Models

### SISA Model (Causality Prediction)
- **Path:** `{state['ml_models']['sisa_model']['path']}`
- **Architecture:** {state['ml_models']['sisa_model']['architecture']}
- **Performance:** {state['ml_models']['sisa_model']['performance']}
- **Status:** {state['ml_models']['sisa_model']['status']}
- **Features:** {', '.join(state['ml_models']['sisa_model']['features'])}

**Purpose:** GDPR-compliant machine unlearning for "right to be forgotten" requests

---

## Recent Session Summary

### Changes Made
"""
        for change in state['recent_changes']['session_summary']:
            prompt += f"- {change}\n"
        
        prompt += f"""
### Key Decisions
"""
        for decision in state['recent_changes']['key_decisions']:
            prompt += f"- {decision}\n"
        
        prompt += f"""
---

## Known Issues

### High Priority
"""
        for issue in state['known_issues']['high_priority']:
            prompt += f"- {issue}\n"
        
        prompt += f"""
### Medium Priority
"""
        for issue in state['known_issues']['medium_priority']:
            prompt += f"- {issue}\n"
        
        prompt += f"""
---

## Next Steps (Priority Order)

### Immediate Tasks
"""
        for task in state['next_steps']['immediate']:
            prompt += f"{task}\n"
        
        prompt += f"""
### Short-Term Goals
"""
        for task in state['next_steps']['short_term']:
            prompt += f"{task}\n"
        
        prompt += f"""
### Long-Term Vision
"""
        for task in state['next_steps']['long_term']:
            prompt += f"{task}\n"
        
        prompt += f"""
---

## Where I Need Help

**Current Focus:** Testing the complete two-stage pipeline (Statistical Signal Detection → ML Causality Assessment) with live FAERS data.

**Specific Questions:**
1. How to verify causality scoring is working correctly?
2. Best approach to train SISA model on large FAERS dataset efficiently?
3. Fixing SHAP explainability tab attribute error

**Context:** The app is running locally, causality integration code is in place, but needs testing after restart.

---

## Quick Start Commands

Project directory
cd C:\Users\koreo\PV_Signal_ML

Activate environment
conda activate base

Run application
streamlit run app_enhanced.py

Key files to check
- app_enhanced.py (line 82-105: causality integration)
- src/ml/causality_scorer.py (new module)
- models/sisa_model.pkl (trained model)


---

## Additional Context

- **Regulatory Compliance:** FDA/EMA pharmacovigilance standards, GDPR for data privacy
- **Portfolio Drugs:** Targeted surveillance for company-specific drugs
- **Signal Detection:** Statistical methods (PRR ≥ 2, Chi² ≥ 4) are regulatory standard
- **ML Enhancement:** Causality assessment prioritizes signals for manual review
- **Machine Unlearning:** SISA sharding allows deletion of specific cases without full retraining

Ready to continue from here! Please help with testing the causality integration.
"""
        
        return prompt
    
    def save_state(self, output_file="project_handoff.json"):
        """Save state as JSON"""
        state = self.capture_state()
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(state, f, indent=2)
        return output_file
    
    def save_prompt(self, output_file="handoff_prompt.md"):
        """Save handoff prompt as Markdown"""
        prompt = self.generate_handoff_prompt()
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(prompt)
        return output_file

# CLI usage
if __name__ == "__main__":
    handoff = ProjectHandoff()
    
    print("=" * 70)
    print("PV Signal ML - Project Handoff Generator")
    print("=" * 70)
    print()
    
    print("Capturing project state...")
    state = handoff.capture_state()
    
    print("Generating handoff documents...")
    json_file = handoff.save_state()
    prompt_file = handoff.save_prompt()
    
    print(f"✅ State saved to: {json_file}")
    print(f"✅ Handoff prompt saved to: {prompt_file}")
    print()
    print("📋 To use in new thread:")
    print(f"   1. Copy content from {prompt_file}")
    print("   2. Paste as first message in new conversation")
    print("   3. Continue development seamlessly!")
    print()
    print("=" * 70)
