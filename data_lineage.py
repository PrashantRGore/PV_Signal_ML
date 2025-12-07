from pathlib import Path
import json
from datetime import datetime
import hashlib
import pandas as pd

class DataLineageTracker:
    def __init__(self, base_dir: Path):
        self.base_dir = base_dir
        self.lineage_dir = base_dir / 'lineage'
        self.lineage_dir.mkdir(exist_ok=True)
    
    def record_dataset(self, dataset_path: Path, source_url: str, script_version: str, description: str):
        """Record full data provenance for regulatory audit"""
        lineage = {
            'dataset': dataset_path.name,
            'source': source_url,
            'extraction_date': datetime.now().isoformat(),
            'script_version': script_version,
            'description': description,
            'file_hash': self._compute_hash(dataset_path),
            'row_count': self._get_row_count(dataset_path),
            'column_count': self._get_column_count(dataset_path)
        }
        
        lineage_path = self.lineage_dir / f"{dataset_path.stem}_lineage.json"
        with open(lineage_path, 'w') as f:
            json.dump(lineage, f, indent=2, default=str)
        print(f"✅ Lineage recorded: {lineage_path}")
        return lineage_path
    
    def _compute_hash(self, file_path: Path) -> str:
        sha256 = hashlib.sha256()
        with open(file_path, 'rb') as f:
            for chunk in iter(lambda: f.read(4096), b""):
                sha256.update(chunk)
        return sha256.hexdigest()
    
    def _get_row_count(self, file_path: Path) -> int:
        if file_path.suffix == '.csv':
            return len(pd.read_csv(file_path))
        return 0
    
    def _get_column_count(self, file_path: Path) -> int:
        if file_path.suffix == '.csv':
            return len(pd.read_csv(file_path, nrows=1).columns)
        return 0

def create_governance_docs():
    """Create regulatory governance documentation"""
    governance = {
        "purpose": "Pharmacovigilance signal detection digital twin",
        "legal_basis": "Public health interest (Art. 9(2)(i) GDPR)",
        "data_categories": ["Aggregated case counts", "Pseudonymised drug-event pairs"],
        "privacy_measures": [
            "No patient identifiers retained",
            "Aggregated counts only (A,B,C,D 2x2 tables)",
            "No free text/narratives processed",
            "FAERS public domain data"
        ],
        "retention_policy": "FAERS quarterly data + 10 years post-product lifecycle",
        "human_in_loop": "ML triage supports stats engine; humans validate signals",
        "audit_trail": "MLflow + data lineage JSONs + Git commit history"
    }
    
    with open('governance_dpia.md', 'w') as f:
        f.write("# PV Signal Detection DPIA / Governance\n\n")
        f.write(f"**Generated:** {datetime.now().isoformat()}\n\n")
        for k, v in governance.items():
            f.write(f"## {k.replace('_', ' ').title()}\n")
            if isinstance(v, list):
                for item in v:
                    f.write(f"- {item}\n")
            else:
                f.write(f"{v}\n\n")
    print("✅ governance_dpia.md created")

if __name__ == "__main__":
    tracker = DataLineageTracker(Path.cwd())
    tracker.record_dataset(
        Path('sar_reports/candidate_signals_2025Q1_stats_engine.csv'),
        "https://fis.fda.gov/sense/app/9524532e-2eb4-490e-b914-0a5f8e970e2d/sheet/7a5acf3b-72d4-4b5d-99a7-4ca3e9ca74d4/state/analysis",
        "v1.0.0",
        "FAERS 2025Q1 candidate signals (PRR>=2, CHISQ>=4, CASES>=3)"
    )
    create_governance_docs()
