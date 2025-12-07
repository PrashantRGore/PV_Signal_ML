import pandas as pd
from pathlib import Path
import ollama
from sentence_transformers import SentenceTransformer
import numpy as np
from typing import List, Dict

class PVSignalRAG:
    def __init__(self):
        self.base_dir = Path.cwd()
        self.model = SentenceTransformer('all-MiniLM-L6-v2')
        self.context_docs = self._load_context_docs()
        self.embeddings = self.model.encode(self.context_docs)
    
    def _load_context_docs(self) -> List[str]:
        """Load regulatory guidance + SAR context for RAG"""
        docs = []
        
        # Load existing SAR bundles as context
        sar_dir = self.base_dir / 'sar_reports' / 'reports'
        for report in sar_dir.glob('*.md'):
            with open(report, 'r', encoding='utf-8') as f:
                docs.append(f"Signal Report: {report.stem}\n" + f.read())
        
        # Regulatory guidance context
        regulatory = [
            "EMA Guideline: PRR ≥ 2, Chi² ≥ 4, min 3 cases indicates potential signal requiring review [web:21]",
            "FDA FAERS: Signals require clinical review + literature search + causality assessment",
            "CIOMS VIII: Statistical signals must be biologically plausible with supporting evidence",
            "Signal triage prioritizes: seriousness, novelty, biological plausibility, case quality"
        ]
        docs.extend(regulatory)
        
        return docs
    
    def explain_signal(self, drug: str, event: str, prr: float, chi2: float, cases: int) -> Dict:
        """Generate regulatory-grade signal explanation via RAG"""
        query = f"{drug} {event} PRR={prr:.1f} Chi2={chi2:.1f} cases={cases}"
        query_emb = self.model.encode([query])
        
        # Semantic search
        scores = np.dot(query_emb, self.embeddings.T)[0]
        top_idx = np.argmax(scores)
        
        explanation = {
            'drug': drug,
            'event': event,
            'prr': prr,
            'chi2': chi2,
            'cases': cases,
            'signal_strength': 'HIGH' if prr > 100 else 'MEDIUM' if prr > 10 else 'LOW',
            'rag_context': self.context_docs[top_idx][:500] + '...',
            'recommendation': self._generate_recommendation(prr, chi2, cases),
            'regulatory_flags': self._check_regulatory_criteria(prr, chi2, cases)
        }
        return explanation
    
    def _generate_recommendation(self, prr: float, chi2: float, cases: int) -> str:
        if prr > 1000 and cases >= 10:
            return "URGENT REVIEW: Extreme disproportionality + sufficient cases"
        elif prr > 10 and cases >= 5:
            return "PRIORITIZED REVIEW: Strong statistical signal"
        elif prr >= 2 and cases >= 3:
            return "STANDARD REVIEW: Meets basic signal criteria"
        return "MONITOR: Weak statistical signal"
    
    def _check_regulatory_criteria(self, prr: float, chi2: float, cases: int) -> List[str]:
        flags = []
        if prr >= 2: flags.append("EMA PRR threshold met")
        if chi2 >= 4: flags.append("Chi-square threshold met") 
        if cases >= 3: flags.append("Minimum cases met")
        if prr > 100: flags.append("HIGH PRIORITY signal")
        return flags

if __name__ == "__main__":
    rag = PVSignalRAG()
    # Test with top signal
    explanation = rag.explain_signal(
        "CALCIUM CHLORIDE\\MAGNESIUM CHLORIDE\\POTASSIUM CHLORIDE",
        "Keratopathy", 1e9, 56601, 15
    )
    print("RAG Signal Explanation:")
    import json
    print(json.dumps(explanation, indent=2))
