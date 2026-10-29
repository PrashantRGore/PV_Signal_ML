"""
RAG Orchestrator - Direct Ollama API (No LangChain)
Graph RAG + Ollama for regulatory SAR generation
"""

from typing import Dict, List, Optional
import json
from pathlib import Path
import pandas as pd
import requests
from src.rag.graph.knowledge_graph_builder import RegulatoryKnowledgeGraph
from src.utils.logger import setup_logger

logger = setup_logger(__name__)


class GraphRAGOrchestrator:
    """Graph RAG with direct Ollama API integration."""
    
    def __init__(
        self,
        ollama_model: str = "llama3.2",
        ollama_url: str = "http://localhost:11434",
        knowledge_graph_path: Optional[Path] = None
    ):
        self.ollama_model = ollama_model
        self.ollama_url = ollama_url
        
        # Load knowledge graph
        self.kg = RegulatoryKnowledgeGraph()
        if knowledge_graph_path and knowledge_graph_path.exists():
            self._load_knowledge_graph(knowledge_graph_path)
        
        logger.info(f"✅ RAG Orchestrator initialized with {ollama_model}")
    
    def _load_knowledge_graph(self, kg_path: Path):
        """Load pre-built knowledge graph."""
        try:
            with open(kg_path, 'r') as f:
                kg_data = json.load(f)
            self.kg.nodes = kg_data['nodes']
            self.kg.edges = kg_data['edges']
            logger.info(f"✅ Loaded knowledge graph: {kg_data['metadata']['node_count']} nodes")
        except Exception as e:
            logger.warning(f"Could not load knowledge graph: {e}")
    
    def generate_sar(self, signal_data: Dict) -> Dict[str, str]:
        """
        Generate SAR using Graph RAG + Ollama.
        
        Args:
            signal_data: Dict with drug_name, event_name, count, prr, chi_square
            
        Returns:
            Dict with 'report', 'evidence_sources', 'model_used'
        """
        logger.info(f"Generating SAR for {signal_data.get('drug_name')} - {signal_data.get('event_name')}")
        
        # Retrieve context from knowledge graph
        drug_context = self.kg.get_drug_context(signal_data.get('drug_name', ''))
        similar_signals = self.kg.query_similar_signals(
            signal_data.get('drug_name', ''),
            signal_data.get('event_name', ''),
            top_k=3
        )
        
        # Build prompt
        prompt = self._build_prompt(signal_data, similar_signals, drug_context)
        
        # Call Ollama API
        try:
            response = requests.post(
                f"{self.ollama_url}/api/generate",
                json={
                    "model": self.ollama_model,
                    "prompt": prompt,
                    "stream": False,
                    "options": {"temperature": 0.3}
                },
                timeout=120
            )
            
            if response.status_code == 200:
                sar_report = response.json().get('response', '')
                logger.info("✅ SAR generated successfully")
                
                return {
                    'report': sar_report,
                    'evidence_sources': self._extract_evidence_sources(drug_context, similar_signals),
                    'graph_context': drug_context,
                    'model_used': self.ollama_model,
                    'success': True
                }
            else:
                error_msg = f"Ollama API error: {response.status_code}"
                logger.error(error_msg)
                return {
                    'report': error_msg,
                    'error': True,
                    'success': False
                }
                
        except requests.exceptions.ConnectionError:
            error_msg = "Cannot connect to Ollama. Is it running? Start with: ollama serve"
            logger.error(error_msg)
            return {
                'report': error_msg,
                'error': True,
                'success': False
            }
        except Exception as e:
            error_msg = f"SAR generation failed: {str(e)}"
            logger.error(error_msg)
            return {
                'report': error_msg,
                'error': True,
                'success': False
            }
    
    def _build_prompt(self, signal_data: Dict, similar_signals: List, drug_context: Dict) -> str:
        """Build SAR generation prompt."""
        similar_text = self._format_similar_signals(similar_signals)
        clinical_text = self._extract_clinical_context(signal_data)
        guidelines_text = self._get_guidelines_text()
        
        prompt = f"""You are a pharmacovigilance expert generating a Signal Assessment Report (SAR) for regulatory submission.

**SIGNAL INFORMATION:**
Drug: {signal_data.get('drug_name', 'Unknown')}
Event: {signal_data.get('event_name', 'Unknown')}
Case Count: {signal_data.get('count', signal_data.get('case_count', 'Unknown'))}
PRR (Proportional Reporting Ratio): {signal_data.get('prr', 0):.2f}
Chi-Square: {signal_data.get('chi_square', signal_data.get('chi2', 0)):.2f}

**SIMILAR KNOWN SIGNALS FROM KNOWLEDGE GRAPH:**
{similar_text}

**REGULATORY GUIDELINES:**
{guidelines_text}

**CLINICAL CONTEXT:**
{clinical_text}

**TASK:**
Generate a comprehensive Signal Assessment Report following CIOMS XIV and EMA GVP Module IX guidelines.

Include these sections:
1. **Signal Summary** (2-3 sentences)
2. **Statistical Evidence** (interpret PRR, Chi², case count)
3. **Clinical Plausibility** (mechanism, biological basis)
4. **Similar Associations** (reference knowledge graph)
5. **Regulatory Considerations** (cite FDA, EMA, CIOMS)
6. **Causality Assessment** (WHO-UMC criteria)
7. **Recommendation** (risk assessment, next steps)

Use factual language. Cite guidelines. Do NOT invent details. Flag missing data explicitly.

Generate the SAR report:"""

        return prompt
    
    def _format_similar_signals(self, similar_signals: List[Dict]) -> str:
        """Format similar signals for prompt."""
        if not similar_signals:
            return "No similar signals found in knowledge graph."
        
        formatted = []
        for i, signal in enumerate(similar_signals, 1):
            formatted.append(
                f"{i}. {signal['drug']} → {signal['event']} "
                f"(PRR: {signal['prr']:.2f}, Similarity: {signal['similarity_reason']})"
            )
        
        return "\n".join(formatted)
    
    def _get_guidelines_text(self) -> str:
        """Get regulatory guidelines text."""
        guidelines = self.kg.nodes.get('guidelines', [])
        
        if not guidelines:
            return "Standard pharmacovigilance guidelines apply (FDA, EMA GVP Module IX, CIOMS XIV)."
        
        formatted = []
        for guideline in guidelines[:3]:
            formatted.append(
                f"- {guideline['name']} ({guideline['authority']}, {guideline['year']})"
            )
        
        return "\n".join(formatted)
    
    def _extract_clinical_context(self, signal_data: Dict) -> str:
        """Extract clinical context, filtering NaNs."""
        items = []
        
        # Check clinical fields
        clinical_fields = {
            'causality_assessment': 'Causality',
            'seriousness': 'Seriousness',
            'outcome': 'Outcome',
            'dechallenge': 'Dechallenge',
            'rechallenge': 'Rechallenge',
            'alt_u_l': 'ALT (U/L)',
            'ast_u_l': 'AST (U/L)',
            'bilirubin_mg_dl': 'Bilirubin (mg/dL)',
            'creatinine_mg_dl': 'Creatinine (mg/dL)'
        }
        
        for field, label in clinical_fields.items():
            value = signal_data.get(field)
            if value and str(value).lower() not in ['nan', 'none', '', 'unknown']:
                if isinstance(value, (int, float)) and not pd.isna(value):
                    items.append(f"{label}: {value}")
                elif isinstance(value, str):
                    items.append(f"{label}: {value}")
        
        if not items:
            return "Limited clinical context available. Assessment based on disproportionality analysis only."
        
        return "; ".join(items)
    
    def _extract_evidence_sources(self, drug_context: Dict, similar_signals: List) -> List[str]:
        """Extract evidence sources for citation."""
        sources = [
            "Knowledge Graph: Regulatory guidelines (FDA, EMA, CIOMS XIV)",
            f"Similar signals: {len(similar_signals)} known associations"
        ]
        
        if drug_context.get('known_events'):
            sources.append(f"Known events for drug: {len(drug_context['known_events'])} documented")
        
        if drug_context.get('signals'):
            sources.append(f"Confirmed signals: {len(drug_context['signals'])} from training data")
        
        return sources
