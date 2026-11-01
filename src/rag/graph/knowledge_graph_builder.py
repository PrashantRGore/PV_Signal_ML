"""
Regulatory Knowledge Graph Builder
Creates Neo4j-compatible graph for Graph RAG retrieval
Nodes: Drugs, Events, Guidelines, Signals
Edges: causes, contraindicated, cited_by, references
"""

import pandas as pd
from typing import Dict, List, Tuple, Optional
from pathlib import Path
import json
from datetime import datetime
from src.utils.logger import setup_logger

logger = setup_logger(__name__)


class RegulatoryKnowledgeGraph:
    """
    Build knowledge graph from:
    - Known drug-event associations
    - FDA guidance documents
    - EMA GVP modules
    - CIOMS guidelines
    - MedDRA hierarchy
    """
    
    def __init__(self, graph_dir: Path = Path("data/knowledge_graph")):
        self.graph_dir = graph_dir
        self.graph_dir.mkdir(parents=True, exist_ok=True)
        
        self.nodes = {
            'drugs': [],
            'events': [],
            'guidelines': [],
            'signals': []
        }
        
        self.edges = []
        
    def add_known_signals(self, signals_df: pd.DataFrame):
        """
        Add known drug-event signals from training data.

        Args:
            signals_df: DataFrame with columns [drug_name, event_term, prr, label]
        """
        # ✅ Validate required columns
        required = ['drug_name', 'event_term']
        missing = [c for c in required if c not in signals_df.columns]
        
        if missing:
            logger.error(f"Missing required columns: {missing}")
            return
        
        logger.info(f"Adding {len(signals_df)} signals to knowledge graph")
        
        for _, row in signals_df.iterrows():
            # ✅ Safe column access
            drug = str(row.get('drug_name', 'UNKNOWN')).upper()
            event = str(row.get('event_term', 'UNKNOWN')).upper()
            
            if drug == 'UNKNOWN' or event == 'UNKNOWN':
                continue
            
            prr = float(row.get('prr', 0))
            label = int(row.get('label', 0))
            
            # Add drug node
            if drug not in [n['name'] for n in self.nodes['drugs']]:
                self.nodes['drugs'].append({
                    'name': drug,
                    'type': 'DRUG',
                    'id': f"DRUG_{len(self.nodes['drugs'])}"
                })
            
            # Add event node
            if event not in [n['name'] for n in self.nodes['events']]:
                self.nodes['events'].append({
                    'name': event,
                    'type': 'EVENT',
                    'id': f"EVENT_{len(self.nodes['events'])}"
                })
            
            # Add signal node if confirmed
            if label == 1:
                self.nodes['signals'].append({
                    'drug': drug,
                    'event': event,
                    'prr': prr,
                    'type': 'SIGNAL',
                    'id': f"SIGNAL_{len(self.nodes['signals'])}"
                })
                
                # Add edge: DRUG -[CAUSES]-> EVENT
                self.edges.append({
                    'from': drug,
                    'to': event,
                    'relationship': 'CAUSES',
                    'strength': prr,
                    'evidence': 'training_data'
                })
        
        logger.info(f"✅ Added {len(self.nodes['signals'])} confirmed signals")
    
    def add_regulatory_guidelines(self):
        """
        Add regulatory guidance documents as nodes.
        
        Guidelines:
        - FDA: Drug Safety Communication, MedWatch
        - EMA: GVP Module IX (Signal Management)
        - CIOMS: CIOMS XIV (AI in PV)
        - WHO-UMC: Causality Assessment
        """
        guidelines = [
            {
                'id': 'FDA_AI_PV_2024',
                'name': 'FDA Guidance on AI in Pharmacovigilance',
                'authority': 'FDA',
                'year': 2024,
                'url': 'https://www.fda.gov/regulatory-information/search-fda-guidance-documents',
                'key_points': [
                    'AI models must be explainable (SHAP, LIME)',
                    'Training data must be documented',
                    'Model versioning required',
                    'Audit trail recommended (Part 11-inspired)'
                ]
            },
            {
                'id': 'EMA_GVP_MODULE_IX',
                'name': 'EMA GVP Module IX - Signal Management',
                'authority': 'EMA',
                'year': 2017,
                'url': 'https://www.ema.europa.eu/en/human-regulatory/post-authorisation/pharmacovigilance',
                'key_points': [
                    'Signal detection thresholds: PRR >= 2, Chi² >= 4',
                    'Causality assessment required',
                    'Medical review mandatory',
                    'Quarterly signal reports'
                ]
            },
            {
                'id': 'CIOMS_XIV_AI',
                'name': 'CIOMS XIV - AI and Machine Learning in Pharmacovigilance',
                'authority': 'CIOMS',
                'year': 2024,
                'url': 'https://cioms.ch/',
                'key_points': [
                    'Missing data handling must follow ICH E9',
                    'Bias mitigation required',
                    'Algorithm validation against gold standard',
                    'Human oversight essential'
                ]
            },
            {
                'id': 'WHO_UMC_CAUSALITY',
                'name': 'WHO-UMC Causality Assessment System',
                'authority': 'WHO',
                'year': 2018,
                'url': 'https://www.who-umc.org/',
                'key_points': [
                    'Certain: Event follows drug, plausible, dechallenge/rechallenge positive',
                    'Probable/Likely: Event follows drug, plausible',
                    'Possible: Event follows drug, other causes possible',
                    'Unlikely: Temporal relationship poor'
                ]
            },
            {
                'id': 'ICH_E2B_R3',
                'name': 'ICH E2B(R3) - ICSR Transmission',
                'authority': 'ICH',
                'year': 2016,
                'url': 'https://www.ich.org/',
                'key_points': [
                    'Standardized ICSR format',
                    'MedDRA coding required',
                    'Minimum data elements defined',
                    'Electronic transmission standards'
                ]
            }
        ]
        
        for guideline in guidelines:
            self.nodes['guidelines'].append(guideline)
            
            # Create edges to relevant drugs/events mentioned in guidelines
            # (This would be enhanced with NLP extraction in production)
        
        logger.info(f"✅ Added {len(guidelines)} regulatory guidelines")
    
    def add_meddra_hierarchy(self):
        """
        Add simplified MedDRA hierarchy for event relationships.
        
        Example: Myocardial Infarction (PT) -> Ischaemic heart disease (HLT) 
                 -> Cardiac disorders (SOC)
        """
        # Simplified examples (production would use full MedDRA)
        meddra_relationships = [
            ('Myocardial infarction', 'Ischaemic heart disease', 'PT_TO_HLT'),
            ('Cardiac failure', 'Heart failures', 'PT_TO_HLT'),
            ('Hepatic failure', 'Liver failures', 'PT_TO_HLT'),
            ('Acute kidney injury', 'Renal failures', 'PT_TO_HLT'),
            ('Rhabdomyolysis', 'Muscle disorders', 'PT_TO_HLT'),
        ]
        
        for pt, hlt, rel_type in meddra_relationships:
            self.edges.append({
                'from': pt,
                'to': hlt,
                'relationship': rel_type,
                'source': 'MedDRA'
            })
        
        logger.info(f"✅ Added {len(meddra_relationships)} MedDRA relationships")
    
    def export_to_json(self):
        """Export knowledge graph to JSON for Neo4j import or local graph DB."""
        graph_data = {
            'metadata': {
                'created': datetime.now().isoformat(),
                'version': '1.0',
                'node_count': sum(len(nodes) for nodes in self.nodes.values()),
                'edge_count': len(self.edges)
            },
            'nodes': self.nodes,
            'edges': self.edges
        }
        
        output_path = self.graph_dir / "regulatory_kg.json"
        with open(output_path, 'w') as f:
            json.dump(graph_data, f, indent=2)
        
        logger.info(f"✅ Knowledge graph exported to {output_path}")
        return output_path
    
    def get_drug_context(self, drug_name: str) -> Dict:
        """
        Retrieve all context for a drug from knowledge graph.
        
        Args:
            drug_name: Drug name to query
            
        Returns:
            Dict with known events, guidelines, signals
        """
        context = {
            'drug': drug_name,
            'known_events': [],
            'guidelines': [],
            'signals': []
        }
        
        # Find edges where drug is source
        for edge in self.edges:
            if edge['from'].upper() == drug_name.upper():
                context['known_events'].append({
                    'event': edge['to'],
                    'relationship': edge['relationship'],
                    'strength': edge.get('strength', 0)
                })
        
        # Find signals
        for signal in self.nodes['signals']:
            if signal['drug'].upper() == drug_name.upper():
                context['signals'].append(signal)
        
        return context
    
    def query_similar_signals(self, drug: str, event: str, top_k: int = 5) -> List[Dict]:
        """
        Find similar known signals using graph traversal.
        
        Args:
            drug: Query drug
            event: Query event
            top_k: Number of similar signals to return
            
        Returns:
            List of similar signals with context
        """
        similar = []
        
        # Find signals with same event
        for signal in self.nodes['signals']:
            if signal['event'].upper() == event.upper():
                similar.append({
                    'drug': signal['drug'],
                    'event': signal['event'],
                    'prr': signal['prr'],
                    'similarity_reason': 'same_event'
                })
        
        # Find signals with same drug
        for signal in self.nodes['signals']:
            if signal['drug'].upper() == drug.upper():
                if signal not in similar:
                    similar.append({
                        'drug': signal['drug'],
                        'event': signal['event'],
                        'prr': signal['prr'],
                        'similarity_reason': 'same_drug'
                    })
        
        # Sort by PRR and return top_k
        similar = sorted(similar, key=lambda x: x['prr'], reverse=True)[:top_k]
        
        return similar
