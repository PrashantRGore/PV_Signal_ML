"""
Lightweight Knowledge Graph for Causal Feature Engineering
Uses NetworkX + JSON/SQLite (no Neo4j required)
Integrates with existing RegulatoryKnowledgeGraph from app_production_v3.py

Purpose: Generate 10 causal plausibility features for signal detection
Author: PV Signal ML Team
Date: January 10, 2026
Version: 3.0
"""

import json
import sqlite3
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
import networkx as nx
import pandas as pd


class LightweightKnowledgeGraph:
    """
    Lightweight knowledge graph for drug-AE causal reasoning
    
    Features:
    - Drug → Metabolite → Organ → AE pathways
    - Drug-Drug Interaction (DDI) database
    - Organ-specific toxicity patterns
    - SmPC/Literature mechanism scores
    
    Storage: NetworkX MultiDiGraph + JSON files
    """
    
    def __init__(self, data_dir: str = None):
        """
        Initialize knowledge graph
        
        Args:
            data_dir: Directory containing KG JSON files
                     Default: src/causal/kg_data/
        """
        if data_dir is None:
            data_dir = Path(__file__).parent / "kg_data"
        
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize empty graph
        self.graph = nx.MultiDiGraph()
        
        # Load knowledge base
        self._load_drug_mechanisms()
        self._load_ddi_database()
        self._load_organ_tropism()
    
    # ========================================================================
    # DATA LOADING METHODS
    # ========================================================================
    
    def _load_drug_mechanisms(self):
        """
        Load drug-metabolite-organ-AE pathways from JSON
        
        JSON Structure:
        {
            "Paracetamol": {
                "metabolites": [
                    {
                        "name": "NAPQI",
                        "targets": ["Liver"],
                        "mechanism": "Glutathione depletion",
                        "evidence": "well-established"
                    }
                ],
                "known_aes": [
                    {
                        "name": "Hepatotoxicity",
                        "meddra_code": "10019851",
                        "plausibility": 0.95,
                        "source": "SmPC"
                    }
                ]
            }
        }
        """
        mechanisms_file = self.data_dir / "drug_mechanisms.json"
        
        if not mechanisms_file.exists():
            # Create default file with examples
            self._create_default_mechanisms()
            return
        
        try:
            with open(mechanisms_file, 'r', encoding='utf-8') as f:
                data = json.load(f)
        except json.JSONDecodeError:
            print(f"Warning: Invalid JSON in {mechanisms_file}. Creating default file.")
            self._create_default_mechanisms()
            return
        
        # Build graph edges from JSON data
        for drug, info in data.items():
            drug_node = f"Drug:{drug}"
            self.graph.add_node(drug_node, type='drug', name=drug)
            
            # Add metabolites
            for metabolite in info.get('metabolites', []):
                met_node = f"Metabolite:{metabolite['name']}"
                self.graph.add_node(met_node, type='metabolite', name=metabolite['name'])
                self.graph.add_edge(
                    drug_node, met_node,
                    relation='metabolizes_to',
                    evidence=metabolite.get('evidence', 'literature')
                )
                
                # Add organ targets
                for organ in metabolite.get('targets', []):
                    organ_node = f"Organ:{organ}"
                    self.graph.add_node(organ_node, type='organ', name=organ)
                    self.graph.add_edge(
                        met_node, organ_node,
                        relation='affects',
                        mechanism=metabolite.get('mechanism', 'unknown')
                    )
            
            # Add direct drug-AE associations (from SmPC/literature)
            for ae in info.get('known_aes', []):
                ae_node = f"AE:{ae['meddra_code']}"
                self.graph.add_node(
                    ae_node, 
                    type='adverse_event', 
                    name=ae['name'], 
                    meddra=ae['meddra_code']
                )
                self.graph.add_edge(
                    drug_node, ae_node,
                    relation='causes',
                    plausibility=ae.get('plausibility', 0.5),
                    source=ae.get('source', 'SmPC')
                )
    
    def _load_ddi_database(self):
        """
        Load drug-drug interaction (DDI) database from JSON
        
        JSON Structure:
        [
            {
                "drug1": "Warfarin",
                "drug2": "Aspirin",
                "severity": "major",
                "mechanism": "Additive bleeding risk",
                "risk_score": 0.8
            }
        ]
        """
        ddi_file = self.data_dir / "ddi_database.json"
        
        if not ddi_file.exists():
            self._create_default_ddi()
            return
        
        try:
            with open(ddi_file, 'r', encoding='utf-8') as f:
                data = json.load(f)
        except json.JSONDecodeError:
            print(f"Warning: Invalid JSON in {ddi_file}. Creating default file.")
            self._create_default_ddi()
            return

        
        for interaction in data:
            drug1 = f"Drug:{interaction['drug1']}"
            drug2 = f"Drug:{interaction['drug2']}"
            
            # Add nodes if they don't exist
            if not self.graph.has_node(drug1):
                self.graph.add_node(drug1, type='drug', name=interaction['drug1'])
            if not self.graph.has_node(drug2):
                self.graph.add_node(drug2, type='drug', name=interaction['drug2'])
            
            # Add bidirectional interaction edge
            self.graph.add_edge(
                drug1, drug2,
                relation='interacts_with',
                severity=interaction.get('severity', 'moderate'),
                mechanism=interaction.get('mechanism', ''),
                risk_score=interaction.get('risk_score', 0.5)
            )
    
    def _load_organ_tropism(self):
        """
        Load drug organ-specific toxicity patterns from JSON
        
        JSON Structure:
        {
            "Paracetamol": [
                {"organ": "Liver", "specificity": 0.9}
            ],
            "Gentamicin": [
                {"organ": "Kidney", "specificity": 0.85}
            ]
        }
        """
        tropism_file = self.data_dir / "organ_tropism.json"
        
        if not tropism_file.exists():
            self._create_default_tropism()
            return
        
        try:
            with open(tropism_file, 'r', encoding='utf-8') as f:
                data = json.load(f)
        except json.JSONDecodeError:
            print(f"Warning: Invalid JSON in {tropism_file}. Creating default file.")
            self._create_default_tropism()
            return
        
        for drug, organs in data.items():
            drug_node = f"Drug:{drug}"
            if not self.graph.has_node(drug_node):
                self.graph.add_node(drug_node, type='drug', name=drug)
            
            for organ_info in organs:
                organ_node = f"Organ:{organ_info['organ']}"
                if not self.graph.has_node(organ_node):
                    self.graph.add_node(organ_node, type='organ', name=organ_info['organ'])
                
                self.graph.add_edge(
                    drug_node, organ_node,
                    relation='targets',
                    specificity=organ_info.get('specificity', 0.5)
                )
    
    # ========================================================================
    # QUERY METHODS (Used by Causal Feature Engineer)
    # ========================================================================
    
    def get_shortest_path(self, drug: str, ae: str, max_depth: int = 4) -> Optional[List]:
        """
        Find shortest path from drug to adverse event
        
        Args:
            drug: Drug name (e.g., "Paracetamol")
            ae: Adverse event name or MedDRA code (e.g., "Hepatotoxicity" or "10019851")
            max_depth: Maximum path length to search (default: 4 edges)
        
        Returns:
            List of nodes in path, or None if no path exists
            
        Example:
            ["Drug:Paracetamol", "Metabolite:NAPQI", "Organ:Liver", "AE:10019851"]
        """
        drug_node = f"Drug:{drug}"
        
        # Find matching AE nodes (by name or code)
        ae_nodes = [
            n for n in self.graph.nodes() 
            if n.startswith('AE:') and (
                ae.lower() in self.graph.nodes[n].get('name', '').lower() or
                ae in n
            )
        ]
        
        if not ae_nodes:
            return None
        
        if not self.graph.has_node(drug_node):
            return None
        
        # Try finding path to any matching AE node
        try:
            for ae_candidate in ae_nodes:
                try:
                    path = nx.shortest_path(self.graph, drug_node, ae_candidate)
                    if len(path) <= max_depth + 1:  # +1 for node count vs edge count
                        return path
                except nx.NetworkXNoPath:
                    continue
            return None
        except Exception:
            return None
    
    def pathway_exists(self, drug: str, ae: str) -> bool:
        """
        Check if any pathway exists from drug to AE
        
        Returns:
            True if mechanistic pathway exists, False otherwise
        """
        path = self.get_shortest_path(drug, ae)
        return path is not None
    
    def get_pathway_length(self, drug: str, ae: str) -> int:
        """
        Get length of shortest pathway (number of edges)
        
        Returns:
            Number of edges in shortest path (1-4)
            Returns 99 if no path exists
        """
        path = self.get_shortest_path(drug, ae)
        if path is None:
            return 99
        return len(path) - 1  # edges = nodes - 1
    
    def get_mechanism_score(self, drug: str, ae: str) -> float:
        """
        Get mechanistic plausibility score (0-1)
        
        Scoring:
        - Direct edge with plausibility annotation: Use that score
        - Path length 1-2: High plausibility (0.7-0.9)
        - Path length 3-4: Moderate plausibility (0.3-0.6)
        - No path: Low plausibility (0.1)
        
        Returns:
            Float between 0.0-1.0 indicating mechanism strength
        """
        path = self.get_shortest_path(drug, ae)
        if path is None:
            return 0.1  # Low plausibility
        
        # Score based on path length (shorter = more plausible)
        length = len(path) - 1
        base_score = 1.0 / (1.0 + length)
        
        # Boost score if direct edge has plausibility annotation
        drug_node = f"Drug:{drug}"
        ae_nodes = [n for n in path if n.startswith('AE:')]
        
        if ae_nodes and self.graph.has_edge(drug_node, ae_nodes[0]):
            edge_data = self.graph.get_edge_data(drug_node, ae_nodes[0])
            if edge_data:
                for key, data in edge_data.items():
                    if 'plausibility' in data:
                        return data['plausibility']
        
        return min(base_score, 0.9)
    
    def get_confounder_drugs(self, drug: str, ae: str, top_k: int = 5) -> List[str]:
        """
        Get list of concomitant drugs that could confound the drug-AE association
        (i.e., drugs that also cause the same AE)
        
        Args:
            drug: Primary drug being evaluated
            ae: Adverse event of interest
            top_k: Maximum number of confounders to return
        
        Returns:
            List of drug names that are known to cause the same AE
        """
        ae_nodes = [
            n for n in self.graph.nodes() 
            if n.startswith('AE:') and 
            ae.lower() in self.graph.nodes[n].get('name', '').lower()
        ]
        
        if not ae_nodes:
            return []
        
        confounders = []
        for ae_node in ae_nodes:
            # Find all drugs that cause this AE
            predecessors = list(self.graph.predecessors(ae_node))
            for pred in predecessors:
                if pred.startswith('Drug:') and pred != f"Drug:{drug}":
                    drug_name = pred.replace('Drug:', '')
                    confounders.append(drug_name)
        
        return list(set(confounders))[:top_k]
    
    def get_ddi_risk(self, drug: str, concomitants: List[str]) -> float:
        """
        Calculate drug-drug interaction risk score
        
        Args:
            drug: Primary drug
            concomitants: List of concomitant drug names
        
        Returns:
            Risk score 0-1 based on known interactions
        """
        if not concomitants:
            return 0.0
        
        drug_node = f"Drug:{drug}"
        if not self.graph.has_node(drug_node):
            return 0.0
        
        risk_scores = []
        for concom in concomitants[:5]:  # Top 5 concomitants only
            concom_node = f"Drug:{concom}"
            if self.graph.has_edge(drug_node, concom_node):
                edge_data = self.graph.get_edge_data(drug_node, concom_node)
                for key, data in edge_data.items():
                    if data.get('relation') == 'interacts_with':
                        risk_scores.append(data.get('risk_score', 0.5))
        
        if not risk_scores:
            return 0.0
        
        # Average risk, capped at 1.0
        return min(sum(risk_scores) / len(concomitants), 1.0)

    
    def check_organ_specificity(self, drug: str, ae_organ: str) -> bool:
        """
        Check if drug has known organ-specific toxicity matching the AE organ
        
        Args:
            drug: Drug name
            ae_organ: Organ affected by AE (e.g., 'Liver', 'Kidney', 'Heart')
        
        Returns:
            True if drug is known to target this organ
        """
        drug_node = f"Drug:{drug}"
        if not self.graph.has_node(drug_node):
            return False
        
        # Check direct organ targeting
        successors = list(self.graph.successors(drug_node))
        for succ in successors:
            if succ.startswith('Organ:'):
                organ_name = succ.replace('Organ:', '')
                if ae_organ.lower() in organ_name.lower():
                    return True
        
        return False
    
    def get_graph_stats(self) -> Dict[str, int]:
        """
        Get knowledge graph statistics for debugging
        
        Returns:
            Dict with node/edge counts by type
        """
        node_types = {}
        for node, data in self.graph.nodes(data=True):
            node_type = data.get('type', 'unknown')
            node_types[node_type] = node_types.get(node_type, 0) + 1
        
        return {
            'total_nodes': self.graph.number_of_nodes(),
            'total_edges': self.graph.number_of_edges(),
            'node_types': node_types
        }
    
    # ========================================================================
    # DEFAULT DATA CREATION METHODS
    # ========================================================================
    
    def _create_default_mechanisms(self):
        """
        Create default drug mechanisms file with common examples
        This is used as starter data - can be expanded with real data
        """
        default_data = {
            "Paracetamol": {
                "metabolites": [
                    {
                        "name": "NAPQI",
                        "targets": ["Liver"],
                        "mechanism": "Glutathione depletion",
                        "evidence": "well-established"
                    }
                ],
                "known_aes": [
                    {
                        "name": "Hepatotoxicity",
                        "meddra_code": "10019851",
                        "plausibility": 0.95,
                        "source": "SmPC"
                    },
                    {
                        "name": "Acute Liver Failure",
                        "meddra_code": "10000804",
                        "plausibility": 0.90,
                        "source": "FDA_Label"
                    }
                ]
            },
            "Metformin": {
                "metabolites": [],
                "known_aes": [
                    {
                        "name": "Lactic Acidosis",
                        "meddra_code": "10023930",
                        "plausibility": 0.85,
                        "source": "SmPC"
                    }
                ]
            },
            "Warfarin": {
                "metabolites": [],
                "known_aes": [
                    {
                        "name": "Hemorrhage",
                        "meddra_code": "10055798",
                        "plausibility": 0.92,
                        "source": "SmPC"
                    },
                    {
                        "name": "Gastrointestinal Hemorrhage",
                        "meddra_code": "10017955",
                        "plausibility": 0.88,
                        "source": "SmPC"
                    }
                ]
            },
            "Gentamicin": {
                "metabolites": [],
                "known_aes": [
                    {
                        "name": "Nephrotoxicity",
                        "meddra_code": "10029095",
                        "plausibility": 0.90,
                        "source": "SmPC"
                    },
                    {
                        "name": "Ototoxicity",
                        "meddra_code": "10033441",
                        "plausibility": 0.85,
                        "source": "SmPC"
                    }
                ]
            },
            "Doxorubicin": {
                "metabolites": [],
                "known_aes": [
                    {
                        "name": "Cardiomyopathy",
                        "meddra_code": "10007554",
                        "plausibility": 0.87,
                        "source": "SmPC"
                    }
                ]
            },
            "Isoniazid": {
                "metabolites": [
                    {
                        "name": "Acetylhydrazine",
                        "targets": ["Liver"],
                        "mechanism": "Reactive metabolite formation",
                        "evidence": "established"
                    }
                ],
                "known_aes": [
                    {
                        "name": "Hepatotoxicity",
                        "meddra_code": "10019851",
                        "plausibility": 0.82,
                        "source": "SmPC"
                    }
                ]
            }
        }
        
        with open(self.data_dir / "drug_mechanisms.json", 'w', encoding='utf-8') as f:
            json.dump(default_data, f, indent=2)
        
        # Reload after creation
        self._load_drug_mechanisms()
    
    def _create_default_ddi(self):
        """
        Create default DDI database with common interactions
        """
        default_data = [
            {
                "drug1": "Warfarin",
                "drug2": "Aspirin",
                "severity": "major",
                "mechanism": "Additive bleeding risk",
                "risk_score": 0.8
            },
            {
                "drug1": "Metformin",
                "drug2": "Contrast Media",
                "severity": "moderate",
                "mechanism": "Increased lactic acidosis risk",
                "risk_score": 0.6
            },
            {
                "drug1": "Warfarin",
                "drug2": "NSAIDs",
                "severity": "major",
                "mechanism": "Increased bleeding risk",
                "risk_score": 0.75
            },
            {
                "drug1": "Digoxin",
                "drug2": "Amiodarone",
                "severity": "major",
                "mechanism": "Digoxin toxicity",
                "risk_score": 0.7
            },
            {
                "drug1": "ACE Inhibitors",
                "drug2": "Potassium Supplements",
                "severity": "moderate",
                "mechanism": "Hyperkalemia risk",
                "risk_score": 0.65
            }
        ]
        
        with open(self.data_dir / "ddi_database.json", 'w', encoding='utf-8') as f:
            json.dump(default_data, f, indent=2)
    
    def _create_default_tropism(self):
        """
        Create default organ tropism file
        """
        default_data = {
            "Paracetamol": [
                {"organ": "Liver", "specificity": 0.9}
            ],
            "Gentamicin": [
                {"organ": "Kidney", "specificity": 0.85},
                {"organ": "Ear", "specificity": 0.80}
            ],
            "Doxorubicin": [
                {"organ": "Heart", "specificity": 0.8}
            ],
            "Methotrexate": [
                {"organ": "Liver", "specificity": 0.75},
                {"organ": "Kidney", "specificity": 0.70}
            ],
            "Isoniazid": [
                {"organ": "Liver", "specificity": 0.82}
            ],
            "Cisplatin": [
                {"organ": "Kidney", "specificity": 0.88}
            ],
            "Amiodarone": [
                {"organ": "Lung", "specificity": 0.78},
                {"organ": "Thyroid", "specificity": 0.75}
            ]
        }
        
        with open(self.data_dir / "organ_tropism.json", 'w', encoding='utf-8') as f:
            json.dump(default_data, f, indent=2)

    
    # ========================================================================
    # UTILITY METHODS
    # ========================================================================
    
    def add_drug_mechanism(self, drug: str, mechanism_data: Dict):
        """
        Add a new drug mechanism to the knowledge graph
        Useful for expanding the KG with project-specific data
        
        Args:
            drug: Drug name
            mechanism_data: Dict with metabolites and known_aes
        
        Example:
            kg.add_drug_mechanism("NewDrug", {
                "metabolites": [],
                "known_aes": [
                    {
                        "name": "Some AE",
                        "meddra_code": "12345678",
                        "plausibility": 0.7,
                        "source": "Literature"
                    }
                ]
            })
        """
        mechanisms_file = self.data_dir / "drug_mechanisms.json"
        
        # Load existing data
        if mechanisms_file.exists():
            with open(mechanisms_file, 'r', encoding='utf-8') as f:
                data = json.load(f)
        else:
            data = {}
        
        # Add new mechanism
        data[drug] = mechanism_data
        
        # Save
        with open(mechanisms_file, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2)
        
        # Reload graph
        self.graph.clear()
        self._load_drug_mechanisms()
        self._load_ddi_database()
        self._load_organ_tropism()
    
    def export_graph(self, output_path: str = None):
        """
        Export knowledge graph to GraphML format for visualization
        Can be opened in Gephi, Cytoscape, or other graph tools
        
        Args:
            output_path: Path to save GraphML file
                        Default: src/causal/kg_data/knowledge_graph.graphml
        """
        if output_path is None:
            output_path = self.data_dir / "knowledge_graph.graphml"
        
        nx.write_graphml(self.graph, output_path)
        return output_path
    
    def query_drug_info(self, drug: str) -> Dict[str, Any]:
        """
        Get all information about a specific drug from the knowledge graph
        
        Args:
            drug: Drug name
        
        Returns:
            Dict with pathways, interactions, and organ targets
        """
        drug_node = f"Drug:{drug}"
        
        if not self.graph.has_node(drug_node):
            return {"exists": False, "message": f"Drug '{drug}' not found in knowledge graph"}
        
        # Get all edges from this drug
        metabolites = []
        known_aes = []
        interactions = []
        organs = []
        
        for successor in self.graph.successors(drug_node):
            edge_data = self.graph.get_edge_data(drug_node, successor)
            
            if successor.startswith('Metabolite:'):
                metabolites.append(successor.replace('Metabolite:', ''))
            elif successor.startswith('AE:'):
                ae_info = self.graph.nodes[successor]
                known_aes.append({
                    'name': ae_info.get('name'),
                    'meddra': ae_info.get('meddra'),
                    'plausibility': list(edge_data.values())[0].get('plausibility', 0.5)
                })
            elif successor.startswith('Organ:'):
                organs.append(successor.replace('Organ:', ''))
            elif successor.startswith('Drug:'):
                interactions.append(successor.replace('Drug:', ''))
        
        return {
            "exists": True,
            "drug": drug,
            "metabolites": metabolites,
            "known_aes": known_aes,
            "organ_targets": organs,
            "interactions": interactions,
            "total_connections": self.graph.degree(drug_node)
        }


# ============================================================================
# CONVENIENCE FUNCTIONS
# ============================================================================

def load_knowledge_graph(data_dir: str = None) -> LightweightKnowledgeGraph:
    """
    Load knowledge graph (convenience function)
    
    Args:
        data_dir: Optional custom data directory
    
    Returns:
        LightweightKnowledgeGraph instance
    """
    return LightweightKnowledgeGraph(data_dir=data_dir)


if __name__ == "__main__":
    # Test script - run with: python -m src.causal.knowledge_graph
    print("Testing Lightweight Knowledge Graph...")
    
    kg = LightweightKnowledgeGraph()
    
    print("\n📊 Knowledge Graph Statistics:")
    stats = kg.get_graph_stats()
    for key, value in stats.items():
        print(f"  {key}: {value}")
    
    print("\n🔍 Test Query: Paracetamol → Hepatotoxicity")
    exists = kg.pathway_exists("Paracetamol", "Hepatotoxicity")
    print(f"  Pathway exists: {exists}")
    
    if exists:
        length = kg.get_pathway_length("Paracetamol", "Hepatotoxicity")
        score = kg.get_mechanism_score("Paracetamol", "Hepatotoxicity")
        print(f"  Pathway length: {length} edges")
        print(f"  Mechanism score: {score:.2f}")
    
    print("\n🔍 Test Query: Warfarin confounders for Hemorrhage")
    confounders = kg.get_confounder_drugs("Warfarin", "Hemorrhage")
    print(f"  Potential confounders: {confounders}")
    
    print("\n✅ Knowledge Graph test complete!")

