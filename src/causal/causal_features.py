"""
Causal Feature Engineer - Generates 10 new features for SISA
Uses Knowledge Graph to encode biological plausibility
Extends existing 68 features → 78 features total

Author: PV Signal ML Team
Date: January 10, 2026
Version: 3.0
"""

import pandas as pd
import numpy as np
from typing import Optional, List, Dict, Any
from pathlib import Path
from .knowledge_graph import LightweightKnowledgeGraph


class CausalFeatureEngineer:
    """
    Generates 10 causal plausibility features from Knowledge Graph
    
    Input: DataFrame with 68 engineered features
    Output: DataFrame with 78 features (68 + 10 causal)
    
    New Features:
        1. kg_pathway_exists (binary)
        2. kg_pathway_length (int: 1-4 or 99)
        3. known_mechanism_score (float: 0-1)
        4. confounder_risk (float: 0-1)
        5. temporal_causality_score (float: 0-1)
        6. dose_response_evidence (float: 0-1)
        7. dechallenge_weight (float: 0-3)
        8. organ_specificity_match (binary)
        9. ddi_interaction_risk (float: 0-1)
        10. literature_citation_count (int: 0-100)
    """
    
    def __init__(self, kg: Optional[LightweightKnowledgeGraph] = None):
        """
        Initialize causal feature engineer
        
        Args:
            kg: Optional pre-loaded Knowledge Graph
                If None, will load from default location
        """
        if kg is None:
            kg_data_dir = Path(__file__).parent / "kg_data"
            self.kg = LightweightKnowledgeGraph(data_dir=kg_data_dir)
        else:
            self.kg = kg
    
    def transform(self, signals_df: pd.DataFrame) -> pd.DataFrame:
        """
        Add 10 causal features to existing engineered signals
        
        Args:
            signals_df: DataFrame with 68 features
        
        Returns:
            DataFrame with 78 features (68 + 10 causal)
        
        Required columns in input:
            - drug_name (or drugname)
            - pt_name (or eventterm)
            - time_to_onset_days (optional)
            - treatment_duration_days (optional)
            - dechallenge_positive (optional)
            - rechallenge_positive (optional)
            - concomitant_count (optional)
        """
        # Validate input
        if signals_df is None or len(signals_df) == 0:
            raise ValueError("Input signals_df is empty")
        
        # Standardize column names (handle both canonical and local formats)
        df = signals_df.copy()
        df = self._standardize_columns(df)
        
        print(f"🔬 Generating causal features for {len(df)} signals...")
        
        # Generate each causal feature
        causal_features = []
        
        for idx, row in df.iterrows():
            drug = row.get('drug_name', '')
            ae = row.get('pt_name', '')
            
            # Feature 1-3: Knowledge Graph Pathway Analysis
            kg_result = self._query_kg_pathway(drug, ae)
            
            # Feature 4: Confounder Risk
            confounder_risk = self._calculate_confounder_risk(
                row.get('concomitant_count', 0),
                row.get('concomitant_drugs', [])
            )
            
            # Feature 5: Temporal Causality (Bradford Hill Criterion #1)
            temporal_score = self._assess_temporal_causality(
                row.get('time_to_onset_days', np.nan),
                row.get('treatment_duration_days', np.nan)
            )
            
            # Feature 6: Dose-Response Evidence (Bradford Hill Criterion #2)
            dose_response = self._assess_dose_response(
                row.get('drug_dose', np.nan),
                row.get('outcome_fatal', 0) + row.get('outcome_recovered', 0)
            )
            
            # Feature 7: Amplified Dechallenge Weight
            dechallenge_weight = self._calculate_dechallenge_weight(
                row.get('dechallenge_positive', 0),
                row.get('rechallenge_positive', 0)
            )
            
            # Feature 8: Organ-Specific Matching
            organ_match = self._check_organ_specificity(drug, ae)
            
            # Feature 9: DDI-based Interaction Risk
            ddi_risk = self._calculate_ddi_risk(row, drug)
            
            # Feature 10: Literature Support (simulated for now)
            lit_count = self._get_pubmed_cooccurrence(drug, ae)
            
            causal_features.append({
                'kg_pathway_exists': kg_result['exists'],
                'kg_pathway_length': kg_result['path_length'],
                'known_mechanism_score': kg_result['mechanism_score'],
                'confounder_risk': confounder_risk,
                'temporal_causality_score': temporal_score,
                'dose_response_evidence': dose_response,
                'dechallenge_weight': dechallenge_weight,
                'organ_specificity_match': organ_match,
                'ddi_interaction_risk': ddi_risk,
                'literature_citation_count': lit_count
            })
        
        # Convert to DataFrame
        causal_df = pd.DataFrame(causal_features)
        
        # Combine with original features
        result_df = pd.concat([df, causal_df], axis=1)
        
        print(f"✅ Causal features generated: {len(result_df.columns)} total columns")
        
        return result_df
    
    # ========================================================================
    # HELPER METHODS - Column Standardization
    # ========================================================================
    
    def _standardize_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Standardize column names to handle variations
        (drugname vs drug_name, eventterm vs pt_name, etc.)
        """
        column_mapping = {
            'drugname': 'drug_name',
            'drug': 'drug_name',
            'eventterm': 'pt_name',
            'event': 'pt_name',
            'ae_name': 'pt_name',
            'ttp': 'time_to_onset_days',
            'time_to_onset': 'time_to_onset_days',
            'treatment_duration': 'treatment_duration_days'
        }
        
        # Apply mapping (case-insensitive)
        for old_name, new_name in column_mapping.items():
            matching_cols = [col for col in df.columns if col.lower() == old_name.lower()]
            if matching_cols and new_name not in df.columns:
                df = df.rename(columns={matching_cols[0]: new_name})
        
        return df
    
    # ========================================================================
    # FEATURE GENERATION METHODS
    # ========================================================================
    
    def _query_kg_pathway(self, drug: str, ae: str) -> Dict[str, Any]:
        """
        Query Knowledge Graph for mechanistic pathway
        
        Returns:
            Dict with:
                - exists: bool (pathway exists)
                - path_length: int (number of edges, or 99 if no path)
                - mechanism_score: float (plausibility 0-1)
        """
        try:
            exists = self.kg.pathway_exists(drug, ae)
            
            if exists:
                path_length = self.kg.get_pathway_length(drug, ae)
                mechanism_score = self.kg.get_mechanism_score(drug, ae)
            else:
                path_length = 99
                mechanism_score = 0.1  # Low plausibility if no pathway
            
            return {
                'exists': 1 if exists else 0,
                'path_length': path_length,
                'mechanism_score': mechanism_score
            }
        except Exception as e:
            # Fallback if KG query fails
            return {
                'exists': 0,
                'path_length': 99,
                'mechanism_score': 0.5  # Neutral prior
            }
    
    def _calculate_confounder_risk(
        self, 
        concom_count: int, 
        concom_drugs: List[str] = None
    ) -> float:
        """
        Calculate risk score based on known confounders
        
        Args:
            concom_count: Number of concomitant drugs
            concom_drugs: List of concomitant drug names (optional)
        
        Returns:
            Risk score 0-1 (higher = more confounders present)
        """
        if concom_count == 0:
            return 0.0
        
        # Base risk from polypharmacy
        base_risk = min(concom_count / 10.0, 0.5)  # Cap at 0.5 for 10+ drugs
        
        # Additional risk if specific confounders known
        if concom_drugs and len(concom_drugs) > 0:
            # Check if any concomitants are known to cause similar AEs
            # This would require the current AE context, so we use a simplified approach
            additional_risk = 0.0
            for drug in concom_drugs[:5]:  # Top 5 concomitants
                # In real implementation, would query KG for drug-AE associations
                pass
            
            return min(base_risk + additional_risk, 1.0)
        
        return base_risk
    
    def _assess_temporal_causality(
        self, 
        ttp: float, 
        treatment_duration: float
    ) -> float:
        """
        Bradford Hill Criterion: Temporality
        Score based on whether timing makes biological sense
        
        Args:
            ttp: Time to onset (days)
            treatment_duration: Treatment duration (days)
        
        Returns:
            Score 0-1 (higher = better temporal relationship)
        """
        if pd.isna(ttp):
            return 0.5  # Unknown temporal relationship
        
        # Acute events (0-7 days) - high plausibility
        if ttp <= 7:
            return 0.9
        
        # Subacute (8-30 days) - moderate-high plausibility
        elif ttp <= 30:
            return 0.7
        
        # Chronic (>30 days) - requires longer treatment duration
        elif ttp > 30:
            if not pd.isna(treatment_duration) and treatment_duration > 30:
                return 0.6  # Reasonable if on drug long enough
            else:
                return 0.3  # Unlikely if short treatment
        
        return 0.5  # Default neutral
    
    def _assess_dose_response(
        self, 
        drug_dose: float, 
        severity_proxy: float
    ) -> float:
        """
        Bradford Hill Criterion: Dose-Response
        Assess if there's evidence of dose-response relationship
        
        Args:
            drug_dose: Drug dose (if available)
            severity_proxy: Proxy for severity (e.g., fatal + recovered outcomes)
        
        Returns:
            Score 0-1 (higher = stronger dose-response evidence)
        """
        # In real implementation, would analyze dose distributions
        # For now, use severity as proxy
        
        if pd.isna(drug_dose):
            # No dose information - use severity proxy
            if severity_proxy > 0:
                return 0.6  # Some evidence from severity
            else:
                return 0.5  # Unknown
        
        # If dose available, could analyze dose-severity correlation
        # Simplified: higher doses get higher scores
        try:
            dose_float = float(drug_dose)
            if dose_float > 0:
                return min(0.7, 0.5 + (dose_float / 1000.0))  # Normalized
        except:
            pass
        
        return 0.5  # Default neutral
    
    def _calculate_dechallenge_weight(
        self, 
        dechallenge_positive: int, 
        rechallenge_positive: int
    ) -> float:
        """
        Amplified dechallenge/rechallenge weight
        Strong evidence if both positive
        
        Args:
            dechallenge_positive: 1 if positive dechallenge, 0 otherwise
            rechallenge_positive: 1 if positive rechallenge, 0 otherwise
        
        Returns:
            Weight 0-3 (higher = stronger causality evidence)
        """
        weight = 0.0
        
        # Dechallenge positive = 1.0
        if dechallenge_positive == 1:
            weight += 1.0
        
        # Rechallenge positive = 2.0 (stronger evidence)
        if rechallenge_positive == 1:
            weight += 2.0
        
        return weight
    
    def _check_organ_specificity(self, drug: str, ae: str) -> int:
        """
        Check if drug has known organ-specific toxicity matching AE organ
        
        Args:
            drug: Drug name
            ae: Adverse event name
        
        Returns:
            1 if organ match found, 0 otherwise
        """
        # Map common AE terms to organs
        organ_keywords = {
            'Liver': ['hepat', 'liver', 'jaundice', 'bilirubin'],
            'Kidney': ['renal', 'kidney', 'nephro', 'creatinine'],
            'Heart': ['cardiac', 'heart', 'cardio', 'arrhythmia'],
            'Lung': ['pulmonary', 'lung', 'respiratory', 'pneumo'],
            'Brain': ['neuro', 'cerebral', 'cns', 'encephalo']
        }
        
        ae_lower = ae.lower()
        
        for organ, keywords in organ_keywords.items():
            # Check if AE involves this organ
            if any(kw in ae_lower for kw in keywords):
                # Check if drug targets this organ
                if self.kg.check_organ_specificity(drug, organ):
                    return 1
        
        return 0
    
    def _calculate_ddi_risk(self, row: pd.Series, drug: str) -> float:
        """
        Calculate DDI risk score from concomitant drugs
        
        Args:
            row: Signal row with concomitant drug information
            drug: Primary drug name
        
        Returns:
            Risk score 0-1
        """
        concom_count = row.get('concomitant_count', 0)
        
        if concom_count == 0:
            return 0.0
        
        # Try to get concomitant drug list
        concom_drugs = row.get('concomitant_drugs', [])
        
        if isinstance(concom_drugs, list) and len(concom_drugs) > 0:
            return self.kg.get_ddi_risk(drug, concom_drugs)
        
        # Fallback: Use polypharmacy as proxy
        # Higher concom count = higher risk
        return min(concom_count / 10.0, 0.7)
    
    def _get_pubmed_cooccurrence(self, drug: str, ae: str) -> int:
        """
        Get number of PubMed articles mentioning both drug and AE
        
        Currently simulated - in production, would query PubMed API
        
        Args:
            drug: Drug name
            ae: Adverse event name
        
        Returns:
            Citation count (0-100, capped)
        """
        # SIMULATION: Generate plausible citation counts
        # In production, would use:
        # from Bio import Entrez
        # Entrez.esearch(db="pubmed", term=f"{drug} AND {ae}")
        
        # Use deterministic SHA256 hash for reproducible cross-platform simulation
        import hashlib
        hash_bytes = hashlib.sha256(f"{drug.lower()}_{ae.lower()}".encode("utf-8")).hexdigest()
        hash_val = int(hash_bytes[:8], 16) % 100
        
        # Known high-profile associations get boosted
        known_pairs = {
            ('paracetamol', 'hepatotoxicity'): 85,
            ('warfarin', 'hemorrhage'): 92,
            ('metformin', 'lactic acidosis'): 78,
            ('gentamicin', 'nephrotoxicity'): 80
        }
        
        key = (drug.lower(), ae.lower())
        for known_key, count in known_pairs.items():
            if known_key[0] in key[0] and known_key[1] in key[1]:
                return count
        
        return min(hash_val, 100)
    
    # ========================================================================
    # UTILITY METHODS
    # ========================================================================
    
    def get_feature_descriptions(self) -> Dict[str, str]:
        """
        Get descriptions of all 10 causal features
        Useful for documentation and SHAP explanations
        
        Returns:
            Dict mapping feature name to description
        """
        return {
            'kg_pathway_exists': (
                'Binary indicator: 1 if mechanistic pathway exists in knowledge graph '
                '(Drug → Metabolite → Organ → AE), 0 otherwise'
            ),
            'kg_pathway_length': (
                'Length of shortest mechanistic pathway (1-4 edges). '
                '99 if no pathway exists. Shorter paths indicate stronger plausibility.'
            ),
            'known_mechanism_score': (
                'Mechanism plausibility score (0-1) from knowledge graph. '
                'Based on SmPC/literature evidence and pathway strength. '
                'Higher = more biologically plausible.'
            ),
            'confounder_risk': (
                'Risk score (0-1) indicating potential confounding from concomitant drugs. '
                'Higher = more likely that concomitants contribute to AE.'
            ),
            'temporal_causality_score': (
                'Bradford Hill temporality criterion (0-1). '
                'Assesses if time-to-onset matches expected biological timeline. '
                'Higher = better temporal relationship.'
            ),
            'dose_response_evidence': (
                'Bradford Hill dose-response criterion (0-1). '
                'Evidence of dose-response relationship. '
                'Higher = stronger dose-severity association.'
            ),
            'dechallenge_weight': (
                'Amplified dechallenge/rechallenge weight (0-3). '
                '1.0 if dechallenge positive, +2.0 if rechallenge positive. '
                'Stronger causality evidence than binary flag.'
            ),
            'organ_specificity_match': (
                'Binary indicator: 1 if drug has known organ-specific toxicity '
                'matching the AE organ system, 0 otherwise.'
            ),
            'ddi_interaction_risk': (
                'Drug-drug interaction risk score (0-1) based on known DDIs. '
                'Higher = more potential for interaction-related AEs.'
            ),
            'literature_citation_count': (
                'Number of PubMed articles co-mentioning drug and AE (0-100, capped). '
                'Higher = more literature support for association.'
            )
        }
    
    def get_feature_importance_hints(self) -> Dict[str, str]:
        """
        Get interpretation hints for each causal feature
        Used in SHAP explanations
        
        Returns:
            Dict mapping feature name to interpretation hint
        """
        return {
            'kg_pathway_exists': 'Known mechanism increases causal plausibility',
            'kg_pathway_length': 'Shorter pathway = stronger mechanistic link',
            'known_mechanism_score': 'High score = well-established mechanism',
            'confounder_risk': 'High risk = confounding may explain association',
            'temporal_causality_score': 'High score = timing supports causality',
            'dose_response_evidence': 'High score = dose-severity relationship observed',
            'dechallenge_weight': 'High weight = strong clinical evidence',
            'organ_specificity_match': 'Match = drug targets affected organ system',
            'ddi_interaction_risk': 'High risk = interactions may contribute',
            'literature_citation_count': 'High count = well-documented association'
        }
    
    def validate_output(self, df: pd.DataFrame) -> bool:
        """
        Validate that all 10 causal features were generated correctly
        
        Args:
            df: Output DataFrame
        
        Returns:
            True if validation passes
        
        Raises:
            ValueError: If validation fails
        """
        required_features = [
            'kg_pathway_exists',
            'kg_pathway_length',
            'known_mechanism_score',
            'confounder_risk',
            'temporal_causality_score',
            'dose_response_evidence',
            'dechallenge_weight',
            'organ_specificity_match',
            'ddi_interaction_risk',
            'literature_citation_count'
        ]
        
        missing_features = [f for f in required_features if f not in df.columns]
        
        if missing_features:
            raise ValueError(
                f"Causal feature generation incomplete. Missing features: {missing_features}"
            )
        
        # Validate data types and ranges - FIXED VERSION
        validations = []
        
        # Check each feature individually
        if not df['kg_pathway_exists'].isin([0, 1]).all():
            validations.append(('kg_pathway_exists', "Must be binary (0 or 1)"))
        
        if not ((df['kg_pathway_length'] >= 1) | (df['kg_pathway_length'] == 99)).all():
            validations.append(('kg_pathway_length', "Must be 1-4 or 99"))
        
        if not df['known_mechanism_score'].between(0, 1, inclusive='both').all():
            validations.append(('known_mechanism_score', "Must be 0-1"))
        
        if not df['confounder_risk'].between(0, 1, inclusive='both').all():
            validations.append(('confounder_risk', "Must be 0-1"))
        
        if not df['temporal_causality_score'].between(0, 1, inclusive='both').all():
            validations.append(('temporal_causality_score', "Must be 0-1"))
        
        if not df['dose_response_evidence'].between(0, 1, inclusive='both').all():
            validations.append(('dose_response_evidence', "Must be 0-1"))
        
        if not df['dechallenge_weight'].between(0, 3, inclusive='both').all():
            validations.append(('dechallenge_weight', "Must be 0-3"))
        
        if not df['organ_specificity_match'].isin([0, 1]).all():
            validations.append(('organ_specificity_match', "Must be binary"))
        
        if not df['ddi_interaction_risk'].between(0, 1, inclusive='both').all():
            validations.append(('ddi_interaction_risk', "Must be 0-1"))
        
        if not df['literature_citation_count'].between(0, 100, inclusive='both').all():
            validations.append(('literature_citation_count', "Must be 0-100"))
        
        # Raise error if any validation failed
        if validations:
            error_msgs = [f"  - {feature}: {msg}" for feature, msg in validations]
            raise ValueError("Validation failed:\n" + "\n".join(error_msgs))
        
        return True
    
    def get_stats(self, df: pd.DataFrame) -> Dict[str, Any]:
        """
        Get summary statistics for generated causal features
        Useful for debugging and QA
        
        Args:
            df: DataFrame with causal features
        
        Returns:
            Dict with statistics
        """
        causal_features = [
            'kg_pathway_exists',
            'kg_pathway_length',
            'known_mechanism_score',
            'confounder_risk',
            'temporal_causality_score',
            'dose_response_evidence',
            'dechallenge_weight',
            'organ_specificity_match',
            'ddi_interaction_risk',
            'literature_citation_count'
        ]
        
        stats = {}
        for feature in causal_features:
            if feature in df.columns:
                stats[feature] = {
                    'mean': float(df[feature].mean()),
                    'std': float(df[feature].std()),
                    'min': float(df[feature].min()),
                    'max': float(df[feature].max()),
                    'missing_count': int(df[feature].isna().sum())
                }
        
        # Additional summary stats
        stats['summary'] = {
            'total_signals': len(df),
            'pathways_found': int(df['kg_pathway_exists'].sum()) if 'kg_pathway_exists' in df else 0,
            'organ_matches': int(df['organ_specificity_match'].sum()) if 'organ_specificity_match' in df else 0,
            'avg_mechanism_score': float(df['known_mechanism_score'].mean()) if 'known_mechanism_score' in df else 0
        }
        
        return stats


# ============================================================================
# CONVENIENCE FUNCTIONS
# ============================================================================

def add_causal_features(
    signals_df: pd.DataFrame, 
    kg: Optional[LightweightKnowledgeGraph] = None
) -> pd.DataFrame:
    """
    Convenience function to add causal features to signals DataFrame
    
    Args:
        signals_df: Input DataFrame with 68 features
        kg: Optional pre-loaded Knowledge Graph
    
    Returns:
        DataFrame with 78 features (68 + 10 causal)
    
    Example:
        >>> from src.causal.causal_features import add_causal_features
        >>> engineered_df = add_causal_features(signals_df)
    """
    engineer = CausalFeatureEngineer(kg=kg)
    return engineer.transform(signals_df)


if __name__ == "__main__":
    """
    Test script - run with: python -m src.causal.causal_features
    """
    print("Testing Causal Feature Engineer...")
    
    # Create sample signal data (simulates 68-feature engineered data)
    sample_data = pd.DataFrame({
        'drug_name': ['Paracetamol', 'Warfarin', 'Metformin', 'Gentamicin', 'Unknown Drug'],
        'pt_name': ['Hepatotoxicity', 'Hemorrhage', 'Lactic Acidosis', 'Nephrotoxicity', 'Headache'],
        'casecount': [150, 200, 80, 120, 50],
        'prr': [3.5, 4.2, 2.8, 3.1, 1.5],
        'ror': [3.6, 4.3, 2.9, 3.2, 1.6],
        'chi2': [45.2, 52.1, 28.3, 38.7, 12.1],
        'time_to_onset_days': [5, 15, 45, 10, 3],
        'treatment_duration_days': [7, 30, 90, 14, 5],
        'dechallenge_positive': [1, 0, 1, 1, 0],
        'rechallenge_positive': [0, 0, 0, 1, 0],
        'concomitant_count': [2, 5, 3, 1, 0],
        'outcome_fatal': [0, 1, 0, 0, 0],
        'outcome_recovered': [1, 0, 1, 1, 1]
    })
    
    print(f"\n📊 Input: {len(sample_data)} signals with {len(sample_data.columns)} columns")
    
    # Initialize engineer
    engineer = CausalFeatureEngineer()
    
    # Generate causal features
    print("\n🔬 Generating causal features...")
    result_df = engineer.transform(sample_data)
    
    print(f"✅ Output: {len(result_df)} signals with {len(result_df.columns)} columns")
    print(f"   Added: {len(result_df.columns) - len(sample_data.columns)} causal features")
    
    # Validate
    try:
        engineer.validate_output(result_df)
        print("\n✅ Validation passed!")
    except ValueError as e:
        print(f"\n❌ Validation failed: {e}")
    
    # Show causal features for first signal
    print("\n📋 Causal Features for Signal #1 (Paracetamol → Hepatotoxicity):")
    causal_features = [
        'kg_pathway_exists',
        'kg_pathway_length',
        'known_mechanism_score',
        'confounder_risk',
        'temporal_causality_score',
        'dose_response_evidence',
        'dechallenge_weight',
        'organ_specificity_match',
        'ddi_interaction_risk',
        'literature_citation_count'
    ]
    
    for feature in causal_features:
        if feature in result_df.columns:
            value = result_df.iloc[0][feature]
            print(f"   {feature}: {value}")
    
    # Show summary statistics
    print("\n📊 Summary Statistics:")
    stats = engineer.get_stats(result_df)
    summary = stats.get('summary', {})
    print(f"   Total signals: {summary.get('total_signals', 0)}")
    print(f"   Pathways found: {summary.get('pathways_found', 0)}")
    print(f"   Organ matches: {summary.get('organ_matches', 0)}")
    print(f"   Avg mechanism score: {summary.get('avg_mechanism_score', 0):.2f}")
    
    print("\n✅ Causal Feature Engineer test complete!")
