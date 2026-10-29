"""
BERT-Based Semantic Column Mapper
Implements Traffic Light Validator (Green/Yellow/Red zones)
As described in LinkedIn post on PV canonicalization
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional
from sentence_transformers import SentenceTransformer
from sklearn.metrics.pairwise import cosine_similarity
from pathlib import Path
import json
from src.utils.logger import setup_logger

logger = setup_logger(__name__)


class BERTSemanticMapper:
    """
    BERT-based semantic column mapper with Traffic Light validation.
    
    Traffic Light Zones:
    🟢 Green (>0.95): Auto-map
    🟡 Yellow (0.75-0.95): Human review
    🔴 Red (<0.75): Reject
    """
    
    def __init__(self, model_name: str = 'all-MiniLM-L6-v2'):
        """
        Initialize BERT semantic mapper.
        
        Args:
            model_name: SentenceTransformer model (lightweight & fast)
        """
        logger.info(f"Loading BERT model: {model_name}")
        self.model = SentenceTransformer(model_name)
        self.golden_schema = self._load_golden_schema()
        self.golden_embeddings = self._compute_golden_embeddings()
        
        # Traffic light thresholds
        self.GREEN_THRESHOLD = 0.95
        self.YELLOW_THRESHOLD = 0.75
        
    def _load_golden_schema(self) -> Dict[str, str]:
        """
        Load 36-column golden schema with descriptions.
        
        Returns:
            Dict mapping column name to semantic description
        """
        return {
            'case_id': 'unique case identifier primary ID patient identifier',
            'receive_date': 'date report received FDA submission date',
            'country': 'country of occurrence reporter country',
            'age': 'patient age in years',
            'age_group': 'age category pediatric adult elderly',
            'sex': 'patient sex gender male female',
            'weight_kg': 'patient weight in kilograms body weight',
            'suspect_drug': 'suspect drug name medication drugname',
            'indication': 'drug indication reason for use',
            'route': 'route of administration oral intravenous',
            'dose': 'dose amount quantity',
            'dose_unit': 'dose unit mg mcg',
            'dose_frequency': 'dosing frequency QD BID',
            'treatment_duration_days': 'treatment duration therapy duration days',
            'event_llt': 'event lowest level term LLT',
            'event_pt': 'event preferred term PT reaction',
            'event_hlt': 'event high level term HLT',
            'event_hlgt': 'event high level group term HLGT',
            'event_soc': 'event system organ class SOC',
            'time_to_onset_days': 'time to onset days latency',
            'event_duration_days': 'event duration days',
            'seriousness': 'seriousness serious criteria death hospitalization',
            'outcome': 'outcome recovered fatal',
            'action_taken': 'action taken drug withdrawn dose reduced',
            'dechallenge': 'dechallenge positive negative',
            'rechallenge': 'rechallenge positive negative',
            'concomitant_medications': 'concomitant medications concomitant drugs',
            'medical_history': 'medical history past medical conditions',
            'reporter_type': 'reporter type qualification physician pharmacist consumer',
            'report_type': 'report type spontaneous clinical trial',
            'causality_assessment': 'causality assessment WHO UMC certain probable possible',
            'alt_u_l': 'ALT alanine aminotransferase liver enzyme U/L',
            'ast_u_l': 'AST aspartate aminotransferase liver enzyme U/L',
            'bilirubin_mg_dl': 'bilirubin total bilirubin mg/dL',
            'creatinine_mg_dl': 'creatinine serum creatinine mg/dL kidney',
            'bun_mg_dl': 'BUN blood urea nitrogen mg/dL kidney'
        }
    
    def _compute_golden_embeddings(self) -> np.ndarray:
        """Compute BERT embeddings for golden schema."""
        logger.info("Computing golden schema embeddings...")
        
        descriptions = list(self.golden_schema.values())
        embeddings = self.model.encode(descriptions, show_progress_bar=False)
        
        logger.info(f"✅ Computed embeddings for {len(embeddings)} golden columns")
        return embeddings
    
    def map_columns(
        self, 
        source_columns: List[str],
        source_name: str = "unknown"
    ) -> Dict[str, Tuple[str, float, str]]:
        """
        Map source columns to golden schema using semantic similarity.
        
        Args:
            source_columns: List of column names from source data
            source_name: Name of data source (for logging)
            
        Returns:
            Dict mapping source_col -> (golden_col, score, zone)
            zone: 'green', 'yellow', or 'red'
        """
        logger.info(f"Mapping {len(source_columns)} columns from {source_name}")
        
        # Compute embeddings for source columns
        source_embeddings = self.model.encode(source_columns, show_progress_bar=False)
        
        # Compute cosine similarity
        similarity_matrix = cosine_similarity(source_embeddings, self.golden_embeddings)
        
        # Map each source column to best match
        mapping = {}
        golden_cols = list(self.golden_schema.keys())
        
        for i, source_col in enumerate(source_columns):
            scores = similarity_matrix[i]
            best_idx = np.argmax(scores)
            best_score = scores[best_idx]
            golden_col = golden_cols[best_idx]
            
            # Determine traffic light zone
            if best_score >= self.GREEN_THRESHOLD:
                zone = 'green'
            elif best_score >= self.YELLOW_THRESHOLD:
                zone = 'yellow'
            else:
                zone = 'red'
            
            mapping[source_col] = (golden_col, best_score, zone)
        
        # Log statistics
        green = sum(1 for _, (_, _, z) in mapping.items() if z == 'green')
        yellow = sum(1 for _, (_, _, z) in mapping.items() if z == 'yellow')
        red = sum(1 for _, (_, _, z) in mapping.items() if z == 'red')
        
        logger.info(f"Mapping complete: 🟢 {green} Green | 🟡 {yellow} Yellow | 🔴 {red} Red")
        
        return mapping
    
    def apply_mapping(
        self,
        df: pd.DataFrame,
        mapping: Dict[str, Tuple[str, float, str]],
        auto_approve_green: bool = True,
        reject_red: bool = True
    ) -> pd.DataFrame:
        """
        Apply semantic mapping to dataframe.
        
        Args:
            df: Source dataframe
            mapping: Mapping from map_columns()
            auto_approve_green: Auto-map green zone columns
            reject_red: Drop red zone columns
            
        Returns:
            DataFrame with golden schema column names
        """
        canonical = pd.DataFrame()
        
        yellow_queue = []
        red_rejected = []
        
        for source_col, (golden_col, score, zone) in mapping.items():
            if zone == 'green' and auto_approve_green:
                canonical[golden_col] = df[source_col]
                logger.debug(f"🟢 Auto-mapped: {source_col} → {golden_col} ({score:.3f})")
                
            elif zone == 'yellow':
                yellow_queue.append((source_col, golden_col, score))
                logger.warning(f"🟡 Human review needed: {source_col} → {golden_col} ({score:.3f})")
                # For now, map anyway (production should queue for review)
                canonical[golden_col] = df[source_col]
                
            elif zone == 'red':
                red_rejected.append((source_col, score))
                if not reject_red:
                    # Keep in original column name
                    canonical[source_col] = df[source_col]
                logger.error(f"🔴 Rejected: {source_col} (score: {score:.3f})")
        
        # Log summary
        if yellow_queue:
            logger.warning(f"⚠️ {len(yellow_queue)} columns need human review")
        if red_rejected:
            logger.error(f"❌ {len(red_rejected)} columns rejected due to low similarity")
        
        return canonical
    
    def save_mapping(self, mapping: Dict, output_path: Path):
        """Save mapping for audit trail."""
        mapping_serializable = {
            source: {
                'golden_column': golden,
                'similarity_score': float(score),
                'zone': zone
            }
            for source, (golden, score, zone) in mapping.items()
        }
        
        with open(output_path, 'w') as f:
            json.dump(mapping_serializable, f, indent=2)
        
        logger.info(f"✅ Mapping saved to {output_path}")
