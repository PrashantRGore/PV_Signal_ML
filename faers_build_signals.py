"""
FAERS Signal Builder v3.0 - NOW WITH 36-COLUMN EXTRACTION
Integrates: FAERSEnrichedExtractor + FeatureEngineering + DAE + BERT Mapper
"""

import pandas as pd
from pathlib import Path
from typing import Optional
import config
from src.data.enrichment.faers_enriched_extractor import FAERSEnrichedExtractor
from src.data.enrichment.bert_semantic_mapper import BERTSemanticMapper
from src.ml.feature_engineering.feature_pipeline import FeatureEngineeringPipeline
from src.ml.nan_handler.denoising_autoencoder import DenoisingAutoencoder
from src.stats_engine.disproportionality import DisproportionalityAnalysis
from src.utils.logger import setup_logger

logger = setup_logger(__name__)


def build_faers_signals_v3(
    quarter: str = "24Q3",
    use_enriched: bool = True,
    use_feature_engineering: bool = True,
    use_bert_mapping: bool = True,
    use_dae: bool = True
) -> Optional[pd.DataFrame]:
    """
    Build FAERS signals with full 36-column extraction pipeline.
    
    Args:
        quarter: FAERS quarter (e.g., "24Q3")
        use_enriched: Use 36-column enriched extractor (vs old 7-column)
        use_feature_engineering: Apply feature engineering (36 → 50+)
        use_bert_mapping: Use BERT semantic mapper for column alignment
        use_dae: Apply DAE for missing data imputation
        
    Returns:
        DataFrame with signals (36-50+ columns depending on options)
    """
    logger.info(f"Building FAERS signals v3.0 for {quarter}")
    logger.info(f"Options: enriched={use_enriched}, features={use_feature_engineering}, "
                f"bert={use_bert_mapping}, dae={use_dae}")
    
    # Step 1: Download FAERS data
    from src.data.faers_downloader import FAERSDownloader
    downloader = FAERSDownloader()
    zip_path = downloader.download_quarter(quarter)
    
    if not zip_path:
        logger.error(f"Failed to download FAERS {quarter}")
        return None
    
    # Step 2: Extract with enriched extractor (36 columns) or legacy (7 columns)
    if use_enriched:
        logger.info("Using FAERSEnrichedExtractor (36 columns)")
        extractor = FAERSEnrichedExtractor()
        icsr_df = extractor.extract_from_zip(zip_path, quarter)
        
        if icsr_df is None or len(icsr_df) == 0:
            logger.error("Enriched extraction failed")
            return None
        
        logger.info(f"✅ Extracted {len(icsr_df):,} ICSRs with {len(icsr_df.columns)} columns")
        
        # Step 3: BERT semantic mapping (align to golden schema)
        if use_bert_mapping:
            logger.info("Applying BERT semantic mapping...")
            bert_mapper = BERTSemanticMapper()
            mapping = bert_mapper.map_columns(
                list(icsr_df.columns),
                source_name=f"FAERS_{quarter}"
            )
            
            # Save mapping for audit
            mapping_path = Path(config.DATA_DIR) / f"faers_{quarter}_mapping.json"
            bert_mapper.save_mapping(mapping, mapping_path)
            
            # Apply mapping
            icsr_df = bert_mapper.apply_mapping(
                icsr_df,
                mapping,
                auto_approve_green=True,
                reject_red=True
            )
            
            logger.info(f"✅ BERT mapping complete: {len(icsr_df.columns)} canonical columns")
    
    else:
        # Legacy 7-column extraction
        logger.warning("Using legacy 7-column extraction (not recommended)")
        from src.data.faers_build_signals import extract_faers_legacy
        icsr_df = extract_faers_legacy(zip_path)
    
    # Step 4: Compute disproportionality signals
    logger.info("Computing disproportionality signals...")
    dispro = DisproportionalityAnalysis(min_cases=3)
    signals_df = dispro.compute_signals(icsr_df)
    
    if len(signals_df) == 0:
        logger.warning("No signals computed")
        return None
    
    logger.info(f"✅ Computed {len(signals_df):,} signals")
    
    # Step 5: Feature engineering (36 → 50+)
    if use_feature_engineering:
        logger.info("Applying feature engineering pipeline...")
        feature_pipeline = FeatureEngineeringPipeline()
        
        # Merge signals with ICSR data to get clinical features
        signals_enriched = signals_df.merge(
            icsr_df,
            left_on=['drug_name', 'event_term'],
            right_on=['suspect_drug', 'event_pt'],
            how='left'
        )
        
        # Engineer features
        signals_enriched = feature_pipeline.transform(signals_enriched)
        
        logger.info(f"✅ Feature engineering complete: {len(signals_enriched.columns)} features")
        signals_df = signals_enriched
    
    # Step 6: DAE imputation (remove NaNs before ML/Ollama)
    if use_dae:
        logger.info("Applying DAE for missing data imputation...")
        dae = DenoisingAutoencoder()
        
        # Load pre-trained DAE or fit on current data
        if not dae.load():
            logger.warning("No pre-trained DAE found. Fitting on current data...")
            # Fit on non-null subset
            clean_subset = signals_df.select_dtypes(include=['number']).dropna()
            if len(clean_subset) > 100:
                dae.fit(clean_subset.sample(min(10000, len(clean_subset))), 
                       validation_dataset=f"FAERS_{quarter}")
                dae.save()
        
        # Impute missing values
        numeric_cols = signals_df.select_dtypes(include=['number']).columns
        signals_df[numeric_cols] = dae.transform(signals_df[numeric_cols])
        
        logger.info("✅ DAE imputation complete: 0 NaNs remaining")
    
    # Step 7: Save to disk
    output_path = Path(config.DATA_DIR) / f"faers_signals_{quarter}_v3.parquet"
    signals_df.to_parquet(output_path, index=False)
    
    logger.info(f"✅ Saved signals to {output_path}")
    logger.info(f"Final dataset: {len(signals_df):,} signals × {len(signals_df.columns)} columns")
    
    return signals_df


def extract_faers_legacy(zip_path: Path) -> pd.DataFrame:
    """Legacy 7-column extraction (for backwards compatibility)."""
    import zipfile
    
    with zipfile.ZipFile(zip_path, 'r') as zf:
        # Find DRUG and REAC tables
        drug_file = [f for f in zf.namelist() if 'drug' in f.lower() and f.endswith('.txt')][0]
        reac_file = [f for f in zf.namelist() if 'reac' in f.lower() and f.endswith('.txt')][0]
        
        # Read tables
        drug_df = pd.read_csv(zf.open(drug_file), sep='$', encoding='latin-1', on_bad_lines='skip')
        reac_df = pd.read_csv(zf.open(reac_file), sep='$', encoding='latin-1', on_bad_lines='skip')
    
    # Normalize columns
    drug_df.columns = drug_df.columns.str.lower()
    reac_df.columns = reac_df.columns.str.lower()
    
    # Create ICSR dataframe
    icsr_df = drug_df[['primaryid', 'drugname']].copy()
    icsr_df = icsr_df.merge(reac_df[['primaryid', 'pt']], on='primaryid', how='inner')
    
    # Rename to canonical
    icsr_df.rename(columns={
        'primaryid': 'case_id',
        'drugname': 'drug_name',
        'pt': 'event_term'
    }, inplace=True)
    
    return icsr_df


if __name__ == "__main__":
    # Test full pipeline
    signals = build_faers_signals_v3(
        quarter="24Q3",
        use_enriched=True,
        use_feature_engineering=True,
        use_bert_mapping=True,
        use_dae=True
    )
    
    if signals is not None:
        print(f"\n✅ SUCCESS: Built {len(signals):,} signals")
        print(f"Columns ({len(signals.columns)}): {list(signals.columns[:10])}...")
        print(f"\nSample signal:")
        print(signals.head(1).T)
