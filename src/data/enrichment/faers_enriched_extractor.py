"""
FAERS Enriched Extractor
Extracts and normalizes FAERS ASCII datasets into the 36-column canonical ICSR schema.
Preserves all reported adverse drug reactions (no .first() truncation) and extracts
concomitant medications from non-primary suspect drug records.

Legal Attribution:
- US Food and Drug Administration (FDA) Adverse Event Reporting System (FAERS).
- Public domain US Government data under 21 CFR Part 314.80 / 312.32.
"""

import zipfile
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional
from src.utils.logger import setup_logger

logger = setup_logger(__name__)


class FAERSEnrichedExtractor:
    """Extracts 36-column golden schema from raw FAERS ASCII quarterly zip files."""

    def __init__(self):
        self.required_tables = ['DEMO', 'DRUG', 'REAC', 'OUTC', 'RPSR', 'THER', 'INDI']

    def extract_from_zip(self, zip_path: Path) -> Optional[pd.DataFrame]:
        """
        Extract all 7 tables from FAERS zip and merge into canonical schema.

        Args:
            zip_path: Path to FAERS quarterly zip file.

        Returns:
            DataFrame with golden schema columns, or None if extraction fails.
        """
        logger.info(f"Extracting from {zip_path.name}...")

        try:
            with zipfile.ZipFile(zip_path, 'r') as zf:
                tables = self._find_table_files(zf)

                missing = [t for t in self.required_tables if t not in tables]
                if missing:
                    logger.warning(f"Missing tables in zip: {missing}")

                # Read tables
                demo_df = self._read_table(zf, tables.get('DEMO'))
                drug_df = self._read_table(zf, tables.get('DRUG'))
                reac_df = self._read_table(zf, tables.get('REAC'))
                outc_df = self._read_table(zf, tables.get('OUTC'))
                rpsr_df = self._read_table(zf, tables.get('RPSR'))
                ther_df = self._read_table(zf, tables.get('THER'))
                indi_df = self._read_table(zf, tables.get('INDI'))

                # Merge tables preserving all reactions and extracting concomitants
                merged = self._merge_tables(demo_df, drug_df, reac_df, outc_df, rpsr_df, ther_df, indi_df)

                # Map to golden schema
                canonical = self._map_to_golden_schema(merged)

                logger.info(f"✅ Extracted {len(canonical):,} ICSR drug-event records")
                return canonical

        except Exception as e:
            logger.error(f"Extraction failed: {e}", exc_info=True)
            return None

    def _find_table_files(self, zf: zipfile.ZipFile) -> Dict[str, str]:
        """Find all 7 FAERS table files inside zip."""
        tables = {}
        for name in zf.namelist():
            name_lower = name.lower()
            for table in self.required_tables:
                if table.lower() in name_lower and name_lower.endswith('.txt'):
                    tables[table] = name
                    break
        logger.info(f"Found tables: {list(tables.keys())}")
        return tables

    def _read_table(self, zf: zipfile.ZipFile, filename: Optional[str]) -> Optional[pd.DataFrame]:
        """Read a single FAERS table from zip file."""
        if not filename:
            return None

        try:
            with zf.open(filename) as f:
                df = pd.read_csv(f, sep='$', encoding='latin-1', on_bad_lines='skip', low_memory=False)
            df.columns = df.columns.str.lower().str.strip()
            return df
        except Exception as e:
            logger.error(f"Failed to read {filename}: {e}")
            return None

    def _merge_tables(self, demo_df, drug_df, reac_df, outc_df, rpsr_df, ther_df, indi_df) -> pd.DataFrame:
        """
        Merge all FAERS tables on primaryid.
        Preserves all reported adverse reactions per case and gathers concomitant drugs.
        """
        if drug_df is None:
            raise ValueError("DRUG table required")

        # 1. Identify primary/secondary suspect drugs and separate concomitant medications
        if 'role_cod' in drug_df.columns:
            suspect_mask = drug_df['role_cod'].astype(str).str.upper().isin(['PS', 'SS'])
            if not suspect_mask.any():
                suspect_drugs = drug_df.copy()
                concom_series = pd.Series(dtype=str)
            else:
                suspect_drugs = drug_df[suspect_mask].copy()
                concom_drugs = drug_df[~suspect_mask].copy()
                # Aggregate concomitant drug names per case
                drug_name_col = 'drugname' if 'drugname' in concom_drugs.columns else ('prod_ai' if 'prod_ai' in concom_drugs.columns else None)
                if drug_name_col:
                    concom_series = concom_drugs.groupby('primaryid')[drug_name_col].apply(
                        lambda s: ", ".join(s.dropna().astype(str).str.upper().unique())
                    )
                else:
                    concom_series = pd.Series(dtype=str)
        else:
            suspect_drugs = drug_df.copy()
            concom_series = pd.Series(dtype=str)

        # 2. Preserve ALL adverse events by merging suspect drugs with reac_df
        # DO NOT truncate using .groupby('primaryid').first()
        if reac_df is not None:
            event_col = 'pt' if 'pt' in reac_df.columns else ('reaction' if 'reaction' in reac_df.columns else None)
            if event_col:
                merged = suspect_drugs.merge(reac_df[['primaryid', event_col]], on='primaryid', how='inner')
                if event_col != 'pt':
                    merged['pt'] = merged[event_col]
            else:
                merged = suspect_drugs.copy()
                merged['pt'] = np.nan
        else:
            merged = suspect_drugs.copy()
            merged['pt'] = np.nan

        # 3. Attach concomitant medications
        if not concom_series.empty:
            merged = merged.merge(concom_series.rename('concomitant_medications'), on='primaryid', how='left')
        else:
            merged['concomitant_medications'] = ''

        # 4. Merge DEMO (demographics)
        if demo_df is not None:
            demo_cols = [c for c in demo_df.columns if c not in merged.columns or c == 'primaryid']
            merged = merged.merge(demo_df[demo_cols], on='primaryid', how='left', suffixes=('', '_demo'))

        # 5. Merge OUTC (patient outcomes - concatenate unique outcome codes per case)
        if outc_df is not None and 'outc_cod' in outc_df.columns:
            outc_agg = outc_df.groupby('primaryid')['outc_cod'].apply(
                lambda s: "; ".join(s.dropna().astype(str).str.upper().unique())
            ).reset_index()
            merged = merged.merge(outc_agg, on='primaryid', how='left', suffixes=('', '_outc'))

        # 6. Merge RPSR (reporter)
        if rpsr_df is not None and 'rpsr_cod' in rpsr_df.columns:
            rpsr_first = rpsr_df.groupby('primaryid').first().reset_index()
            merged = merged.merge(rpsr_first, on='primaryid', how='left', suffixes=('', '_rpsr'))

        # 7. Merge THER (therapy dates)
        if ther_df is not None and 'drug_seq' in ther_df.columns and 'drug_seq' in merged.columns:
            ther_first = ther_df.groupby(['primaryid', 'drug_seq']).first().reset_index()
            merged = merged.merge(ther_first, on=['primaryid', 'drug_seq'], how='left', suffixes=('', '_ther'))

        # 8. Merge INDI (indications)
        if indi_df is not None and 'drug_seq' in indi_df.columns and 'drug_seq' in merged.columns:
            indi_first = indi_df.groupby(['primaryid', 'drug_seq']).first().reset_index()
            merged = merged.merge(indi_first, on=['primaryid', 'drug_seq'], how='left', suffixes=('', '_indi'))

        return merged

    def _map_to_golden_schema(self, faers_df: pd.DataFrame) -> pd.DataFrame:
        """Map merged FAERS dataframe to the 36-column canonical ICSR schema."""
        c = pd.DataFrame()

        # Basic identification
        c['case_id'] = faers_df.get('primaryid', '').astype(str)
        c['receive_date'] = pd.to_datetime(faers_df.get('fda_dt'), errors='coerce')
        c['country'] = faers_df.get('reporter_country', faers_df.get('occr_country', '')).fillna('')

        # Demographics
        c['age'] = pd.to_numeric(faers_df.get('age'), errors='coerce')
        c['age_group'] = c['age'].apply(self._categorize_age)
        c['sex'] = faers_df.get('sex', 'UNK').astype(str).str.upper().map({'M': 'Male', 'F': 'Female'}).fillna('Unknown')
        c['weight_kg'] = pd.to_numeric(faers_df.get('wt'), errors='coerce')

        # Drug information
        drug_name_col = 'drugname' if 'drugname' in faers_df.columns else ('prod_ai' if 'prod_ai' in faers_df.columns else None)
        c['suspect_drug'] = faers_df.get(drug_name_col, '').astype(str).str.upper().str.strip() if drug_name_col else ''
        c['indication'] = faers_df.get('indi_pt', '').fillna('')
        c['route'] = faers_df.get('route', '').fillna('')
        c['dose'] = pd.to_numeric(faers_df.get('dose_amt'), errors='coerce')
        c['dose_unit'] = faers_df.get('dose_unit', '').fillna('')
        c['dose_frequency'] = faers_df.get('dose_freq', '').fillna('')
        c['treatment_duration_days'] = np.nan

        # Adverse Event information
        c['event_llt'] = faers_df.get('pt', '').fillna('')
        c['event_pt'] = faers_df.get('pt', '').astype(str).str.upper().str.strip()
        c['event_hlt'] = ''
        c['event_hlgt'] = ''
        c['event_soc'] = ''
        c['time_to_onset_days'] = np.nan
        c['event_duration_days'] = np.nan

        # Outcomes and Severity
        c['seriousness'] = faers_df.get('serious', '').fillna('')
        c['outcome'] = faers_df.get('outc_cod', '').fillna('')
        c['action_taken'] = ''
        c['dechallenge'] = ''
        c['rechallenge'] = ''

        # Clinical context (populated from non-primary suspect records)
        c['concomitant_medications'] = faers_df.get('concomitant_medications', '').fillna('')
        c['medical_history'] = ''

        # Reporter
        c['reporter_type'] = faers_df.get('rpsr_cod', faers_df.get('reporter_qualification', '')).fillna('')
        c['report_type'] = 'Spontaneous'

        # Regulatory & Lab markers
        c['causality_assessment'] = 'Unassessed'
        c['alt_u_l'] = np.nan
        c['ast_u_l'] = np.nan
        c['bilirubin_mg_dl'] = np.nan
        c['creatinine_mg_dl'] = np.nan
        c['bun_mg_dl'] = np.nan

        return c

    def _categorize_age(self, age):
        """Categorize age into regulatory demographic groups."""
        if pd.isna(age):
            return 'Unknown'
        if age < 0.08:
            return 'neonate'
        elif age < 2:
            return 'infant'
        elif age < 12:
            return 'child'
        elif age < 18:
            return 'adolescent'
        elif age < 65:
            return 'adult'
        else:
            return 'elderly'
