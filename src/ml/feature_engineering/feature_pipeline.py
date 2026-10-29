"""
Feature Engineering Pipeline - 36 Columns → 50+ Features
Generates derived features for qualitative signal analysis

STRICT VERSION: No silent defaults, explicit NaN + data quality flags
"""

import pandas as pd
import numpy as np
from typing import List
from src.utils.logger import setup_logger

logger = setup_logger(__name__)

class FeatureEngineeringPipeline:
    """Engineer 50+ features from 36-column golden schema with strict column handling."""
    
    def __init__(self):
        self.feature_names = []
    
    def transform(self, df: pd.DataFrame) -> pd.DataFrame:
        logger.info(f"Engineering features from {len(df.columns)} base columns")
        engineered = df.copy()
        
        engineered = self._engineer_demographics(engineered)
        engineered = self._engineer_temporal(engineered)
        engineered = self._engineer_clinical(engineered)
        engineered = self._engineer_lab_values(engineered)
        engineered = self._engineer_causality(engineered)
        engineered = self._engineer_polypharmacy(engineered)
        engineered = self._engineer_geographic(engineered)
        
        logger.info(f"✅ Engineered {len(engineered.columns)} total features")
        return engineered
    
    def _engineer_demographics(self, df: pd.DataFrame) -> pd.DataFrame:
        if 'age' not in df.columns:
            logger.critical("⚠️⚠️⚠️ MISSING COLUMN: 'age' - setting age features to NaN")
            df['age_bin_pediatric'] = np.nan
            df['age_bin_adult'] = np.nan
            df['age_bin_elderly'] = np.nan
            df['_data_quality_age_missing'] = 1
        else:
            df['age_bin_pediatric'] = (df['age'] < 18).astype(int)
            df['age_bin_adult'] = ((df['age'] >= 18) & (df['age'] < 65)).astype(int)
            df['age_bin_elderly'] = (df['age'] >= 65).astype(int)
            df['_data_quality_age_missing'] = 0
        
        weight_col = 'weight' if 'weight' in df.columns else ('weight_kg' if 'weight_kg' in df.columns else None)
        height_col = next((c for c in ['height', 'height_m', 'height_cm', 'height_in'] if c in df.columns), None)

        if not weight_col:
            logger.critical("⚠️⚠️⚠️ MISSING COLUMN: 'weight' - setting BMI features to NaN")
            df['bmi_estimated'] = np.nan
            df['bmi_underweight'] = np.nan
            df['bmi_normal'] = np.nan
            df['bmi_overweight'] = np.nan
            df['bmi_obese'] = np.nan
            df['_data_quality_weight_missing'] = 1
            df['_data_quality_height_imputed'] = 0
        else:
            df['_data_quality_weight_missing'] = 0
            # Calculate height in meters if available
            if height_col:
                h = pd.to_numeric(df[height_col], errors='coerce')
                # If values > 3, assume height is in cm
                h_m = np.where(h > 3.0, h / 100.0, h)
                height_missing = pd.isna(h_m) | (h_m <= 0)
                h_m_filled = np.where(height_missing, 1.7, h_m)
                df['bmi_estimated'] = df[weight_col] / (h_m_filled ** 2)
                df['_data_quality_height_imputed'] = height_missing.astype(int)
            else:
                # Height unavailable: use reference standard (1.7m) but flag as imputed
                df['bmi_estimated'] = df[weight_col] / (1.7 ** 2)
                df['_data_quality_height_imputed'] = 1

            df['bmi_underweight'] = (df['bmi_estimated'] < 18.5).astype(int)
            df['bmi_normal'] = ((df['bmi_estimated'] >= 18.5) & (df['bmi_estimated'] < 25)).astype(int)
            df['bmi_overweight'] = ((df['bmi_estimated'] >= 25) & (df['bmi_estimated'] < 30)).astype(int)
            df['bmi_obese'] = (df['bmi_estimated'] >= 30).astype(int)
        
        return df
    
    def _engineer_temporal(self, df: pd.DataFrame) -> pd.DataFrame:
        ttp_col = 'ttp_days' if 'ttp_days' in df.columns else ('time_to_onset_days' if 'time_to_onset_days' in df.columns else None)
        if not ttp_col:
            logger.critical("⚠️⚠️⚠️ MISSING COLUMN: 'time_to_onset' - setting TTP features to NaN")
            df['ttp_acute'] = np.nan
            df['ttp_subacute'] = np.nan
            df['ttp_chronic'] = np.nan
            df['ttp_missing'] = 1
        else:
            df['ttp_acute'] = (df[ttp_col] <= 7).astype(int)
            df['ttp_subacute'] = ((df[ttp_col] > 7) & (df[ttp_col] <= 90)).astype(int)
            df['ttp_chronic'] = (df[ttp_col] > 90).astype(int)
            df['ttp_missing'] = df[ttp_col].isna().astype(int)
        
        dur_col = 'therapy_duration' if 'therapy_duration' in df.columns else ('treatment_duration_days' if 'treatment_duration_days' in df.columns else None)
        if not dur_col:
            logger.critical("⚠️⚠️⚠️ MISSING COLUMN: 'treatment_duration' - setting duration features to NaN")
            df['treatment_short'] = np.nan
            df['treatment_long'] = np.nan
            df['_data_quality_duration_missing'] = 1
        else:
            df['treatment_short'] = (df[dur_col] <= 7).astype(int)
            df['treatment_long'] = (df[dur_col] > 90).astype(int)
            df['_data_quality_duration_missing'] = 0
        
        return df
    
    def _engineer_clinical(self, df: pd.DataFrame) -> pd.DataFrame:
        serious_col = 'serious' if 'serious' in df.columns else ('seriousness' if 'seriousness' in df.columns else None)
        if not serious_col:
            logger.critical("⚠️⚠️⚠️ MISSING COLUMN: 'seriousness' - setting seriousness features to NaN")
            df['serious_death'] = np.nan
            df['serious_hospitalization'] = np.nan
            df['serious_lifethreat'] = np.nan
            df['_data_quality_serious_missing'] = 1
        else:
            df['serious_death'] = df[serious_col].astype(str).str.contains('Death', case=False, na=False).astype(int)
            df['serious_hospitalization'] = df[serious_col].astype(str).str.contains('Hosp', case=False, na=False).astype(int)
            df['serious_lifethreat'] = df[serious_col].astype(str).str.contains('Life', case=False, na=False).astype(int)
            df['_data_quality_serious_missing'] = 0
        
        if 'outcome' not in df.columns:
            logger.critical("⚠️⚠️⚠️ MISSING COLUMN: 'outcome' - setting outcome features to NaN")
            df['outcome_recovered'] = np.nan
            df['outcome_fatal'] = np.nan
            df['outcome_unknown'] = np.nan
            df['_data_quality_outcome_missing'] = 1
        else:
            df['outcome_recovered'] = df['outcome'].astype(str).str.contains('Recover', case=False, na=False).astype(int)
            df['outcome_fatal'] = df['outcome'].astype(str).str.contains('Death|Fatal', case=False, na=False).astype(int)
            df['outcome_unknown'] = df['outcome'].astype(str).str.contains('Unknown', case=False, na=False).astype(int)
            df['_data_quality_outcome_missing'] = 0
        
        action_col = 'action_drug' if 'action_drug' in df.columns else ('action_taken' if 'action_taken' in df.columns else None)
        if not action_col:
            logger.critical("⚠️⚠️⚠️ MISSING COLUMN: 'action_taken' - setting action features to NaN")
            df['action_withdrawn'] = np.nan
            df['action_dose_reduced'] = np.nan
            df['action_no_change'] = np.nan
            df['_data_quality_action_missing'] = 1
        else:
            df['action_withdrawn'] = df[action_col].astype(str).str.contains('Withdrawn|Stopped', case=False, na=False).astype(int)
            df['action_dose_reduced'] = df[action_col].astype(str).str.contains('Reduced', case=False, na=False).astype(int)
            df['action_no_change'] = df[action_col].astype(str).str.contains('No change', case=False, na=False).astype(int)
            df['_data_quality_action_missing'] = 0
        
        return df
    
    def _engineer_lab_values(self, df: pd.DataFrame) -> pd.DataFrame:
        for col, hi, severe in [('alt', 40, 120), ('ast', 40, 120)]:
            if col not in df.columns:
                logger.critical(f"⚠️⚠️⚠️ MISSING COLUMN: '{col}' - setting lab features to NaN")
                df[f'{col}_elevated'] = np.nan
                df[f'{col}_severe'] = np.nan
            else:
                df[f'{col}_elevated'] = (df[col] > hi).astype(int)
                df[f'{col}_severe'] = (df[col] > severe).astype(int)
        
        if 'bilirubin' not in df.columns:
            logger.critical("⚠️⚠️⚠️ MISSING COLUMN: 'bilirubin'")
            df['bilirubin_elevated'] = np.nan
        else:
            df['bilirubin_elevated'] = (df['bilirubin'] > 1.2).astype(int)
        
        if 'creatinine' not in df.columns:
            logger.critical("⚠️⚠️⚠️ MISSING COLUMN: 'creatinine'")
            df['creatinine_elevated'] = np.nan
            df['creatinine_severe'] = np.nan
        else:
            df['creatinine_elevated'] = (df['creatinine'] > 1.2).astype(int)
            df['creatinine_severe'] = (df['creatinine'] > 2.0).astype(int)
        
        if 'bun' not in df.columns:
            logger.critical("⚠️⚠️⚠️ MISSING COLUMN: 'bun'")
            df['bun_elevated'] = np.nan
        else:
            df['bun_elevated'] = (df['bun'] > 20).astype(int)
        
        df['hepatic_injury_score'] = df[['alt_elevated', 'ast_elevated', 'bilirubin_elevated']].sum(axis=1)
        df['renal_injury_score'] = df[['creatinine_elevated', 'bun_elevated']].sum(axis=1)
        return df
    
    def _engineer_causality(self, df: pd.DataFrame) -> pd.DataFrame:
        causality_col = 'causality' if 'causality' in df.columns else ('causality_assessment' if 'causality_assessment' in df.columns else None)
        if not causality_col:
            logger.critical("⚠️⚠️⚠️ MISSING COLUMN: 'causality'")
            df['causality_certain'] = np.nan
            df['causality_probable'] = np.nan
            df['causality_possible'] = np.nan
            df['causality_unlikely'] = np.nan
        else:
            df['causality_certain'] = df[causality_col].astype(str).str.contains('Certain', case=False, na=False).astype(int)
            df['causality_probable'] = df[causality_col].astype(str).str.contains('Probable', case=False, na=False).astype(int)
            df['causality_possible'] = df[causality_col].astype(str).str.contains('Possible', case=False, na=False).astype(int)
            df['causality_unlikely'] = df[causality_col].astype(str).str.contains('Unlikely', case=False, na=False).astype(int)
        return df
    
    def _engineer_polypharmacy(self, df: pd.DataFrame) -> pd.DataFrame:
        concom_col = 'concomitant_drugs' if 'concomitant_drugs' in df.columns else ('concomitant_medications' if 'concomitant_medications' in df.columns else None)
        if not concom_col:
            logger.critical("⚠️⚠️⚠️ MISSING COLUMN: 'concomitant_drugs'")
            df['concomitant_count'] = np.nan
            df['polypharmacy_risk'] = np.nan
        else:
            def _parse_concom(val):
                if pd.isna(val) or val is None:
                    return 0
                if isinstance(val, (list, tuple, set)):
                    return len(val)
                s = str(val).strip()
                if not s or s.lower() in ('none', 'nan', 'null', 'unknown', '[]'):
                    return 0
                tokens = [t.strip() for t in s.split(',') if t.strip() and t.strip().lower() not in ('none', 'nan', 'null', 'unknown')]
                return len(tokens)

            df['concomitant_count'] = df[concom_col].apply(_parse_concom)
            df['polypharmacy_risk'] = (df['concomitant_count'] >= 5).astype(int)
        return df
    
    def _engineer_geographic(self, df: pd.DataFrame) -> pd.DataFrame:
        reporter_col = 'reporter_qualification' if 'reporter_qualification' in df.columns else ('reporter_type' if 'reporter_type' in df.columns else None)
        if not reporter_col:
            logger.critical("⚠️⚠️⚠️ MISSING COLUMN: 'reporter'")
            df['reporter_hcp'] = np.nan
            df['reporter_consumer'] = np.nan
        else:
            df['reporter_hcp'] = df[reporter_col].astype(str).str.contains('Physician|Pharmacist|Health', case=False, na=False).astype(int)
            df['reporter_consumer'] = df[reporter_col].astype(str).str.contains('Consumer', case=False, na=False).astype(int)
        
        if 'country' not in df.columns:
            logger.critical("⚠️⚠️⚠️ MISSING COLUMN: 'country'")
            df['country_developed'] = np.nan
        else:
            developed = ['USA', 'GBR', 'DEU', 'FRA', 'CAN', 'AUS', 'JPN']
            df['country_developed'] = df['country'].isin(developed).astype(int)
        return df
