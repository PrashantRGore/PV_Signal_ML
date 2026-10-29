# src/ml/enhanced_feature_pipeline.py

import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler

class EnhancedFeaturePipeline:
    """
    Enhanced feature pipeline utilizing all 35 columns from synthetic FAERS dataset
    """
    
    def __init__(self):
        self.feature_names = []
        self.scaler = StandardScaler()
        
    def transform(self, df):
        """Transform signal data into ML features using all available columns"""
        
        features = pd.DataFrame(index=df.index)
        
        # ===== STATISTICAL FEATURES (Original) =====
        features['case_count'] = df.get('case_count', 0)
        features['prr'] = df.get('prr', 1.0)
        features['chi2'] = df.get('chi2', 0.0)
        features['ror'] = df.get('ror', 1.0)
        
        # ===== TEMPORAL FEATURES (NEW) =====
        if 'ttp_days' in df.columns:
            features['ttp_days'] = df['ttp_days'].fillna(df['ttp_days'].median())
            features['ttp_log'] = np.log1p(features['ttp_days'])
            
            # Time to onset categories (critical for causality)
            features['ttp_acute'] = (features['ttp_days'] <= 7).astype(int)
            features['ttp_subacute'] = ((features['ttp_days'] > 7) & (features['ttp_days'] <= 30)).astype(int)
            features['ttp_chronic'] = (features['ttp_days'] > 30).astype(int)
        
        if 'therapy_duration' in df.columns:
            features['therapy_duration'] = df['therapy_duration'].fillna(0)
            if 'ttp_days' in features:
                features['therapy_ttp_ratio'] = features['therapy_duration'] / (features['ttp_days'] + 1)
        
        if 'event_duration' in df.columns:
            features['event_duration'] = df['event_duration'].fillna(0)
        
        # ===== CAUSALITY FEATURES (NEW - CRITICAL) =====
        if 'dechallenge' in df.columns:
            features['dechallenge_positive'] = (df['dechallenge'] == 'Positive').astype(int)
            features['dechallenge_negative'] = (df['dechallenge'] == 'Negative').astype(int)
        
        if 'rechallenge' in df.columns:
            features['rechallenge_positive'] = (df['rechallenge'] == 'Positive').astype(int)
            features['rechallenge_negative'] = (df['rechallenge'] == 'Negative').astype(int)
        
        if 'action_drug' in df.columns:
            features['drug_withdrawn'] = (df['action_drug'] == 'Drug withdrawn').astype(int)
            features['dose_reduced'] = (df['action_drug'] == 'Dose reduced').astype(int)
        
        if 'causality' in df.columns:
            causality_score_map = {
                'Certain': 1.0,
                'Probable/Likely': 0.8,
                'Possible': 0.5,
                'Unlikely': 0.2,
                'Conditional/Unclassified': 0.1
            }
            features['causality_score'] = df['causality'].map(causality_score_map).fillna(0.5)
        
        # ===== CONCOMITANT MEDICATION FEATURES (NEW) =====
        if 'concomitant_drugs' in df.columns:
            features['has_concomitant'] = (df['concomitant_drugs'] != 'None').astype(int)
            features['num_concomitant'] = df['concomitant_drugs'].apply(
                lambda x: len(str(x).split(';')) if pd.notna(x) and x != 'None' else 0
            )
        
        # ===== MEDICAL HISTORY FEATURES (NEW) =====
        if 'medical_history' in df.columns:
            features['has_medical_history'] = (df['medical_history'] != 'None').astype(int)
            features['num_conditions'] = df['medical_history'].apply(
                lambda x: len(str(x).split(';')) if pd.notna(x) and x != 'None' else 0
            )
        
        # ===== SERIOUSNESS FEATURES =====
        if 'serious' in df.columns:
            features['is_death'] = (df['serious'] == 'Death').astype(int)
            features['is_life_threatening'] = (df['serious'] == 'Life-threatening').astype(int)
            features['is_hospitalization'] = (df['serious'].str.contains('Hospitalization', na=False)).astype(int)
            features['is_serious'] = (~df['serious'].isin(['Non-serious', 'Not serious'])).astype(int)
        
        # ===== OUTCOME FEATURES =====
        if 'outcome' in df.columns:
            features['outcome_fatal'] = (df['outcome'] == 'Fatal').astype(int)
            features['outcome_recovered'] = (df['outcome'] == 'Recovered/Resolved').astype(int)
            features['outcome_recovering'] = (df['outcome'] == 'Recovering/Resolving').astype(int)
        
        # ===== DEMOGRAPHIC FEATURES =====
        if 'age' in df.columns:
            features['age'] = df['age'].fillna(df['age'].median())
            features['age_elderly'] = (features['age'] >= 65).astype(int)
            features['age_pediatric'] = (features['age'] < 18).astype(int)
        
        if 'sex' in df.columns:
            features['sex_male'] = (df['sex'] == 'Male').astype(int)
            features['sex_female'] = (df['sex'] == 'Female').astype(int)
        
        if 'weight' in df.columns:
            features['weight'] = df['weight'].fillna(df['weight'].median())
            if 'age' in features:
                features['bmi_proxy'] = features['weight'] / (features['age'] + 1)
        
        
        # ===== USE AGGREGATED FEATURES FROM SIGNALS_DF =====
# Add after line 106 (after demographics section):

        # Use aggregated features if base columns are missing
        if 'age_mean' in df.columns and 'age' not in features:
            features['age'] = df['age_mean']
        
        if 'weight_mean' in df.columns and 'weight' not in features:
            features['weight'] = df['weight_mean']
        
        # Use percentage features for demographics
        if 'elderly_pct' in df.columns:
            features['age_elderly'] = df['elderly_pct']
        
        if 'pediatric_pct' in df.columns:
            features['age_pediatric'] = df['pediatric_pct']
        
        if 'male_pct' in df.columns:
            features['sex_male'] = df['male_pct']
        
        if 'female_pct' in df.columns:
            features['sex_female'] = df['female_pct']
        
        # Use temporal aggregates
        if 'ttp_mean' in df.columns and 'ttp_days' not in features:
            features['ttp_days'] = df['ttp_mean']
        
        # Use causality aggregates
        if 'dechallenge_positive_pct' in df.columns:
            features['dechallenge_positive'] = df['dechallenge_positive_pct']
        
        if 'dechallenge_negative_pct' in df.columns:
            features['dechallenge_negative'] = df['dechallenge_negative_pct']
        
        if 'rechallenge_positive_pct' in df.columns:
            features['rechallenge_positive'] = df['rechallenge_positive_pct']
        
        if 'rechallenge_negative_pct' in df.columns:
            features['rechallenge_negative'] = df['rechallenge_negative_pct']
        
        if 'causality_score_mean' in df.columns:
            features['causality_score'] = df['causality_score_mean']
        
        # Use seriousness aggregates
        if 'death_pct' in df.columns:
            features['is_death'] = df['death_pct']
        
        if 'life_threatening_pct' in df.columns:
            features['is_life_threatening'] = df['life_threatening_pct']
        
        if 'hospitalization_pct' in df.columns:
            features['is_hospitalization'] = df['hospitalization_pct']
        
        if 'fatal_pct' in df.columns:
            features['outcome_fatal'] = df['fatal_pct']
        
        if 'recovered_pct' in df.columns:
            features['outcome_recovered'] = df['recovered_pct']
        
        if 'serious_pct' in df.columns:
            features['is_serious'] = df['serious_pct']
        
        # Use concomitant aggregates
        if 'has_concomitant_pct' in df.columns:
            features['has_concomitant'] = df['has_concomitant_pct']
        
        if 'has_medical_history_pct' in df.columns:
            features['has_medical_history'] = df['has_medical_history_pct']
        
        # Use reporter aggregates
        if 'hcp_reporter_pct' in df.columns:
            features['reporter_hcp'] = df['hcp_reporter_pct']
        
        if 'spontaneous_pct' in df.columns:
            features['report_spontaneous'] = df['spontaneous_pct']
        
        # Use drug action aggregates
        if 'drug_withdrawn_pct' in df.columns:
            features['drug_withdrawn'] = df['drug_withdrawn_pct']
        
        if 'dose_reduced_pct' in df.columns:
            features['dose_reduced'] = df['dose_reduced_pct']



# ===== LAB VALUE FEATURES (NEW - for hepatic/renal signals) =====
        if 'alt' in df.columns:
            features['alt'] = df['alt'].fillna(df['alt'].median())
            features['alt_elevated'] = (features['alt'] > 55).astype(int)
        
        if 'ast' in df.columns:
            features['ast'] = df['ast'].fillna(df['ast'].median())
            features['ast_elevated'] = (features['ast'] > 48).astype(int)
        
        if 'bilirubin' in df.columns:
            features['bilirubin'] = df['bilirubin'].fillna(df['bilirubin'].median())
            features['bilirubin_elevated'] = (features['bilirubin'] > 1.2).astype(int)
        
        if 'creatinine' in df.columns:
            features['creatinine'] = df['creatinine'].fillna(df['creatinine'].median())
            features['creatinine_elevated'] = (features['creatinine'] > 1.3).astype(int)
        
        if 'bun' in df.columns:
            features['bun'] = df['bun'].fillna(df['bun'].median())
            features['bun_elevated'] = (features['bun'] > 20).astype(int)
        
        # Composite scores
        if all(col in features for col in ['alt_elevated', 'ast_elevated', 'bilirubin_elevated']):
            features['hepatic_injury_score'] = (
                features['alt_elevated'] + 
                features['ast_elevated'] + 
                features['bilirubin_elevated']
            )
        
        if all(col in features for col in ['creatinine_elevated', 'bun_elevated']):
            features['renal_injury_score'] = (
                features['creatinine_elevated'] + 
                features['bun_elevated']
            )
        
        # ===== MEDDRA HIERARCHY FEATURES (NEW) =====
        if 'soc' in df.columns:
            soc_dummies = pd.get_dummies(df['soc'], prefix='soc')
            for col in soc_dummies.columns:
                features[col] = soc_dummies[col]
        
        # ===== REPORTER QUALITY FEATURES =====
        if 'reporter_qualification' in df.columns:
            features['reporter_hcp'] = df['reporter_qualification'].isin([
                'Physician', 'Pharmacist', 'Other health professional'
            ]).astype(int)
        
        if 'report_source' in df.columns:
            features['report_spontaneous'] = (df['report_source'] == 'Spontaneous').astype(int)
            features['report_clinical_trial'] = (df['report_source'] == 'Clinical trial').astype(int)
        
        # ===== INTERACTION FEATURES (Advanced) =====
        if 'dechallenge_positive' in features and 'rechallenge_positive' in features:
            features['dechallenge_rechallenge_both_positive'] = (
                features['dechallenge_positive'] * features['rechallenge_positive']
            )
        
        if 'ttp_acute' in features and 'dechallenge_positive' in features:
            features['acute_with_positive_dechallenge'] = (
                features['ttp_acute'] * features['dechallenge_positive']
            )
        
        if 'is_serious' in features and 'hepatic_injury_score' in features:
            features['serious_hepatic'] = features['is_serious'] * features['hepatic_injury_score']
        
        # Store feature names
        self.feature_names = list(features.columns)
        
        # Fill any remaining NaN
        features = features.fillna(0)
        
        return features
    

    def fit_transform(self, df):
        """Fit and transform - for sklearn compatibility"""
        # EnhancedFeaturePipeline doesn't need fitting, just transform
        return self.transform(df)
    
    def fit(self, df):
        """Fit method - for sklearn compatibility (no-op)"""
        return self


    def save(self, filepath):
        """Save the feature pipeline to disk"""
        import pickle
        from pathlib import Path
        
        filepath = Path(filepath)
        filepath.parent.mkdir(parents=True, exist_ok=True)
        
        with open(filepath, 'wb') as f:
            pickle.dump(self, f)
    
    @staticmethod
    def load(filepath):
        """Load the feature pipeline from disk"""
        import pickle
        
        with open(filepath, 'rb') as f:
            return pickle.load(f)

    def get_feature_importance_categories(self):
        """Return feature categories for interpretation"""
        return {
            'statistical': ['case_count', 'prr', 'chi2', 'ror'],
            'temporal': [f for f in self.feature_names if 'ttp_' in f or 'duration' in f],
            'causality': [f for f in self.feature_names if any(x in f for x in ['dechallenge', 'rechallenge', 'causality'])],
            'demographic': [f for f in self.feature_names if any(x in f for x in ['age', 'sex', 'weight'])],
            'lab_values': [f for f in self.feature_names if any(x in f for x in ['alt', 'ast', 'bilirubin', 'creatinine', 'bun', 'hepatic', 'renal'])],
            'concomitant': [f for f in self.feature_names if 'concomitant' in f],
            'medical_history': [f for f in self.feature_names if 'medical_history' in f or 'conditions' in f],
            'seriousness': [f for f in self.feature_names if 'serious' in f or 'outcome' in f or 'death' in f],
            'meddra': [f for f in self.feature_names if 'soc_' in f],
            'interactions': [f for f in self.feature_names if 'both_positive' in f or 'acute_with' in f or 'serious_' in f]
        }



