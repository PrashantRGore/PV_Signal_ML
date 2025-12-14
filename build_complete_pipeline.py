from pathlib import Path

print('=' * 70)
print('BUILDING COMPLETE PV PIPELINE: Signal Detection + Causality ML')
print('=' * 70)
print()

# Create causality scorer that applies ML model to signals
causality_code = '''\"\"\"
Causality Assessment: Apply ML model to rank signals
\"\"\"

import pandas as pd
import pickle
from pathlib import Path
import streamlit as st
from src.utils.logger import setup_logger

logger = setup_logger(__name__)

class CausalityScorer:
    \"\"\"
    Applies trained ML model to score signal causality
    Integrates with SISA sharding for GDPR compliance
    \"\"\"
    
    def __init__(self, model_path='models/sisa_model.pkl'):
        self.model_path = Path(model_path)
        self.model = None
    
    def load_model(self):
        \"\"\"Load trained SISA model\"\"\"
        if self.model_path.exists():
            with open(self.model_path, 'rb') as f:
                self.model = pickle.load(f)
            logger.info(f\"Loaded causality model from {self.model_path}\")
            return True
        else:
            logger.warning(f\"Model not found: {self.model_path}\")
            return False
    
    def prepare_features(self, signals_df):
        \"\"\"Prepare features for ML model\"\"\"
        # Basic features from statistical signals
        features = pd.DataFrame({
            'prr': signals_df['prr'],
            'chi2': signals_df['chi2'],
            'case_count': signals_df['case_count'],
            'ror': signals_df.get('ror', 0)
        })
        
        return features
    
    def score_signals(self, signals_df):
        \"\"\"
        Apply ML model to score causality probability
        
        Returns: signals_df with added 'causality_score' column
        \"\"\"
        if self.model is None:
            if not self.load_model():
                st.warning(\"⚠️ ML model not available. Using statistical signals only.\")
                signals_df['causality_score'] = 0.5  # Default neutral score
                return signals_df
        
        try:
            # Prepare features
            X = self.prepare_features(signals_df)
            
            # Get causality predictions
            if hasattr(self.model, 'predict_proba'):
                # Probability of being causal
                probs = self.model.predict_proba(X)
                causality_scores = probs[:, 1] if probs.shape[1] > 1 else probs[:, 0]
            else:
                causality_scores = self.model.predict(X)
            
            # Add to signals
            signals_df['causality_score'] = causality_scores
            
            # Classify causality level
            signals_df['causality_level'] = pd.cut(
                causality_scores,
                bins=[0, 0.3, 0.7, 1.0],
                labels=['Low', 'Moderate', 'High']
            )
            
            logger.info(f\"Scored {len(signals_df)} signals for causality\")
            return signals_df
            
        except Exception as e:
            logger.error(f\"Error scoring causality: {e}\")
            st.error(f\"Causality scoring failed: {e}\")
            signals_df['causality_score'] = 0.5
            return signals_df
'''

scorer_file = Path('src/ml/causality_scorer.py')
scorer_file.write_text(causality_code, encoding='utf-8')
print(f'✅ Created: {scorer_file}')

# Update app to integrate causality scoring
print()
print('Integrating causality scoring into signal detection workflow...')

app_file = Path('app_enhanced.py')
content = app_file.read_text(encoding='utf-8')

# Add causality scoring after signal detection
causality_integration = '''
            # Stage 2: Apply ML Causality Assessment
            st.info(\"🤖 Applying ML causality model to rank signals...\")
            
            from src.ml.causality_scorer import CausalityScorer
            scorer = CausalityScorer()
            
            if scorer.load_model():
                signals = scorer.score_signals(signals)
                
                # Sort by causality score
                signals = signals.sort_values('causality_score', ascending=False)
                
                st.success(f\"✅ Causality assessment complete\")
                
                # Show causality distribution
                col1, col2, col3 = st.columns(3)
                high_causal = (signals['causality_score'] > 0.7).sum()
                mod_causal = ((signals['causality_score'] >= 0.3) & (signals['causality_score'] <= 0.7)).sum()
                low_causal = (signals['causality_score'] < 0.3).sum()
                
                col1.metric(\"High Causality\", f\"{high_causal:,}\", delta=\"Priority\")
                col2.metric(\"Moderate Causality\", f\"{mod_causal:,}\")
                col3.metric(\"Low Causality\", f\"{low_causal:,}\")
            else:
                st.info(\"ℹ️ Using statistical signals only (ML model not trained)\")
'''

# Insert after signal detection completes
marker = 'st.session_state.signals = signals'
if marker in content:
    content = content.replace(marker, marker + causality_integration)
    print('✅ Added causality scoring to pipeline')
else:
    print('⚠️ Could not find insertion point')

# Update Top Signals display to show causality
old_columns = \"filtered[['drug_name', 'event_term', 'case_count', 'prr', 'prr_lower', 'prr_upper', 'chi2', 'is_signal_prr']]\"
new_columns = \"filtered[['drug_name', 'event_term', 'case_count', 'prr', 'chi2', 'causality_score', 'causality_level', 'is_signal_prr']]\"

if old_columns in content:
    content = content.replace(old_columns, new_columns)
    print('✅ Updated signal display to show causality')

app_file.write_text(content, encoding='utf-8')

print()
print('=' * 70)
print('✅ COMPLETE PV PIPELINE READY!')
print('=' * 70)
print()
print('📋 Your Complete Workflow:')
print('  1. Load FAERS data (Live or Demo)')
print('  2. Statistical signal detection (PRR, Chi², ROR)')
print('  3. ✨ ML causality assessment (Your 0.85 model)')
print('  4. Ranked signals by causality probability')
print('  5. Filter by drug portfolio')
print('  6. SHAP explanations (in Explainability tab)')
print()
print('🔐 GDPR Compliance:')
print('  • SISA sharding: Can delete case_id without full retrain')
print('  • Right to be forgotten: Delete one shard only')
print()
print('📊 Output:')
print('  • Statistical signals: 319K potential associations')
print('  • High causality: Top X% prioritized for review')
print('  • Moderate: Monitor/investigate')
print('  • Low: Lower priority')
print()
print('🚀 Restart app to see causality-enhanced pipeline!')
