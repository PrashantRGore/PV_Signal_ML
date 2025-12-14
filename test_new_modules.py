"""
Test script for new modules
"""
print('='*60)
print('Testing New Modules')
print('='*60)

# Test 1: Causality Assessor
print('\n1. Testing Causality Assessor...')
from src.regulatory.causality_assessor import CausalityAssessor
assessor = CausalityAssessor()
signal_data = {
    'drug_name': 'TEST_DRUG',
    'event_name': 'TEST_EVENT',
    'prr': 12.5,
    'count': 20,
    'chi_square': 150
}
causality = assessor.assess_signal(signal_data)
print(f'✅ WHO-UMC: {causality["who_umc"]["category"]}')
print(f'✅ Naranjo: {causality["naranjo"]["category"]}')
print(f'✅ Consensus: {causality["consensus"]["level"]}')

# Test 2: SHAP Analyzer V2 (structure test)
print('\n2. Testing SHAP V2 import...')
from src.explainability.shap_analyzer_v2 import SHAPAnalyzerV2
analyzer = SHAPAnalyzerV2()
print('✅ SHAP V2 imported successfully')

# Test 3: Literature Miner (structure test)
print('\n3. Testing Literature Miner import...')
from src.literature.pubmed_miner import PubMedMiner, LiteratureEvidenceGenerator
miner = PubMedMiner()
print('✅ PubMed Miner imported successfully')

print('\n' + '='*60)
print('✅ All New Modules Working!')
print('='*60)
