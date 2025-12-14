"""
Causality Assessment Engine - WHO-UMC + Naranjo
"""
from enum import Enum

class CausalityCategory(Enum):
    CERTAIN = 'Certain/Definite'
    PROBABLE = 'Probable/Likely'
    POSSIBLE = 'Possible'
    UNLIKELY = 'Unlikely'
    CONDITIONAL = 'Conditional/Unclassified'
    UNASSESSABLE = 'Unassessable/Unclassifiable'

class WHOUMCCausality:
    @staticmethod
    def assess(signal_data, case_data=None):
        prr = signal_data.get('prr', 0)
        case_count = signal_data.get('count', 0)
        
        if prr >= 10 and case_count >= 10:
            category = CausalityCategory.PROBABLE
            rationale = 'Strong statistical signal with adequate case count'
        elif prr >= 5 and case_count >= 5:
            category = CausalityCategory.POSSIBLE
            rationale = 'Moderate statistical signal'
        elif prr >= 2 and case_count >= 3:
            category = CausalityCategory.POSSIBLE
            rationale = 'Weak but detectable signal'
        else:
            category = CausalityCategory.UNLIKELY
            rationale = 'Insufficient evidence'
        
        return {
            'category': category.value,
            'rationale': rationale,
            'algorithm': 'WHO-UMC'
        }

class NaranjoCausality:
    @staticmethod
    def assess(signal_data, answers=None):
        if answers is None:
            answers = [0, 2, 0, 0, 0, 0, 0, 0, 0, 1]
        
        total_score = sum(answers)
        
        if total_score >= 9:
            category = 'Definite'
        elif total_score >= 5:
            category = 'Probable'
        elif total_score >= 1:
            category = 'Possible'
        else:
            category = 'Doubtful'
        
        return {
            'score': total_score,
            'category': category,
            'interpretation': f'Naranjo score: {total_score} ({category})',
            'algorithm': 'Naranjo'
        }

class CausalityAssessor:
    def __init__(self):
        self.who_umc = WHOUMCCausality()
        self.naranjo = NaranjoCausality()
    
    def assess_signal(self, signal_data, method='both'):
        results = {}
        
        if method in ['who-umc', 'both']:
            results['who_umc'] = self.who_umc.assess(signal_data)
        
        if method in ['naranjo', 'both']:
            results['naranjo'] = self.naranjo.assess(signal_data)
        
        if method == 'both':
            results['consensus'] = self._get_consensus(results)
        
        return results
    
    def _get_consensus(self, assessments):
        categories = []
        
        if 'who_umc' in assessments:
            categories.append(assessments['who_umc']['category'])
        if 'naranjo' in assessments:
            categories.append(assessments['naranjo']['category'])
        
        if any('Definite' in cat or 'Certain' in cat for cat in categories):
            consensus = 'High Causality'
        elif any('Probable' in cat or 'Likely' in cat for cat in categories):
            consensus = 'Moderate Causality'
        elif any('Possible' in cat for cat in categories):
            consensus = 'Possible Causality'
        else:
            consensus = 'Low Causality'
        
        return {
            'level': consensus,
            'rationale': 'Consensus from WHO-UMC and Naranjo assessments'
        }
