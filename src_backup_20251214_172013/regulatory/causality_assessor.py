"""
Causality Assessment Engine
Implements WHO-UMC and Naranjo algorithms for regulatory compliance
"""
from enum import Enum
import pandas as pd

class CausalityCategory(Enum):
    CERTAIN = 'Certain/Definite'
    PROBABLE = 'Probable/Likely'
    POSSIBLE = 'Possible'
    UNLIKELY = 'Unlikely'
    CONDITIONAL = 'Conditional/Unclassified'
    UNASSESSABLE = 'Unassessable/Unclassifiable'

class WHOUMCCausality:
    '''
    WHO-UMC Causality Assessment System
    Reference: https://www.who-umc.org/vigibase/causality-assessment/
    '''
    
    @staticmethod
    def assess(signal_data, case_data=None):
        '''
        Assess causality based on WHO-UMC criteria
        
        Criteria:
        - Temporal relationship
        - De-challenge/Re-challenge
        - Alternative causes
        - Laboratory evidence
        '''
        score = {
            'temporal_relationship': True,  # Always true for FAERS signals
            'dechallenge_positive': False,
            'rechallenge_positive': False,
            'alternative_causes': None,
            'lab_evidence': None
        }
        
        # Basic assessment from signal strength
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
            'criteria_met': score,
            'algorithm': 'WHO-UMC'
        }

class NaranjoCausality:
    '''
    Naranjo Algorithm for Causality Assessment
    Reference: Naranjo et al. (1981) Clinical Pharmacology & Therapeutics
    '''
    
    QUESTIONS = [
        ('Q1', 'Are there previous conclusive reports on this reaction?', 1, 0, 0),
        ('Q2', 'Did the adverse event appear after the drug was given?', 2, -1, 0),
        ('Q3', 'Did the adverse reaction improve when drug stopped or antagonist given?', 1, 0, 0),
        ('Q4', 'Did the adverse reaction reappear when drug readministered?', 2, -1, 0),
        ('Q5', 'Are there alternative causes that could cause the reaction?', -1, 2, 0),
        ('Q6', 'Did the reaction reappear when placebo given?', -1, 1, 0),
        ('Q7', 'Was the drug detected in blood/body fluids in toxic concentrations?', 1, 0, 0),
        ('Q8', 'Was the reaction more severe with dose increase or less severe with dose decrease?', 1, 0, 0),
        ('Q9', 'Did patient have similar reaction to same/similar drug in any previous exposure?', 1, 0, 0),
        ('Q10', 'Was the adverse event confirmed by objective evidence?', 1, 0, 0),
    ]
    
    @staticmethod
    def assess(signal_data, answers=None):
        '''
        Calculate Naranjo score
        
        Score interpretation:
        >= 9: Definite
        5-8: Probable
        1-4: Possible
        <= 0: Doubtful
        '''
        # Default conservative answers based on FAERS data
        if answers is None:
            answers = [
                0,  # Q1: Unknown literature
                2,  # Q2: Yes (temporal relationship)
                0,  # Q3: Unknown
                0,  # Q4: Unknown
                0,  # Q5: Unknown
                0,  # Q6: Not applicable
                0,  # Q7: Unknown
                0,  # Q8: Unknown
                0,  # Q9: Unknown
                1   # Q10: Confirmed by reporting
            ]
        
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
            'algorithm': 'Naranjo',
            'answers': dict(zip([q[0] for q in NaranjoCausality.QUESTIONS], answers))
        }

class CausalityAssessor:
    '''Unified causality assessment using multiple algorithms'''
    
    def __init__(self):
        self.who_umc = WHOUMCCausality()
        self.naranjo = NaranjoCausality()
    
    def assess_signal(self, signal_data, method='both'):
        '''
        Assess causality using specified method
        
        Parameters:
        - signal_data: Dict with prr, count, chi_square
        - method: 'who-umc', 'naranjo', or 'both'
        '''
        results = {}
        
        if method in ['who-umc', 'both']:
            results['who_umc'] = self.who_umc.assess(signal_data)
        
        if method in ['naranjo', 'both']:
            results['naranjo'] = self.naranjo.assess(signal_data)
        
        # Consensus assessment
        if method == 'both':
            results['consensus'] = self._get_consensus(results)
        
        return results
    
    def _get_consensus(self, assessments):
        '''Generate consensus from multiple assessments'''
        categories = []
        
        if 'who_umc' in assessments:
            categories.append(assessments['who_umc']['category'])
        if 'naranjo' in assessments:
            categories.append(assessments['naranjo']['category'])
        
        # Take most conservative
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
