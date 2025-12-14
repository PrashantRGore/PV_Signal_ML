"""
Enhanced SAR Generator with Causality Assessment and Literature Evidence
Production-ready version
"""
import requests
import json
from datetime import datetime
import os
from pathlib import Path

# Import new modules
try:
    from src.regulatory.causality_assessor import CausalityAssessor
    CAUSALITY_AVAILABLE = True
except:
    CAUSALITY_AVAILABLE = False

try:
    from src.literature.pubmed_miner import LiteratureEvidenceGenerator
    PUBMED_AVAILABLE = True
except:
    PUBMED_AVAILABLE = False

class SARGenerator:
    def __init__(self, model_name='llama3.2'):
        self.model_name = model_name
        self.ollama_url = 'http://localhost:11434'
        self.output_dir = Path('outputs/sar_reports')
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize enhanced features if available
        self.causality_assessor = CausalityAssessor() if CAUSALITY_AVAILABLE else None
        self.lit_generator = LiteratureEvidenceGenerator() if PUBMED_AVAILABLE else None
    
    def check_service(self):
        try:
            response = requests.get(f'{self.ollama_url}/api/tags', timeout=3)
            return response.status_code == 200, 'Service running'
        except:
            return False, 'Ollama not running'
    
    def check_model(self):
        try:
            response = requests.get(f'{self.ollama_url}/api/tags', timeout=3)
            if response.status_code == 200:
                models = response.json().get('models', [])
                for model in models:
                    if self.model_name in model.get('name', ''):
                        return True, f'Model available'
                return False, 'Model not found'
        except:
            return False, 'Check failed'
    
    def generate_sar(self, signal_data, model_choice):
        '''Generate comprehensive SAR with causality and literature'''
        self.model_name = model_choice
        
        # Service check
        service_ok, _ = self.check_service()
        if not service_ok:
            raise Exception('Ollama service not running')
        
        model_ok, _ = self.check_model()
        if not model_ok:
            raise Exception(f'Model {self.model_name} not available')
        
        # Get causality assessment
        causality = None
        causality_text = ''
        if self.causality_assessor:
            try:
                causality = self.causality_assessor.assess_signal(signal_data, method='both')
                causality_text = self._format_causality(causality)
            except Exception as e:
                causality_text = f'Causality assessment unavailable: {str(e)}'
        
        # Get literature evidence
        literature = None
        literature_text = ''
        if self.lit_generator:
            try:
                literature = self.lit_generator.generate_evidence_section(
                    signal_data['drug_name'],
                    signal_data['event_name'],
                    max_papers=5
                )
                literature_text = literature['evidence_text']
            except Exception as e:
                literature_text = f'Literature search unavailable: {str(e)}'
        
        # Create enhanced prompt
        prompt = self._create_enhanced_prompt(signal_data, causality_text, literature_text)
        
        # Generate with LLM
        try:
            response = requests.post(
                f'{self.ollama_url}/api/generate',
                json={
                    'model': self.model_name,
                    'prompt': prompt,
                    'stream': False,
                    'options': {
                        'temperature': 0.7,
                        'num_predict': 1500
                    }
                },
                timeout=180
            )
            
            if response.status_code == 200:
                result = response.json()
                report_text = result.get('response', '')
                
                return {
                    'report': report_text,
                    'signal_data': signal_data,
                    'causality_assessment': causality,
                    'literature_evidence': literature,
                    'model_used': self.model_name,
                    'generated_at': datetime.now().isoformat(),
                    'evidence_sources': self._get_evidence_sources(causality, literature)
                }
            else:
                raise Exception(f'API error {response.status_code}')
        
        except Exception as e:
            raise Exception(f'Generation failed: {str(e)}')
    
    def _format_causality(self, causality):
        '''Format causality assessment for prompt'''
        lines = ['**CAUSALITY ASSESSMENT:**', '']
        
        if 'who_umc' in causality:
            who = causality['who_umc']
            lines.append(f'WHO-UMC Classification: {who["category"]}')
            lines.append(f'Rationale: {who["rationale"]}')
            lines.append('')
        
        if 'naranjo' in causality:
            naranjo = causality['naranjo']
            lines.append(f'Naranjo Algorithm Score: {naranjo["score"]} ({naranjo["category"]})')
            lines.append(f'Interpretation: {naranjo["interpretation"]}')
            lines.append('')
        
        if 'consensus' in causality:
            consensus = causality['consensus']
            lines.append(f'Consensus Assessment: {consensus["level"]}')
            lines.append('')
        
        return '\n'.join(lines)
    
    def _create_enhanced_prompt(self, signal_data, causality_text, literature_text):
        '''Create comprehensive prompt'''
        prompt = f'''Generate a comprehensive Signal Assessment Report (SAR) for regulatory submission.

**SIGNAL INFORMATION:**
- Drug: {signal_data['drug_name']}
- Adverse Event: {signal_data['event_name']}
- Case Count: {signal_data['count']}
- PRR: {signal_data['prr']:.2f}
- Chi-Square: {signal_data['chi_square']:.2f}

{causality_text}

{literature_text}

**REQUIRED REPORT STRUCTURE:**

1. **Executive Summary** (150 words)
   - Signal overview with causality assessment
   - Key findings and recommendations

2. **Signal Description** (200 words)
   - Drug pharmacology and mechanism
   - Adverse event clinical presentation
   - Biological plausibility

3. **Statistical Evidence** (150 words)
   - Disproportionality analysis interpretation
   - PRR and Chi-square significance
   - Signal strength evaluation

4. **Causality Assessment** (200 words)
   - WHO-UMC classification with rationale
   - Naranjo algorithm score interpretation
   - Overall causality conclusion

5. **Literature Review** (200 words)
   - Summary of published evidence
   - Consistency with current signal
   - Evidence quality assessment

6. **Clinical Significance** (150 words)
   - Patient risk implications
   - Clinical management considerations
   - Risk-benefit evaluation

7. **Regulatory Context** (150 words)
   - FDA/EMA/WHO guidelines applicability
   - Required regulatory actions
   - Reporting obligations

8. **Recommendations** (200 words)
   - Immediate actions required
   - Further investigations needed
   - Label update considerations
   - Risk minimization measures

9. **Limitations** (100 words)
   - Data source constraints
   - Methodological limitations
   - Confounding factors

10. **Conclusion** (100 words)
    - Causality summary
    - Next steps

**COMPLIANCE REQUIREMENTS:**
- Follow ICH E2A guidelines
- Include all causality evidence
- Reference literature sources
- Maintain regulatory tone
- Provide actionable recommendations

Generate the complete SAR report now:'''
        
        return prompt
    
    def _get_evidence_sources(self, causality, literature):
        '''Compile evidence sources'''
        sources = ['FAERS Database', 'Statistical Disproportionality Analysis']
        
        if causality:
            sources.extend(['WHO-UMC Causality Assessment', 'Naranjo Algorithm'])
        
        if literature and literature.get('n_papers', 0) > 0:
            sources.append(f'PubMed Literature ({literature["n_papers"]} papers)')
        
        return sources
    
    def save_report(self, sar_result):
        '''Save comprehensive SAR report'''
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        drug_name = sar_result['signal_data']['drug_name'].replace(' ', '_')[:30]
        filename = f'SAR_{drug_name}_{timestamp}.txt'
        filepath = self.output_dir / filename
        
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write('='*100 + '\n')
            f.write('SIGNAL ASSESSMENT REPORT (SAR)\n')
            f.write('Regulatory Standard: ICH E2A\n')
            f.write('='*100 + '\n\n')
            
            f.write(f'Generated: {sar_result["generated_at"]}\n')
            f.write(f'LLM Model: {sar_result["model_used"]}\n')
            f.write(f'Model Version: {sar_result.get("model_version", 1)}\n\n')
            
            f.write('-'*100 + '\n')
            f.write('SIGNAL DATA\n')
            f.write('-'*100 + '\n')
            for key, value in sar_result['signal_data'].items():
                f.write(f'{key}: {value}\n')
            
            # Causality section
            if sar_result.get('causality_assessment'):
                f.write('\n' + '-'*100 + '\n')
                f.write('CAUSALITY ASSESSMENT\n')
                f.write('-'*100 + '\n\n')
                
                causality = sar_result['causality_assessment']
                if 'who_umc' in causality:
                    f.write(f'WHO-UMC: {causality["who_umc"]["category"]}\n')
                    f.write(f'Rationale: {causality["who_umc"]["rationale"]}\n\n')
                
                if 'naranjo' in causality:
                    f.write(f'Naranjo Score: {causality["naranjo"]["score"]} ({causality["naranjo"]["category"]})\n\n')
                
                if 'consensus' in causality:
                    f.write(f'Consensus: {causality["consensus"]["level"]}\n\n')
            
            # Literature section
            if sar_result.get('literature_evidence'):
                lit = sar_result['literature_evidence']
                if lit.get('n_papers', 0) > 0:
                    f.write('-'*100 + '\n')
                    f.write('LITERATURE EVIDENCE\n')
                    f.write('-'*100 + '\n\n')
                    f.write(lit['evidence_text'])
                    f.write('\n\n')
            
            # Generated report
            f.write('-'*100 + '\n')
            f.write('GENERATED REPORT\n')
            f.write('-'*100 + '\n\n')
            f.write(sar_result['report'])
            
            # Evidence sources
            f.write('\n\n' + '-'*100 + '\n')
            f.write('EVIDENCE SOURCES\n')
            f.write('-'*100 + '\n')
            for i, source in enumerate(sar_result['evidence_sources'], 1):
                f.write(f'{i}. {source}\n')
            
            # References
            if sar_result.get('literature_evidence') and sar_result['literature_evidence'].get('references'):
                f.write('\n' + '-'*100 + '\n')
                f.write('REFERENCES\n')
                f.write('-'*100 + '\n')
                for i, ref in enumerate(sar_result['literature_evidence']['references'], 1):
                    f.write(f'[{i}] {ref["citation"]}\n')
                    f.write(f'    PMID: {ref["pmid"]} | {ref["url"]}\n\n')
        
        return str(filepath)
