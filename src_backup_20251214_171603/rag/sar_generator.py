"""
Enhanced SAR Generator with robust model checking
"""
import requests
import json
from datetime import datetime
import os
from pathlib import Path
import time

class SARGenerator:
    def __init__(self, model_name='llama3.2'):
        self.model_name = model_name
        self.ollama_url = 'http://localhost:11434'
        self.output_dir = Path('outputs/sar_reports')
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def check_service(self):
        '''Check if Ollama service is available'''
        try:
            response = requests.get(f'{self.ollama_url}/api/tags', timeout=3)
            return response.status_code == 200, 'Service running'
        except requests.exceptions.ConnectionError:
            return False, 'Ollama not running. Start with: ollama serve'
        except Exception as e:
            return False, f'Error: {str(e)}'
    
    def check_model(self):
        '''Check if model is downloaded and ready'''
        try:
            response = requests.get(f'{self.ollama_url}/api/tags', timeout=3)
            if response.status_code == 200:
                data = response.json()
                models = data.get('models', [])
                
                # Check for exact or partial match
                for model in models:
                    model_name = model.get('name', '')
                    # Handle versioned names like llama3.2:latest
                    if self.model_name in model_name or model_name.startswith(self.model_name):
                        return True, f'Model {model_name} available'
                
                # Model not found
                available = [m['name'] for m in models]
                return False, f'Model not found. Available: {available[:5]}'
        except Exception as e:
            return False, f'Check failed: {str(e)}'
        return False, 'Unknown error'
    
    def pull_model(self):
        '''Pull model using Ollama API'''
        try:
            response = requests.post(
                f'{self.ollama_url}/api/pull',
                json={'name': self.model_name},
                stream=True,
                timeout=600
            )
            
            if response.status_code == 200:
                # Stream pull progress
                for line in response.iter_lines():
                    if line:
                        data = json.loads(line)
                        status = data.get('status', '')
                        if 'success' in status.lower():
                            return True
                return True
            else:
                return False
        except Exception as e:
            return False
    
    def generate_sar(self, signal_data, model_name=None):
        '''Generate SAR with comprehensive checks'''
        if model_name:
            self.model_name = model_name
        
        # Health check
        service_ok, service_msg = self.check_service()
        if not service_ok:
            raise Exception(service_msg)
        
        model_ok, model_msg = self.check_model()
        if not model_ok:
            raise Exception(f'{model_msg}. Pull with: ollama pull {self.model_name}')
        
        # Generate prompt
        prompt = self._create_sar_prompt(signal_data)
        
        # Call Ollama API
        try:
            response = requests.post(
                f'{self.ollama_url}/api/generate',
                json={
                    'model': self.model_name,
                    'prompt': prompt,
                    'stream': False,
                    'options': {
                        'temperature': 0.7,
                        'num_predict': 1000
                    }
                },
                timeout=120
            )
            
            if response.status_code == 200:
                result = response.json()
                report_text = result.get('response', '')
                
                return {
                    'report': report_text,
                    'signal_data': signal_data,
                    'model_used': self.model_name,
                    'generated_at': datetime.now().isoformat(),
                    'evidence_sources': ['FAERS Database', 'Statistical Analysis', 'Disproportionality Metrics']
                }
            else:
                raise Exception(f'API error {response.status_code}: {response.text}')
        
        except requests.exceptions.Timeout:
            raise Exception('Request timeout. Model loading may take time, try again.')
        except Exception as e:
            raise Exception(f'Generation failed: {str(e)}')
    
    def _create_sar_prompt(self, signal_data):
        '''Create structured SAR prompt'''
        return f'''You are a pharmacovigilance expert. Generate a Signal Assessment Report (SAR) for:

Drug: {signal_data['drug_name']}
Adverse Event: {signal_data['event_name']}
Case Count: {signal_data['count']}
PRR: {signal_data['prr']:.2f}
Chi-Square: {signal_data['chi_square']:.2f}

Provide a structured report with:
1. Executive Summary
2. Signal Description
3. Statistical Evidence
4. Clinical Significance
5. Regulatory Context
6. Recommendations

Keep it professional and concise (800-1000 words).'''
    
    def save_report(self, sar_result):
        '''Save SAR to file'''
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        drug_name = sar_result['signal_data']['drug_name'].replace(' ', '_')[:30]
        filename = f'SAR_{drug_name}_{timestamp}.txt'
        filepath = self.output_dir / filename
        
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write('='*80 + '\n')
            f.write('SIGNAL ASSESSMENT REPORT (SAR)\n')
            f.write('='*80 + '\n\n')
            f.write(f"Generated: {sar_result['generated_at']}\n")
            f.write(f"Model: {sar_result['model_used']}\n")
            f.write(f"Model Version: {sar_result.get('model_version', 1)}\n\n")
            f.write('-'*80 + '\n')
            f.write(sar_result['report'])
            f.write('\n\n' + '-'*80 + '\n')
            f.write('Evidence Sources:\n')
            for src in sar_result['evidence_sources']:
                f.write(f'- {src}\n')
        
        return str(filepath)
