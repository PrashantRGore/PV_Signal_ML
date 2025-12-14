"""
Enhanced Service Health Checker with model pull capability
"""
import requests
import subprocess
import time
import json

class ServiceHealthChecker:
    def __init__(self):
        self.ollama_url = 'http://localhost:11434'
        self.mlflow_url = 'http://localhost:5000'
    
    def check_ollama(self):
        try:
            response = requests.get(f'{self.ollama_url}/api/tags', timeout=2)
            return response.status_code == 200
        except:
            return False
    
    def check_ollama_model(self, model_name):
        '''Check if specific model exists'''
        try:
            response = requests.get(f'{self.ollama_url}/api/tags', timeout=2)
            if response.status_code == 200:
                models = response.json().get('models', [])
                for model in models:
                    name = model.get('name', '')
                    if model_name in name or name.startswith(model_name):
                        return True
        except:
            pass
        return False
    
    def get_available_models(self):
        '''Get list of downloaded models'''
        try:
            response = requests.get(f'{self.ollama_url}/api/tags', timeout=2)
            if response.status_code == 200:
                models = response.json().get('models', [])
                return [m['name'] for m in models]
        except:
            pass
        return []
    
    def start_ollama(self):
        try:
            subprocess.Popen(
                ['powershell', '-Command', 'Start-Process', 'ollama', 'serve'],
                creationflags=subprocess.CREATE_NO_WINDOW
            )
            time.sleep(3)
            return self.check_ollama()
        except:
            return False
    
    def pull_ollama_model(self, model_name):
        '''Pull Ollama model using subprocess'''
        try:
            result = subprocess.run(
                ['ollama', 'pull', model_name],
                capture_output=True,
                text=True,
                timeout=600
            )
            return result.returncode == 0
        except:
            return False
    
    def check_mlflow(self):
        try:
            response = requests.get(self.mlflow_url, timeout=2)
            return response.status_code == 200
        except:
            return False
    
    def start_mlflow_ui(self, port=5000):
        try:
            subprocess.Popen(
                ['mlflow', 'ui', '--port', str(port)],
                creationflags=subprocess.CREATE_NO_WINDOW,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
            )
            time.sleep(2)
            return True
        except:
            return False
    
    def get_service_status(self):
        return {
            'ollama': self.check_ollama(),
            'mlflow': self.check_mlflow()
        }
