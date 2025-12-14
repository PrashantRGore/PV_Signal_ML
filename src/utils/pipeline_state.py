"""
Pipeline State Manager
Manages workflow state and dependencies between tabs
Ensures proper execution order: Signals → SISA → SHAP → SAR
"""

class PipelineState:
    def __init__(self):
        self.stages = {
            'data_loaded': False,
            'signals_detected': False,
            'sisa_trained': False,
            'shap_computed': False,
            'mlflow_run_active': False
        }
        
        self.data = {
            'raw_data': None,
            'signals': None,
            'sisa_model': None,
            'sisa_results': None,
            'shap_results': None,
            'model_version': 1
        }
        
        self.mlflow_run_id = None
    
    def can_proceed_to(self, stage):
        \"\"\"Check if prerequisites are met for a stage\"\"\"
        dependencies = {
            'signals': ['data_loaded'],
            'sisa': ['data_loaded', 'signals_detected'],
            'shap': ['data_loaded', 'signals_detected', 'sisa_trained'],
            'sar': ['data_loaded', 'signals_detected'],  # Can work without ML
            'unlearning': ['data_loaded', 'signals_detected', 'sisa_trained']
        }
        
        if stage not in dependencies:
            return True
        
        required = dependencies[stage]
        return all(self.stages.get(dep, False) for dep in required)
    
    def mark_complete(self, stage):
        \"\"\"Mark a stage as complete\"\"\"
        self.stages[stage] = True
    
    def get_status_message(self, stage):
        \"\"\"Get user-friendly status message\"\"\"
        if self.can_proceed_to(stage):
            return "✅ Ready"
        
        # Find missing prerequisites
        dependencies = {
            'signals': 'Load data first',
            'sisa': 'Run signal detection first',
            'shap': 'Train SISA model first',
            'sar': 'Run signal detection first',
            'unlearning': 'Train SISA model first'
        }
        
        return f"⚠️ {dependencies.get(stage, 'Prerequisites not met')}"
    
    def reset_downstream(self, from_stage):
        \"\"\"Reset all stages downstream of a given stage (for unlearning)\"\"\"
        cascade = {
            'sisa_trained': ['shap_computed'],
            'signals_detected': ['sisa_trained', 'shap_computed']
        }
        
        if from_stage in cascade:
            for stage in cascade[from_stage]:
                self.stages[stage] = False
    
    def increment_model_version(self):
        \"\"\"Increment model version after unlearning\"\"\"
        self.data['model_version'] += 1
        return self.data['model_version']
