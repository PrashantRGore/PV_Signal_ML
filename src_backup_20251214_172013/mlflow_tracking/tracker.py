import mlflow
import mlflow.sklearn
from datetime import datetime
import json
import os

class PVMLflowTracker:
    def __init__(self, tracking_uri='file:./mlruns', experiment_name='PV_Signal_Detection'):
        mlflow.set_tracking_uri(tracking_uri)
        mlflow.set_experiment(experiment_name)
        self.experiment_name = experiment_name
        self.current_run_id = None
        self.model_version = 1
    
    def start_pipeline_run(self, run_name=None):
        if run_name is None:
            run_name = f'Pipeline_{datetime.now().strftime("%Y%m%d_%H%M%S")}'
        mlflow.start_run(run_name=run_name)
        self.current_run_id = mlflow.active_run().info.run_id
        mlflow.log_param('pipeline_start', datetime.now().isoformat())
        mlflow.log_param('model_version', self.model_version)
        return self.current_run_id
    
    def log_signal_detection(self, signals_df, parameters):
        mlflow.log_params(parameters)
        mlflow.log_metric('total_signals', len(signals_df))
        mlflow.log_metric('avg_prr', float(signals_df['prr'].mean()))
        mlflow.log_metric('max_chi2', float(signals_df['chi2'].max()))
        return len(signals_df)
    
    def log_sisa_training(self, training_results, model=None):
        mlflow.log_metric('sisa_auc', training_results.get('auc', 0))
        mlflow.log_metric('sisa_accuracy', training_results.get('accuracy', 0))
        mlflow.log_metric('n_train_samples', training_results.get('n_train', 0))
        mlflow.set_tag('sisa_trained', 'true')
        return training_results.get('auc', 0)
    
    def log_shap_analysis(self, shap_results):
        mlflow.log_dict(shap_results.get('top_features', {}), 'shap/top_features.json')
        mlflow.set_tag('shap_computed', 'true')
        mlflow.log_param('shap_timestamp', datetime.now().isoformat())
        return True
    
    def log_sar_generation(self, sar_result):
        with mlflow.start_run(nested=True, run_name=f'SAR_{sar_result["signal_data"]["drug_name"][:20]}'):
            mlflow.log_param('drug', sar_result['signal_data']['drug_name'])
            mlflow.log_param('event', sar_result['signal_data']['event_name'])
            mlflow.log_param('llm_model', sar_result['model_used'])
            mlflow.log_metric('prr', float(sar_result['signal_data']['prr']))
            mlflow.log_metric('chi_square', float(sar_result['signal_data']['chi_square']))
            mlflow.log_metric('case_count', int(sar_result['signal_data']['count']))
            sar_run_id = mlflow.active_run().info.run_id
        return sar_run_id
    
    def log_unlearning_event(self, unlearn_result):
        new_run_name = f'Unlearning_{unlearn_result.get("case_id", "unknown")}_{datetime.now().strftime("%H%M%S")}'
        with mlflow.start_run(run_name=new_run_name):
            mlflow.log_param('unlearning_event', 'true')
            mlflow.log_param('case_id_removed', unlearn_result.get('case_id', 'unknown'))
            mlflow.log_param('timestamp', datetime.now().isoformat())
            self.model_version += 1
            mlflow.log_param('model_version_after', self.model_version)
            mlflow.set_tag('unlearning_completed', 'true')
            unlearn_run_id = mlflow.active_run().info.run_id
        return unlearn_run_id
    
    def end_pipeline_run(self):
        if mlflow.active_run():
            mlflow.log_param('pipeline_end', datetime.now().isoformat())
            mlflow.end_run()
            self.current_run_id = None
    
    def get_all_runs(self):
        client = mlflow.tracking.MlflowClient()
        experiment = client.get_experiment_by_name(self.experiment_name)
        if experiment:
            runs = client.search_runs(experiment.experiment_id, order_by=['start_time DESC'])
            return runs
        return []
    
    def get_model_version_history(self):
        runs = self.get_all_runs()
        history = []
        for run in runs:
            if run.data.params.get('model_version'):
                history.append({
                    'version': run.data.params.get('model_version'),
                    'timestamp': datetime.fromtimestamp(run.info.start_time/1000),
                    'run_id': run.info.run_id[:8],
                    'unlearning_event': run.data.params.get('unlearning_event', 'false'),
                    'auc': run.data.metrics.get('sisa_auc', 'N/A')
                })
        return history
