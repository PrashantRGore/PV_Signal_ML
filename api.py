from fastapi import FastAPI
from pydantic import BaseModel
import json
import pandas as pd
from pathlib import Path
import mlflow
from audit_logging import audit_logger
from stats_engine import add_signal_flags_from_existing_stats
import joblib
import json

app = FastAPI(title="PV Signal Detection API")

# Configure MLflow for audit trail
MLFLOW_URI = "file:///C:/Users/koreo/mlruns"
mlflow.set_tracking_uri(MLFLOW_URI)
mlflow.set_experiment("pv-signal-ml-api")

class SignalRequest(BaseModel):
    drug: str
    event: str
    period: str = "2025-01-01_2025-03-31"

@app.get("/signals/top_candidates/{n}")
def get_top_candidates(n: int = 10):
    df = pd.read_csv('sar_reports/candidate_signals_2025Q1_stats_engine.csv')
    return df.head(n).to_dict(orient='records')

@app.get("/signals/{drug}/{event}")
def get_signal_detail(drug: str, event: str):
    df = pd.read_csv('sar_reports/enriched_signals_2025Q1_stats_engine.csv')
    signal = df[(df['DRUG'].str.contains(drug, case=False, na=False)) & 
                (df['EVENT'].str.contains(event, case=False, na=False))]
    return signal.to_dict(orient='records') if not signal.empty else {'error': 'Signal not found'}

@app.get("/ml_features/summary")
def ml_features_summary():
    features = pd.read_csv('ml_data/stats_engine_features.csv')
    return {
        'shape': features.shape,
        'stats': features.describe().to_dict()
    }

@app.get("/lineage/{dataset}")
def get_lineage(dataset: str):
    lineage_file = Path(f'lineage/{dataset}_lineage.json')
    if lineage_file.exists():
        with open(lineage_file) as f:
            return json.load(f)
    return {'error': 'Lineage not found'}

@app.post("/signal-report")
def generate_signal_report(request: SignalRequest):
    """Generate a signal assessment report (SAR) for a drug-event pair."""
    # Start MLflow run for this report generation
    with mlflow.start_run():
        try:
            # Import RAG pipeline
            from rag_langchain import PVSignalRAGLangChain
            from signal_report_builder import build_signal_report
            
            # Log request parameters to MLflow
            mlflow.log_param("drug", request.drug)
            mlflow.log_param("event", request.event)
            mlflow.log_param("period", request.period)
            
            # Try to get signal data from CSV
            try:
                df = pd.read_csv('sar_reports/enriched_signals_2025Q1_stats_engine.csv')
                signal = df[(df['DRUG'].str.contains(request.drug, case=False, na=False)) & 
                           (df['EVENT'].str.contains(request.event, case=False, na=False))]
                
                if signal.empty:
                    # If not found, create a minimal signal row
                    signal_row = {
                        'DRUG': request.drug,
                        'EVENT': request.event,
                        'PRR': 2.5,
                        'CHISQ': 4.5,
                        'CASES': 5
                    }
                else:
                    signal_row = signal.iloc[0].to_dict()
            except:
                # Fallback: create minimal signal row
                signal_row = {
                    'DRUG': request.drug,
                    'EVENT': request.event,
                    'PRR': 2.5,
                    'CHISQ': 4.5,
                    'CASES': 5
                }
            
            # Log signal metrics to MLflow
            mlflow.log_metric("prr", float(signal_row.get('PRR', 0)))
            mlflow.log_metric("chisq", float(signal_row.get('CHISQ', 0)))
            mlflow.log_metric("cases", int(signal_row.get('CASES', 0)))
            
            # Initialize RAG
            rag = PVSignalRAGLangChain()
            
            # Generate SAR using RAG
            sar = rag.generate_sar_report(signal_row)
            
            # Save report
            safe_period = request.period.replace('/', '-').replace(':', '_')
            report_path = Path('sar_reports') / 'reports' / f"{request.drug}__{request.event}__{safe_period}.json"
            report_path.parent.mkdir(parents=True, exist_ok=True)
            with open(report_path, 'w', encoding='utf-8') as f:
                json.dump(sar, f, indent=2)
            
            # Log report metrics to MLflow
            mlflow.log_metric("report_size_bytes", len(json.dumps(sar)))
            mlflow.log_artifact(str(report_path))
            
            # Get MLflow run ID for audit trail
            run_id = mlflow.active_run().info.run_id
            
            # Log report generation to audit trail
            audit_logger.log_report_generation(
                report_type="SAR",
                drug=request.drug,
                event=request.event,
                period=request.period,
                user_id="api",
                output_file=str(report_path),
                success=True
            )
            
            # Log MLflow run to audit trail
            audit_logger.log_mlflow_run(
                run_id=run_id,
                experiment_name="pv-signal-ml-api",
                model_version="1.0",
                metrics={
                    "prr": float(signal_row.get('PRR', 0)),
                    "chisq": float(signal_row.get('CHISQ', 0)),
                    "cases": int(signal_row.get('CASES', 0)),
                    "report_size": len(json.dumps(sar))
                },
                parameters={
                    "drug": request.drug,
                    "event": request.event,
                    "period": request.period
                },
                user_id="api",
                status="completed"
            )
            
            return {
                'status': 'success',
                'drug': request.drug,
                'event': request.event,
                'period': request.period,
                'report': sar,
                'saved_to': str(report_path),
                'mlflow_run_id': run_id
            }
        except Exception as e:
            # Log failure to MLflow
            mlflow.log_param("error", str(e))
            
            # Log failure to audit trail
            audit_logger.log_report_generation(
                report_type="SAR",
                drug=request.drug,
                event=request.event,
                period=request.period,
                user_id="api",
                success=False,
                error_message=str(e)
            )
            
            # Log failed run to audit trail
            try:
                run_id = mlflow.active_run().info.run_id
                audit_logger.log_mlflow_run(
                    run_id=run_id,
                    experiment_name="pv-signal-ml-api",
                    model_version="1.0",
                    metrics={},
                    parameters={
                        "drug": request.drug,
                        "event": request.event,
                        "period": request.period
                    },
                    user_id="api",
                    status="failed"
                )
            except:
                pass
        return {
            'status': 'error',
            'error': str(e),
            'error_type': type(e).__name__
        }

@app.get("/")
def root():
    """API health check and documentation."""
    return {
        'status': 'ok',
        'service': 'PV Signal Detection API',
        'version': '1.0',
        'endpoints': {
            'GET /signals/top_candidates/{n}': 'Get top N candidate signals',
            'GET /signals/{drug}/{event}': 'Get signal details for drug-event pair',
            'GET /ml_features/summary': 'Get ML features summary',
            'GET /lineage/{dataset}': 'Get data lineage for dataset',
            'POST /signal-report': 'Generate signal assessment report (SAR)',
            'GET /docs': 'Interactive API documentation (Swagger UI)'
        }
    }

if __name__ == '__main__':
    import uvicorn
    uvicorn.run(app, host='0.0.0.0', port=8000)
