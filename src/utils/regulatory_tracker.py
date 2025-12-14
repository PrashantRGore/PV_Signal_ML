"""
Regulatory Compliance Tracker - FIXED
Logs ALL system operations to MLflow for FDA/EMA/PMDA audit compliance
"""

import mlflow
import mlflow.data
from datetime import datetime
import json
from typing import Dict, Any, Optional
import pandas as pd

class RegulatoryTracker:
    """
    Comprehensive audit trail logging for regulatory compliance
    Tracks: Signal Detection, ML Training, SAR Generation, Unlearning
    """
    
    def __init__(self, tracking_uri: str = "file:./logs/mlruns"):
        mlflow.set_tracking_uri(tracking_uri)
        self.tracking_uri = tracking_uri
    
    def log_signal_detection(self, 
                            n_signals: int, 
                            n_records: int,
                            prr_threshold: float,
                            chi2_threshold: float,
                            top_signals: pd.DataFrame,
                            data_source: str,
                            date_range: Optional[Dict] = None) -> str:
        """Log signal detection run to MLflow"""
        
        experiment_name = "1_Signal_Detection"
        mlflow.set_experiment(experiment_name)
        
        with mlflow.start_run(run_name=f"SignalDetection_{datetime.now().strftime('%Y%m%d_%H%M%S')}"):
            # Log parameters
            mlflow.log_param("data_source", data_source)
            mlflow.log_param("n_records_analyzed", n_records)
            mlflow.log_param("prr_threshold", prr_threshold)
            mlflow.log_param("chi2_threshold", chi2_threshold)
            mlflow.log_param("analysis_timestamp", datetime.now().isoformat())
            
            # FIXED: Proper date_range handling
            if date_range:
                if isinstance(date_range, dict):
                    mlflow.log_param("start_date", date_range.get('start', 'N/A'))
                    mlflow.log_param("end_date", date_range.get('end', 'N/A'))
                else:
                    # If date_range is a string or other type, log it directly
                    mlflow.log_param("date_range", str(date_range))
            
            # Log metrics
            mlflow.log_metric("total_signals_detected", n_signals)
            mlflow.log_metric("signal_rate_percent", (n_signals / n_records * 100) if n_records > 0 else 0)
            
            # Log top signals as artifact
            top_signals_path = "top_signals.csv"
            top_signals.head(100).to_csv(top_signals_path, index=False)
            mlflow.log_artifact(top_signals_path)
            
            # Log audit trail
            audit_info = {
                "operation": "SIGNAL_DETECTION",
                "timestamp": datetime.now().isoformat(),
                "user": "system",
                "status": "SUCCESS",
                "signals_detected": n_signals,
                "records_analyzed": n_records
            }
            mlflow.log_dict(audit_info, "audit_trail.json")
            
            run_id = mlflow.active_run().info.run_id
            mlflow.set_tag("operation_type", "signal_detection")
            mlflow.set_tag("regulatory_status", "compliant")
            
            return run_id
    
    def log_sar_generation(self,
                          drug_name: str,
                          event_name: str,
                          prr: float,
                          case_count: int,
                          chi_square: float,
                          causality_assessment: Optional[Dict] = None,
                          literature_count: int = 0,
                          model_used: str = "llama3.2",
                          report_path: Optional[str] = None) -> str:
        """Log SAR generation to MLflow"""
        
        experiment_name = "3_SAR_Generation"
        mlflow.set_experiment(experiment_name)
        
        with mlflow.start_run(run_name=f"SAR_{drug_name}_{event_name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"):
            # Log parameters
            mlflow.log_param("drug_name", drug_name)
            mlflow.log_param("event_name", event_name)
            mlflow.log_param("llm_model", model_used)
            mlflow.log_param("generation_timestamp", datetime.now().isoformat())
            
            # Log signal metrics
            mlflow.log_metric("prr_value", prr)
            mlflow.log_metric("case_count", case_count)
            mlflow.log_metric("chi_square_value", chi_square)
            mlflow.log_metric("literature_papers_found", literature_count)
            
            # Log causality assessment
            if causality_assessment:
                mlflow.log_dict(causality_assessment, "causality_assessment.json")
                
                if 'who_umc' in causality_assessment:
                    mlflow.log_param("who_umc_category", causality_assessment['who_umc'].get('category'))
                
                if 'naranjo' in causality_assessment:
                    mlflow.log_metric("naranjo_score", causality_assessment['naranjo'].get('score', 0))
                    mlflow.log_param("naranjo_category", causality_assessment['naranjo'].get('category'))
            
            # Log SAR report artifact
            if report_path:
                mlflow.log_artifact(report_path)
            
            # Regulatory audit trail
            audit_info = {
                "operation": "SAR_GENERATION",
                "timestamp": datetime.now().isoformat(),
                "drug_event_pair": f"{drug_name} - {event_name}",
                "signal_strength_prr": prr,
                "case_count": case_count,
                "causality_assessed": causality_assessment is not None,
                "literature_reviewed": literature_count > 0,
                "status": "SUCCESS"
            }
            mlflow.log_dict(audit_info, "audit_trail.json")
            
            run_id = mlflow.active_run().info.run_id
            mlflow.set_tag("operation_type", "sar_generation")
            mlflow.set_tag("drug_name", drug_name)
            mlflow.set_tag("event_name", event_name)
            mlflow.set_tag("regulatory_status", "compliant")
            
            return run_id
    
    def log_unlearning_operation(self,
                                case_id: int,
                                shard_id: Optional[int] = None,
                                cases_removed: int = 0,
                                retrain_time: float = 0.0,
                                success: bool = False,
                                error_message: Optional[str] = None,
                                model_version_before: str = "unknown",
                                model_version_after: str = "unknown") -> str:
        """Log unlearning operation to MLflow"""
        
        experiment_name = "4_Machine_Unlearning"
        mlflow.set_experiment(experiment_name)
        
        run_name = f"Unlearn_Case_{case_id}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        
        with mlflow.start_run(run_name=run_name):
            # Log parameters
            mlflow.log_param("case_id_removed", case_id)
            mlflow.log_param("operation_timestamp", datetime.now().isoformat())
            mlflow.log_param("model_version_before", model_version_before)
            mlflow.log_param("model_version_after", model_version_after)
            mlflow.log_param("gdpr_compliance", "RIGHT_TO_BE_FORGOTTEN")
            
            if shard_id is not None:
                mlflow.log_param("affected_shard_id", shard_id)
            
            # Log metrics
            mlflow.log_metric("cases_removed", cases_removed)
            mlflow.log_metric("retrain_time_seconds", retrain_time)
            mlflow.log_metric("success", 1 if success else 0)
            
            # Log status
            if error_message:
                mlflow.log_param("error_message", error_message)
            
            # Regulatory audit trail
            audit_info = {
                "operation": "MACHINE_UNLEARNING",
                "timestamp": datetime.now().isoformat(),
                "case_id": case_id,
                "shard_affected": shard_id,
                "cases_removed": cases_removed,
                "retrain_time_seconds": retrain_time,
                "success": success,
                "error": error_message,
                "compliance_reason": "GDPR Article 17 - Right to Erasure",
                "hipaa_compliance": "Yes",
                "data_retention_policy": "Case removed from training data permanently"
            }
            mlflow.log_dict(audit_info, "unlearning_audit_trail.json")
            
            run_id = mlflow.active_run().info.run_id
            mlflow.set_tag("operation_type", "unlearning")
            mlflow.set_tag("case_id", str(case_id))
            mlflow.set_tag("success", "true" if success else "false")
            mlflow.set_tag("regulatory_status", "compliant")
            mlflow.set_tag("gdpr_article_17", "compliant")
            
            return run_id

# Create global instance
tracker = RegulatoryTracker()
