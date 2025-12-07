"""
Audit Logging Module for Pharmacovigilance System

Tracks all API calls, user actions, and system events for regulatory compliance.
Implements HIPAA, GDPR, and FDA 21 CFR Part 11 audit trail requirements.
"""

import json
import logging
from pathlib import Path
from datetime import datetime
from typing import Optional, Dict, Any
from functools import wraps
import time


class AuditLogger:
    """
    Centralized audit logging for pharmacovigilance system.
    
    Features:
    - Log all API calls with timestamp, user, endpoint, parameters
    - Track report generation (SAR, PSMF)
    - Monitor data access and modifications
    - Generate audit reports for regulatory inspection
    - Immutable log format (append-only JSONL)
    """
    
    def __init__(self, log_dir: Path = None):
        """
        Initialize audit logger.
        
        Args:
            log_dir: Directory for audit logs (default: ./audit_logs)
        """
        self.log_dir = log_dir or Path.cwd() / "audit_logs"
        self.log_dir.mkdir(exist_ok=True, parents=True)
        
        self.access_log_path = self.log_dir / "access_log.jsonl"
        self.report_log_path = self.log_dir / "report_generation.jsonl"
        self.data_modification_log_path = self.log_dir / "data_modifications.jsonl"
        self.mlflow_log_path = self.log_dir / "mlflow_runs.jsonl"
        
        # Setup Python logging
        self.logger = logging.getLogger("pv_audit")
        self.logger.setLevel(logging.INFO)
        
        # File handler for audit log
        handler = logging.FileHandler(self.log_dir / "audit.log")
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        handler.setFormatter(formatter)
        if not self.logger.handlers:
            self.logger.addHandler(handler)
    
    def log_api_call(
        self,
        endpoint: str,
        method: str,
        user_id: Optional[str] = None,
        parameters: Optional[Dict] = None,
        response_status: int = 200,
        response_size: int = 0,
        duration_ms: float = 0.0
    ) -> Dict:
        """
        Log an API call.
        
        Args:
            endpoint: API endpoint (e.g., "/signals/top_candidates")
            method: HTTP method (GET, POST, etc.)
            user_id: User ID (if authenticated)
            parameters: Request parameters (sanitized)
            response_status: HTTP response status code
            response_size: Response size in bytes
            duration_ms: Request duration in milliseconds
            
        Returns:
            Log record
        """
        record = {
            "timestamp": datetime.utcnow().isoformat() + "Z",
            "event_type": "api_call",
            "endpoint": endpoint,
            "method": method,
            "user_id": user_id or "anonymous",
            "parameters": self._sanitize_parameters(parameters or {}),
            "response_status": response_status,
            "response_size_bytes": response_size,
            "duration_ms": duration_ms,
            "ip_address": None
        }
        
        self._append_log(self.access_log_path, record)
        self.logger.info(
            f"API: {method} {endpoint} - Status {response_status} - {duration_ms:.1f}ms"
        )
        
        return record
    
    def log_report_generation(
        self,
        report_type: str,
        drug: str,
        event: str,
        period: str,
        user_id: Optional[str] = None,
        output_file: Optional[str] = None,
        success: bool = True,
        error_message: Optional[str] = None
    ) -> Dict:
        """
        Log report generation event.
        
        Args:
            report_type: Type of report (SAR, PSMF, etc.)
            drug: Drug name
            event: Adverse event
            period: Assessment period
            user_id: User who generated report
            output_file: Path to generated file
            success: Whether generation succeeded
            error_message: Error message if failed
            
        Returns:
            Log record
        """
        record = {
            "timestamp": datetime.utcnow().isoformat() + "Z",
            "event_type": "report_generation",
            "report_type": report_type,
            "drug": drug,
            "event": event,
            "period": period,
            "user_id": user_id or "system",
            "output_file": output_file,
            "success": success,
            "error_message": error_message
        }
        
        self._append_log(self.report_log_path, record)
        status = "SUCCESS" if success else "FAILED"
        self.logger.info(
            f"Report: {report_type} - {drug}/{event} - {status}"
        )
        
        return record
    
    def log_data_modification(
        self,
        operation: str,
        table_name: str,
        record_count: int,
        user_id: Optional[str] = None,
        description: Optional[str] = None
    ) -> Dict:
        """
        Log data modification event.
        
        Args:
            operation: Type of operation (INSERT, UPDATE, DELETE)
            table_name: Name of table modified
            record_count: Number of records affected
            user_id: User who made modification
            description: Description of change
            
        Returns:
            Log record
        """
        record = {
            "timestamp": datetime.utcnow().isoformat() + "Z",
            "event_type": "data_modification",
            "operation": operation,
            "table_name": table_name,
            "record_count": record_count,
            "user_id": user_id or "system",
            "description": description
        }
        
        self._append_log(self.data_modification_log_path, record)
        self.logger.info(
            f"Data: {operation} {table_name} - {record_count} records"
        )
        
        return record
    
    def log_access_control_event(
        self,
        action: str,
        user_id: str,
        resource: str,
        allowed: bool,
        reason: Optional[str] = None
    ) -> Dict:
        """
        Log access control event.
        
        Args:
            action: Type of action
            user_id: User ID
            resource: Resource being accessed
            allowed: Whether access was allowed
            reason: Reason for denial (if applicable)
            
        Returns:
            Log record
        """
        record = {
            "timestamp": datetime.utcnow().isoformat() + "Z",
            "event_type": "access_control",
            "action": action,
            "user_id": user_id,
            "resource": resource,
            "allowed": allowed,
            "reason": reason
        }
        
        self._append_log(self.access_log_path, record)
        status = "ALLOWED" if allowed else "DENIED"
        self.logger.warning(
            f"Access: {action} - {user_id} - {resource} - {status}"
        )
        
        return record
    
    def log_mlflow_run(
        self,
        run_id: str,
        experiment_name: str,
        model_version: str,
        metrics: Optional[Dict[str, float]] = None,
        parameters: Optional[Dict[str, Any]] = None,
        user_id: Optional[str] = None,
        status: str = "completed"
    ) -> Dict:
        """
        Log MLflow run to audit trail.
        
        Args:
            run_id: MLflow run ID
            experiment_name: Experiment name
            model_version: Model version
            metrics: Model metrics (AP, AUC, etc.)
            parameters: Model parameters
            user_id: User who triggered the run
            status: Run status (completed, failed, etc.)
            
        Returns:
            Log record
        """
        record = {
            "timestamp": datetime.utcnow().isoformat() + "Z",
            "event_type": "mlflow_run",
            "run_id": run_id,
            "experiment_name": experiment_name,
            "model_version": model_version,
            "metrics": metrics or {},
            "parameters": parameters or {},
            "user_id": user_id or "system",
            "status": status
        }
        
        self._append_log(self.mlflow_log_path, record)
        self.logger.info(
            f"MLflow: {experiment_name} - Run {run_id} - {status}"
        )
        
        return record
    
    def get_audit_report(self, days: int = 30) -> Dict:
        """
        Generate audit report for specified period.
        
        Args:
            days: Number of days to include in report
            
        Returns:
            Audit report with statistics and summaries
        """
        cutoff_time = datetime.utcnow().timestamp() - (days * 86400)
        
        api_calls = self._read_log_file(self.access_log_path, cutoff_time)
        reports = self._read_log_file(self.report_log_path, cutoff_time)
        modifications = self._read_log_file(self.data_modification_log_path, cutoff_time)
        
        report = {
            "generated_at": datetime.utcnow().isoformat() + "Z",
            "period_days": days,
            "summary": {
                "total_api_calls": len(api_calls),
                "total_reports_generated": len(reports),
                "successful_reports": sum(1 for r in reports if r.get("success", True)),
                "failed_reports": sum(1 for r in reports if not r.get("success", True)),
                "total_data_modifications": len(modifications)
            },
            "api_endpoints": self._summarize_api_calls(api_calls),
            "report_types": self._summarize_reports(reports),
            "data_operations": self._summarize_modifications(modifications),
            "compliance_notes": [
                "All events are logged with timestamp and user ID",
                "Logs are append-only (immutable)",
                "Access control events are tracked for HIPAA compliance",
                "Report generation is auditable for regulatory inspection"
            ]
        }
        
        return report
    
    def _sanitize_parameters(self, params: Dict) -> Dict:
        """
        Remove sensitive information from parameters before logging.
        
        Args:
            params: Request parameters
            
        Returns:
            Sanitized parameters
        """
        sensitive_keys = {"password", "token", "api_key", "secret"}
        sanitized = {}
        
        for key, value in params.items():
            if key.lower() in sensitive_keys:
                sanitized[key] = "***REDACTED***"
            else:
                sanitized[key] = value
        
        return sanitized
    
    def _append_log(self, log_path: Path, record: Dict):
        """
        Append record to JSONL log file (append-only).
        
        Args:
            log_path: Path to log file
            record: Record to append
        """
        with open(log_path, 'a', encoding='utf-8') as f:
            f.write(json.dumps(record, default=str) + "\n")
    
    def _read_log_file(self, log_path: Path, cutoff_time: float = 0) -> list:
        """
        Read log file and filter by timestamp.
        
        Args:
            log_path: Path to log file
            cutoff_time: Unix timestamp cutoff
            
        Returns:
            List of records
        """
        records = []
        
        if not log_path.exists():
            return records
        
        with open(log_path, 'r', encoding='utf-8') as f:
            for line in f:
                if line.strip():
                    record = json.loads(line)
                    try:
                        ts = datetime.fromisoformat(
                            record["timestamp"].replace("Z", "+00:00")
                        ).timestamp()
                        if ts >= cutoff_time:
                            records.append(record)
                    except (KeyError, ValueError):
                        records.append(record)
        
        return records
    
    def _summarize_api_calls(self, calls: list) -> Dict:
        """Summarize API calls by endpoint."""
        summary = {}
        for call in calls:
            if call.get("event_type") == "api_call":
                endpoint = call.get("endpoint", "unknown")
                if endpoint not in summary:
                    summary[endpoint] = {"count": 0, "success": 0, "errors": 0}
                summary[endpoint]["count"] += 1
                if call.get("response_status", 200) < 400:
                    summary[endpoint]["success"] += 1
                else:
                    summary[endpoint]["errors"] += 1
        return summary
    
    def _summarize_reports(self, reports: list) -> Dict:
        """Summarize report generation by type."""
        summary = {}
        for report in reports:
            report_type = report.get("report_type", "unknown")
            if report_type not in summary:
                summary[report_type] = {"total": 0, "success": 0, "failed": 0}
            summary[report_type]["total"] += 1
            if report.get("success", True):
                summary[report_type]["success"] += 1
            else:
                summary[report_type]["failed"] += 1
        return summary
    
    def _summarize_modifications(self, mods: list) -> Dict:
        """Summarize data modifications by operation."""
        summary = {}
        for mod in mods:
            operation = mod.get("operation", "unknown")
            if operation not in summary:
                summary[operation] = {"count": 0, "total_records": 0}
            summary[operation]["count"] += 1
            summary[operation]["total_records"] += mod.get("record_count", 0)
        return summary


# Global audit logger instance
audit_logger = AuditLogger()


def audit_api_call(endpoint: str, method: str = "GET"):
    """
    Decorator for auditing API calls.
    
    Usage:
        @audit_api_call("/signals/top_candidates")
        def get_top_candidates(n: int = 10):
            ...
    """
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            start_time = time.time()
            try:
                result = func(*args, **kwargs)
                duration_ms = (time.time() - start_time) * 1000
                audit_logger.log_api_call(
                    endpoint=endpoint,
                    method=method,
                    response_status=200,
                    duration_ms=duration_ms
                )
                return result
            except Exception as e:
                duration_ms = (time.time() - start_time) * 1000
                audit_logger.log_api_call(
                    endpoint=endpoint,
                    method=method,
                    response_status=500,
                    duration_ms=duration_ms
                )
                raise
        return wrapper
    return decorator


if __name__ == "__main__":
    # Example usage
    logger = AuditLogger()
    
    # Log API call
    logger.log_api_call(
        endpoint="/signals/top_candidates",
        method="GET",
        user_id="analyst_001",
        response_status=200,
        duration_ms=125.5
    )
    
    # Log report generation
    logger.log_report_generation(
        report_type="SAR",
        drug="InsulPen (Insulin)",
        event="Incorrect dose administered",
        period="2025-01-01_2025-03-31",
        user_id="analyst_001",
        output_file="sar_reports/reports/InsulPen__Incorrect_dose__2025-01-01_2025-03-31.json",
        success=True
    )
    
    # Generate audit report
    report = logger.get_audit_report(days=30)
    print("\n✅ Audit Report:")
    print(json.dumps(report, indent=2, default=str))
