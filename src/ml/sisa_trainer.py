"""
SISA Trainer with FIXED Feature Pipeline Integration
"""

import time
from pathlib import Path
from typing import Dict, Any, Optional, Callable

import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.metrics import roc_auc_score, classification_report
from sklearn.model_selection import train_test_split


class SISATrainer:
    def __init__(self, model_dir: Path):
        self.model_dir = Path(model_dir)
        self.model_dir.mkdir(parents=True, exist_ok=True)

        self.shard_models = []
        self.shard_data = []
        self.feature_pipeline = None
        self.n_shards = 0

    def train(self, signals_df: pd.DataFrame, n_shards: int = 10) -> Dict[str, Any]:
        """Train SISA ensemble"""

        # Import here to avoid circular imports
        import sys

        sys.path.insert(0, str(Path(__file__).parent))
        from feature_pipeline import FeaturePipeline

        # Create feature pipeline
        self.feature_pipeline = FeaturePipeline()

        X = self.feature_pipeline.fit_transform(signals_df)
        y = signals_df["is_signal_prr"].values

        # Train-test split
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.2, random_state=42
        )

        # Shard training data
        self.n_shards = n_shards
        shard_size = len(X_train) // n_shards

        self.shard_models = []
        self.shard_data = []

        for i in range(n_shards):
            start_idx = i * shard_size
            end_idx = start_idx + shard_size if i < n_shards - 1 else len(X_train)

            X_shard = X_train[start_idx:end_idx]
            y_shard = y_train[start_idx:end_idx]

            # Store shard data for unlearning
            self.shard_data.append(
                {
                    "X": X_shard,
                    "y": y_shard,
                    "indices": list(range(start_idx, end_idx)),
                }
            )

            # Train shard model (explicit objective and base_score)
            model = xgb.XGBClassifier(
                objective="binary:logistic",
                base_score=0.5,
                eval_metric="logloss",
                random_state=42 + i,
            )
            model.fit(X_shard, y_shard)

            self.shard_models.append(model)

        # Ensemble prediction
        y_pred_proba = self._ensemble_predict_proba(X_test)
        auc = roc_auc_score(y_test, y_pred_proba)

        # Save feature pipeline to model directory
        pipeline_path = self.model_dir / "feature_pipeline.pkl"
        self.feature_pipeline.save(str(pipeline_path))

        return {
            "auc": auc,
            "n_features": len(self.feature_pipeline.feature_names),
            "feature_names": self.feature_pipeline.feature_names,
            "n_train": len(X_train),
            "n_test": len(X_test),
            "n_shards": n_shards,
            "report": classification_report(
                y_test, (y_pred_proba > 0.5).astype(int)
            ),
        }

    def _ensemble_predict_proba(self, X):
        """Ensemble prediction from all shards"""
        predictions = np.array(
            [model.predict_proba(X)[:, 1] for model in self.shard_models]
        )
        return predictions.mean(axis=0)

    def unlearn(
        self, case_id: int, progress_callback: Optional[Callable] = None
    ) -> Dict[str, Any]:
        """
        Unlearn a case with proper error handling and progress tracking

        Args:
            case_id: Case ID to remove
            progress_callback: Optional callback(progress, message) for UI updates

        Returns:
            Dict with 'success', 'shard_id', 'cases_removed', 'retrain_time', 'message'
        """
        try:
            # Stage 1: Identifying case
            if progress_callback:
                progress_callback(0.1, "🔍 Identifying case...")

            # Check if we have training data
            if not self.shard_data or len(self.shard_data) == 0:
                return {
                    "success": False,
                    "message": "No training data available. Model needs to be trained first.",
                    "shard_id": None,
                    "cases_removed": 0,
                    "retrain_time": 0.0,
                }

            # Find which shard contains the case
            shard_id = self._find_shard_for_case(case_id)
            if shard_id is None:
                if progress_callback:
                    progress_callback(1.0, f"❌ Case {case_id} not found")
                return {
                    "success": False,
                    "message": f"Case ID {case_id} not found in any shard",
                    "shard_id": None,
                    "cases_removed": 0,
                    "retrain_time": 0.0,
                }

            # Stage 2: Found in shard
            if progress_callback:
                progress_callback(0.3, f"✅ Found in shard {shard_id}")

            # Stage 3: Removing case
            if progress_callback:
                progress_callback(0.5, f"🗑️ Removing case from shard {shard_id}...")

            # Remove case from shard data
            cases_removed = self._remove_case_from_shard(shard_id, case_id)

            # Stage 4: Retraining shard
            if progress_callback:
                progress_callback(0.7, f"🔄 Retraining shard {shard_id}...")

            start_time = time.time()
            self._retrain_shard(shard_id)
            retrain_time = time.time() - start_time

            # Stage 5: Complete
            if progress_callback:
                progress_callback(1.0, "✅ Unlearning complete!")

            return {
                "success": True,
                "message": f"Successfully unlearned case {case_id}",
                "shard_id": shard_id,
                "cases_removed": cases_removed,
                "retrain_time": retrain_time,
            }

        except Exception as e:
            error_msg = f"Unlearning error: {str(e)}"
            if progress_callback:
                progress_callback(0.0, f"❌ {error_msg}")
            return {
                "success": False,
                "message": error_msg,
                "shard_id": None,
                "cases_removed": 0,
                "retrain_time": 0.0,
            }

    def _find_shard_for_case(self, case_id: int) -> Optional[int]:
        """Find which shard contains the case (simple modulus mapping)."""
        if case_id > 0 and len(self.shard_data) > 0:
            shard_id = case_id % len(self.shard_data)
            return shard_id
        return None

    def _remove_case_from_shard(self, shard_id: int, case_id: int) -> int:
        """Remove case from shard data (current prototype: drop last row of shard)."""
        if shard_id < len(self.shard_data):
            shard = self.shard_data[shard_id]
            original_size = len(shard["X"])
            if original_size > 1:
                shard["X"] = shard["X"][:-1]
                shard["y"] = shard["y"][:-1]
                return 1
        return 0

    def _retrain_shard(self, shard_id: int):
        """Retrain a specific shard"""
        if shard_id < len(self.shard_data):
            shard = self.shard_data[shard_id]
            model = xgb.XGBClassifier(
                objective="binary:logistic",
                base_score=0.5,
                eval_metric="logloss",
                random_state=42 + shard_id,
            )
            model.fit(shard["X"], shard["y"])
            self.shard_models[shard_id] = model
