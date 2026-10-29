"""
Denoising Autoencoder & Multivariate Imputation for Missing Pharmacovigilance Data
==================================================================================
Replaces mock identity scaling with authentic multivariate missing data imputation.

Supported Imputation Engines:
1. MICE (Multivariate Imputation by Chained Equations via IterativeImputer):
   Statistically rigorous, fully deterministic, models each feature conditionally on all others.
2. Tabular DAE (PyTorch Neural Denoising Autoencoder):
   Trained with corruption noise (dropout/Gaussian) and reconstruction MSE loss.

Legal & Scientific Attribution:
- van Buuren S, Groothuis-Oudshoorn K. mice: Multivariate Imputation by Chained
  Equations in R. Journal of Statistical Software 2011; 45(3): 1-67. DOI: 10.18637/jss.v045.i03.
- Vincent P, Larochelle H, Bengio Y, Manzagol PA. Extracting and composing robust
  features with denoising autoencoders. In ICML 2008; 1096-1103. DOI: 10.1145/1390156.1390294.
- ICH E9 Statistical Principles for Clinical Trials (Section 5.3 Missing Values).
"""

from pathlib import Path
from typing import List, Optional, Union
import pickle
import numpy as np
import pandas as pd
from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import IterativeImputer, SimpleImputer
from sklearn.preprocessing import StandardScaler

try:
    import torch
    import torch.nn as nn
    import torch.optim as optim
    HAS_TORCH = True
except ImportError:
    HAS_TORCH = False

from src.utils.logger import setup_logger

logger = setup_logger(__name__)


if HAS_TORCH:
    class TabularDAENet(nn.Module):
        """PyTorch Feedforward Denoising Autoencoder for tabular features."""

        def __init__(self, input_dim: int, hidden_dim: Optional[int] = None, dropout_p: float = 0.2):
            super().__init__()
            hidden_dim = hidden_dim or max(4, input_dim // 2)
            self.encoder = nn.Sequential(
                nn.Dropout(p=dropout_p),
                nn.Linear(input_dim, hidden_dim),
                nn.BatchNorm1d(hidden_dim),
                nn.GELU()
            )
            self.decoder = nn.Sequential(
                nn.Linear(hidden_dim, input_dim)
            )

        def forward(self, x: torch.Tensor) -> torch.Tensor:
            encoded = self.encoder(x)
            decoded = self.decoder(encoded)
            return decoded


class DenoisingAutoencoder:
    """
    Multivariate Missing Data Imputer for Pharmacovigilance Signals.
    Supports MICE (IterativeImputer) and authentic PyTorch Tabular DAE.
    """

    def __init__(
        self,
        model_path: Optional[Union[Path, str]] = None,
        strategy: str = "mice",
        random_state: int = 42
    ):
        """
        Initialize imputer.

        Args:
            model_path: Path for saving or loading trained weights.
            strategy: 'mice' (IterativeImputer) or 'dae' (PyTorch Denoising Autoencoder).
            random_state: Random seed for deterministic reproducibility.
        """
        self.model_path = Path(model_path) if model_path else Path("models/dae/imputer_model.pkl")
        self.strategy = strategy.lower()
        self.random_state = random_state
        self.is_fitted = False
        self.feature_names: List[str] = []
        self.scaler = StandardScaler()
        self.model_version = "2.0.0"

        if self.strategy == "mice":
            self.imputer = IterativeImputer(
                max_iter=10,
                random_state=random_state,
                initial_strategy="median",
                skip_complete=True
            )
            self.dae_net = None
        elif self.strategy == "dae":
            if not HAS_TORCH:
                logger.warning("PyTorch unavailable; falling back to MICE (IterativeImputer)")
                self.strategy = "mice"
                self.imputer = IterativeImputer(max_iter=10, random_state=random_state, initial_strategy="median")
                self.dae_net = None
            else:
                self.imputer = SimpleImputer(strategy="median")
                self.dae_net = None
        else:
            self.imputer = SimpleImputer(strategy="median")
            self.dae_net = None

    def fit(self, X: pd.DataFrame, epochs: int = 40, batch_size: int = 64, lr: float = 1e-3) -> "DenoisingAutoencoder":
        """
        Fit imputer on clean or partially missing feature matrix.

        Args:
            X: Training DataFrame (numeric columns).
            epochs: Training epochs (for DAE).
            batch_size: Batch size (for DAE).
            lr: Learning rate (for DAE).
        """
        # Filter to numeric columns only
        numeric_df = X.select_dtypes(include=[np.number])
        if numeric_df.empty:
            raise ValueError("Input DataFrame contains no numeric columns for imputation")

        self.feature_names = numeric_df.columns.tolist()
        logger.info(f"Fitting {self.strategy.upper()} imputer on {len(numeric_df):,} samples ({len(self.feature_names)} features)")

        if self.strategy == "mice":
            self.imputer.fit(numeric_df)
            self.is_fitted = True
            logger.info("MICE IterativeImputer successfully fitted")
            return self

        # Strategy: DAE with PyTorch
        imputed_init = self.imputer.fit_transform(numeric_df)
        scaled_init = self.scaler.fit_transform(imputed_init)

        input_dim = len(self.feature_names)
        self.dae_net = TabularDAENet(input_dim=input_dim, dropout_p=0.2)
        optimizer = optim.Adam(self.dae_net.parameters(), lr=lr, weight_decay=1e-5)
        criterion = nn.MSELoss()

        tensor_x = torch.tensor(scaled_init, dtype=torch.float32)
        dataset = torch.utils.data.TensorDataset(tensor_x)
        loader = torch.utils.data.DataLoader(dataset, batch_size=min(batch_size, len(tensor_x)), shuffle=True)

        self.dae_net.train()
        for epoch in range(epochs):
            total_loss = 0.0
            for batch in loader:
                bx = batch[0]
                # Corrupt batch with additive Gaussian noise
                noise = torch.randn_like(bx) * 0.15
                corrupted = bx + noise

                optimizer.zero_grad()
                reconstructed = self.dae_net(corrupted)
                loss = criterion(reconstructed, bx)
                loss.backward()
                optimizer.step()
                total_loss += loss.item()

        self.dae_net.eval()
        self.is_fitted = True
        logger.info(f"PyTorch Tabular DAE fitted ({epochs} epochs, final loss: {total_loss / max(len(loader), 1):.5f})")
        return self

    def transform(self, X: pd.DataFrame) -> pd.DataFrame:
        """
        Impute missing values in feature matrix.

        Args:
            X: DataFrame with potential NaNs.

        Returns:
            DataFrame with imputed values.
        """
        if not self.is_fitted:
            raise ValueError("Imputer must be fitted before transform. Call fit() first.")

        res = X.copy()
        numeric_cols = [c for c in self.feature_names if c in res.columns]
        if not numeric_cols:
            return res

        missing_before = res[numeric_cols].isna().sum().sum()
        if missing_before == 0:
            return res

        # Ensure all fitted features are present
        res_numeric = res[numeric_cols].copy()
        for col in self.feature_names:
            if col not in res_numeric.columns:
                res_numeric[col] = np.nan
        res_numeric = res_numeric[self.feature_names]

        if self.strategy == "mice":
            imputed_array = self.imputer.transform(res_numeric)
            imputed_df = pd.DataFrame(imputed_array, columns=self.feature_names, index=res.index)
        else:
            # PyTorch DAE transform
            init_imputed = self.imputer.transform(res_numeric)
            scaled = self.scaler.transform(init_imputed)
            with torch.no_grad():
                t_in = torch.tensor(scaled, dtype=torch.float32)
                t_out = self.dae_net(t_in).numpy()
            unscaled = self.scaler.inverse_transform(t_out)
            # Replace only the originally missing entries
            was_missing = res_numeric.isna().values
            final_array = np.where(was_missing, unscaled, res_numeric.values)
            # Fallback for any lingering NaNs
            nan_mask = np.isnan(final_array)
            if nan_mask.any():
                final_array = np.where(nan_mask, init_imputed, final_array)
            imputed_df = pd.DataFrame(final_array, columns=self.feature_names, index=res.index)

        for col in numeric_cols:
            res[col] = imputed_df[col]

        missing_after = res[numeric_cols].isna().sum().sum()
        logger.info(f"Imputation complete: {missing_before:,} NaNs -> {missing_after:,} NaNs remaining")
        return res

    def fit_transform(self, X: pd.DataFrame) -> pd.DataFrame:
        """Fit on X and transform X."""
        return self.fit(X).transform(X)

    def save(self, path: Optional[Union[Path, str]] = None) -> None:
        """Save fitted model to disk."""
        target = Path(path) if path else self.model_path
        target.parent.mkdir(parents=True, exist_ok=True)
        with open(target, "wb") as f:
            pickle.dump(self, f)
        logger.info(f"Saved DAE model to {target}")

    @classmethod
    def load(cls, path: Union[Path, str]) -> "DenoisingAutoencoder":
        """Load fitted model from disk."""
        with open(path, "rb") as f:
            model = pickle.load(f)
        logger.info(f"Loaded DAE model from {path}")
        return model
