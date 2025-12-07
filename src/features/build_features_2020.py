import pandas as pd
import numpy as np
from pathlib import Path

BASE_DIR = Path(r"C:\Users\koreo\Downloads\pv-signal-ml")
INTERIM_DIR = BASE_DIR / "data" / "interim"
PROCESSED_DIR = BASE_DIR / "data" / "processed"
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

def main():
    print("[INFO] Loading pairs...")
    pairs_path = INTERIM_DIR / "pairs_2020.parquet"
    df = pd.read_parquet(pairs_path)
    print(f"[INFO] Loaded {len(df):,} pairs")

    # Normalize drug column
    drug_col = "DRUGNAME_NORM" if "DRUGNAME_NORM" in df.columns else "DRUGNAME_RAW"
    df = df.rename(columns={drug_col: "DRUGNAME_NORM"})
    
    df["REPORT_DATE"] = pd.to_datetime(df["REPORT_DATE"], errors="coerce")
    df = df.dropna(subset=["REPORT_DATE"])
    print(f"[INFO] After date cleaning: {len(df):,} pairs")

    # Monthly reference date
    df["reference_date"] = df["REPORT_DATE"].dt.to_period('M').dt.to_timestamp()
    
    print("[INFO] Aggregating monthly counts...")
    grp = df.groupby(["DRUGNAME_NORM", "pt", "reference_date"])
    monthly = grp.size().reset_index(name="A")
    
    total_month = monthly.groupby("reference_date")["A"].sum().rename("N")
    monthly = monthly.merge(total_month, on="reference_date", how="left")
    
    monthly = monthly.sort_values(["DRUGNAME_NORM", "pt", "reference_date"]).reset_index(drop=True)
    print(f"[INFO] Monthly shape: {monthly.shape}")

    # Rolling windows
    print("[INFO] Computing rolling windows...")
    for months, col in [(1, "A_1m"), (3, "A_3m"), (6, "A_6m"), (12, "A_12m")]:
        monthly[col] = monthly.groupby(["DRUGNAME_NORM", "pt"])["A"].rolling(
            window=months, min_periods=1
        ).sum().values

    for months, col in [(1, "N_1m"), (3, "N_3m"), (6, "N_6m"), (12, "N_12m")]:
        monthly[col] = monthly.groupby("pt")["N"].rolling(
            window=months, min_periods=1
        ).sum().values

    # Features
    print("[INFO] Computing features...")
    eps = 1e-6
    monthly["prr_6m"] = monthly["A_6m"] / (monthly["N_6m"] + eps)

    for w in ["1m", "3m", "6m", "12m"]:
        monthly[f"log1p_reports_count_{w}"] = np.log1p(monthly[f"A_{w}"])

    monthly["growth_3to6"] = monthly["A_6m"] - monthly["A_3m"]
    monthly["growth_6to12"] = monthly["A_12m"] - monthly["A_6m"]
    monthly["ps_frac_6m"] = monthly["A_6m"] / (monthly["N_6m"] + eps)

    # Future fields (optional)
    future = monthly.copy()
    future["reference_date"] = future["reference_date"] - pd.DateOffset(months=6)
    future = future.rename(columns={
        "A_6m": "future_A_6m", "N_6m": "future_N_6m", "prr_6m": "future_prr_6m"
    })[["DRUGNAME_NORM", "pt", "reference_date", "future_A_6m", "future_N_6m", "future_prr_6m"]]

    features = monthly.merge(future, on=["DRUGNAME_NORM", "pt", "reference_date"], how="left")

    # FIXED LABELS: Data-driven thresholds for your PRR distribution
    prr_thr = 0.1      # Top ~1% of PRR values (realistic signals)
    min_A = 3          # Minimum statistical reliability
    
    features["label_prr"] = (
        (features["prr_6m"] >= prr_thr) & (features["A_6m"] >= min_A)
    ).astype(int)

    out_path = PROCESSED_DIR / "features_2020_prr_labeled.parquet"
    features.to_parquet(out_path, index=False)
    
    print("Saved features to:", out_path)
    print("Shape:", features.shape)
    print("Label distribution:\n", features["label_prr"].value_counts())

if __name__ == "__main__":
    main()
