from pathlib import Path
import sys
import zipfile
import requests
import pandas as pd
import numpy as np

BASE_DIR = Path(r"C:\Users\koreo\Downloads\pv-signal-ml")
FAERS_DIR = BASE_DIR / "faers_raw"
SAR_DIR = BASE_DIR / "sar_reports"

FAERS_DIR.mkdir(exist_ok=True, parents=True)


def _quarter_from_date(date_str: str):
    """Return (year, quarter_str) for 'YYYY-MM-DD'."""
    year, month, _ = map(int, date_str.split("-"))
    if month <= 3:
        q = "Q1"
    elif month <= 6:
        q = "Q2"
    elif month <= 9:
        q = "Q3"
    else:
        q = "Q4"
    return year, q


def _faers_quarters_between(start_date: str, end_date: str):
    """
    Simple implementation: return all quarters between start and end years.
    You can tighten this later if needed.
    """
    sy, _ = _quarter_from_date(start_date)
    ey, _ = _quarter_from_date(end_date)
    qs = []
    for y in range(sy, ey + 1):
        for q in ["Q1", "Q2", "Q3", "Q4"]:
            qs.append((y, q))
    return qs


def _download_quarter(year: int, quarter: str):
    """
    Download FAERS quarterly ZIP for year, quarter into FAERS_DIR if not present.

    Uses FDA public URLs (ASCII exports).[web:98][web:100]
    """
    fname = f"faers_ascii_{year}{quarter}.zip"
    out_path = FAERS_DIR / fname
    if out_path.exists():
        return out_path

    url = f"https://fis.fda.gov/content/Exports/{fname}"
    print(f"Downloading {url} ...")
    r = requests.get(url, timeout=300)
    r.raise_for_status()
    with open(out_path, "wb") as f:
        f.write(r.content)
    print(f"✅ Saved {out_path}")
    return out_path


def _load_drug_reac_for_quarters(quarters, start_date: str, end_date: str):
    """
    Load DRUG and REAC tables for given quarters and filter by REPT_DT date window.
    Uses a tolerant parser for FAERS ASCII (malformed quotes, etc.).
    """
    all_drug = []
    all_reac = []

    def _read_faers_table(zf, member_name):
        # Read raw bytes, decode loosely, clean, then parse with '$' separator.
        raw = zf.read(member_name)
        text = raw.decode("latin-1", errors="ignore")
        text = (
            text.replace("\r\n", "\n")
            .replace("\r", "\n")
            .replace("\x00", "")
        )
        from io import StringIO

        buf = StringIO(text)
        df = pd.read_csv(
            buf,
            sep="$",
            dtype=str,
            engine="python",
            on_bad_lines="skip",  # skip malformed lines instead of raising
        )
        return df

    for year, q in quarters:
        zip_path = _download_quarter(year, q)
        with zipfile.ZipFile(zip_path, "r") as zf:
            names = zf.namelist()
            drug_name = [n for n in names if "DRUG" in n.upper()][0]
            reac_name = [n for n in names if "REAC" in n.upper()][0]
            demo_name = [n for n in names if "DEMO" in n.upper()][0]

            df_drug = _read_faers_table(zf, drug_name)
            df_reac = _read_faers_table(zf, reac_name)
            df_demo = _read_faers_table(zf, demo_name)

            # Normalize date and filter by REPT_DT
            if "rept_dt" in df_demo.columns:
                df_demo["REPT_DT"] = pd.to_datetime(
                    df_demo["rept_dt"], errors="coerce", format="%Y%m%d"
                )
            elif "REPT_DT" in df_demo.columns:
                df_demo["REPT_DT"] = pd.to_datetime(
                    df_demo["REPT_DT"], errors="coerce", format="%Y%m%d"
                )
            else:
                df_demo["REPT_DT"] = pd.NaT

            mask = (df_demo["REPT_DT"] >= start_date) & (
                df_demo["REPT_DT"] <= end_date
            )
            demo_ids = df_demo.loc[mask, "caseid"].dropna().unique()

            if "caseid" in df_drug.columns:
                df_drug = df_drug[df_drug["caseid"].isin(demo_ids)]
            if "caseid" in df_reac.columns:
                df_reac = df_reac[df_reac["caseid"].isin(demo_ids)]

            all_drug.append(df_drug)
            all_reac.append(df_reac)

    if not all_drug or not all_reac:
        return pd.DataFrame(), pd.DataFrame()

    drug = pd.concat(all_drug, ignore_index=True)
    reac = pd.concat(all_reac, ignore_index=True)
    return drug, reac


def _compute_prr(drug_reac):
    """
    Compute PRR and chi-square for DRUG–EVENT pairs with exact contingency counts.

    drug_reac columns: DRUG, EVENT, CASEID.
    """
    total_cases = drug_reac["CASEID"].nunique()

    # a: cases with drug & event
    a_df = (
        drug_reac.groupby(["DRUG", "EVENT"])["CASEID"]
        .nunique()
        .reset_index(name="a")
    )

    # Cases with drug (any event)
    drug_cases = (
        drug_reac.groupby("DRUG")["CASEID"].nunique().reset_index(name="drug_cases")
    )

    # Cases with event (any drug)
    event_cases = (
        drug_reac.groupby("EVENT")["CASEID"].nunique().reset_index(name="event_cases")
    )

    df = a_df.merge(drug_cases, on="DRUG").merge(event_cases, on="EVENT")

    # b: cases with drug & other events
    df["b"] = df["drug_cases"] - df["a"]
    # c: cases with event & other drugs
    df["c"] = df["event_cases"] - df["a"]
    # d: all other cases
    df["d"] = total_cases - (df["a"] + df["b"] + df["c"])

    a = df["a"].astype(float)
    b = df["b"].astype(float)
    c = df["c"].astype(float)
    d = df["d"].astype(float)

    num = a / (a + b + 1e-9)
    den = c / (c + d + 1e-9)
    prr = num / (den + 1e-9)

    chisq = (a + b + c + d) * (a * d - b * c) ** 2 / (
        (a + b) * (c + d) * (a + c) * (b + d) + 1e-9
    )

    df["PRR"] = prr.replace([np.inf, -np.inf], np.nan).fillna(0.0)
    df["CHISQ"] = chisq.fillna(0.0)
    df["CASES"] = df["a"].astype(int)

    # Basic threshold like common signal rules
    df = df[df["CASES"] >= 3]

    return df[["DRUG", "EVENT", "CASES", "PRR", "CHISQ"]]


def build_signals_for_period(start_date: str, end_date: str):
    quarters = _faers_quarters_between(start_date, end_date)
    drug, reac = _load_drug_reac_for_quarters(quarters, start_date, end_date)
    if drug.empty or reac.empty:
        raise ValueError(f"No FAERS records for {start_date} to {end_date}")

    # join DRUG and REAC by CASEID to get DRUG–EVENT
    drug_col = "prod_ai" if "prod_ai" in drug.columns else "drugname"
    event_col = "pt" if "pt" in reac.columns else "PT"

    dr = drug[["caseid", drug_col]].rename(columns={drug_col: "DRUG"})
    re = reac[["caseid", event_col]].rename(columns={event_col: "EVENT"})
    long = dr.merge(re, on="caseid", how="inner").dropna(subset=["DRUG", "EVENT"])
    long["CASEID"] = long["caseid"]

    signals = _compute_prr(long)

    period_tag = f"{start_date}_{end_date}"
    out_path = SAR_DIR / f"full_signals_{period_tag}.csv"
    signals.to_csv(out_path, index=False)
    print(f"✅ Saved FAERS period signals to {out_path}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python faers_build_signals.py YYYY-MM-DD YYYY-MM-DD")
        sys.exit(1)
    start, end = sys.argv[1], sys.argv[2]
    build_signals_for_period(start, end)
