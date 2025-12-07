from pathlib import Path
import pickle
import numpy as np
import pandas as pd
from sentence_transformers import SentenceTransformer
from sklearn.metrics.pairwise import cosine_similarity
import requests

BASE_DIR = Path(r"C:\Users\koreo\Downloads\pv-signal-ml")
SAR_DIR = BASE_DIR / "sar_reports"
EMB_DIR = BASE_DIR / "rag_embeds"
EMB_DIR.mkdir(exist_ok=True, parents=True)

EMB_MODEL_NAME = "all-MiniLM-L6-v2"


def build_signal_embeddings():
    """Embed all signals from full_signals_1M.csv once."""
    signals_path = SAR_DIR / "full_signals_1M.csv"
    df = pd.read_csv(signals_path)

    texts = [
        f"{row.DRUG} | {row.EVENT} | PRR={row.PRR:.2f} | CASES={row.CASES}"
        for _, row in df.iterrows()
    ]

    model = SentenceTransformer(EMB_MODEL_NAME)
    embeds = model.encode(texts)

    np.save(EMB_DIR / "signal_embeddings.npy", embeds)
    df.to_csv(EMB_DIR / "signals_index.csv", index=False)
    with open(EMB_DIR / "texts.pkl", "wb") as f:
        pickle.dump(texts, f)

    print(f"✅ Stored {len(df)} signal embeddings at {EMB_DIR}")


def enrich_evidence_with_similar_signals(
    drug: str, event: str, top_k: int = 3
):
    """Return related internal signals as evidence items."""
    signals = pd.read_csv(EMB_DIR / "signals_index.csv")
    embeds = np.load(EMB_DIR / "signal_embeddings.npy")
    with open(EMB_DIR / "texts.pkl", "rb") as f:
        texts = pickle.load(f)

    model = SentenceTransformer(EMB_MODEL_NAME)
    query_text = f"{drug} | {event}"
    q_emb = model.encode([query_text])

    sims = cosine_similarity(q_emb, embeds)[0]
    top_idx = sims.argsort()[::-1][:top_k]

    items = []
    for idx in top_idx:
        row = signals.iloc[idx]
        items.append(
            {
                "type": "Related signal",
                "source": "Internal disproportionality analysis",
                "citation": f"{row.DRUG} / {row.EVENT}",
                "summary": (
                    f"Related signal: {row.DRUG} – {row.EVENT} "
                    f"(PRR={row.PRR:.2f}, CASES={row.CASES})."
                ),
            }
        )
    return items


def fetch_pubmed_snippets(
    drug: str,
    event: str,
    start_date: str,
    end_date: str,
    max_results: int = 5,
):
    """
    Fetch live PubMed abstracts for drug+event within date range.

    Dates: 'YYYY/MM/DD'. Filters to titles containing the event term and
    supports simple drug synonyms for better recall.[web:105][web:108]
    """
    synonyms = {
        "InsulPen (Insulin)": ["insulin pen", "insulin", "insulin injection"],
        "CardioFlow (Atenolol)": ["atenolol", "beta blocker"],
        "PainAway (Oxycodone)": ["oxycodone"],
        "ChemoGold (Doxorubicin)": ["doxorubicin", "anthracycline"],
        "StatShield (Atorvastatin)": ["atorvastatin", "statin"],
        # extend for remaining products as needed
    }
    extra_terms = synonyms.get(drug, [])
    drug_terms = [drug] + extra_terms
    drug_query = " OR ".join([f'"{t}"[Title/Abstract]' for t in drug_terms])

    query = f'({drug_query}) AND "{event}"[Title/Abstract]'
    params = {
        "db": "pubmed",
        "term": query,
        "retmax": max_results,
        "retmode": "json",
        "mindate": start_date,
        "maxdate": end_date,
        "datetype": "pdat",
    }
    r = requests.get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
        params=params,
        timeout=30,
    )
    r.raise_for_status()
    ids = r.json().get("esearchresult", {}).get("idlist", [])
    if not ids:
        return []

    id_str = ",".join(ids)
    r2 = requests.get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
        params={"db": "pubmed", "id": id_str, "rettype": "abstract", "retmode": "text"},
        timeout=30,
    )
    r2.raise_for_status()
    blocks = r2.text.split("\n\n")

    event_lower = event.lower()
    items = []
    for pmid, block in zip(ids, blocks):
        lines = [ln.strip() for ln in block.splitlines() if ln.strip()]
        if not lines:
            continue
        title = lines[0][:300]
        # Require event term in title to reduce noise
        if event_lower not in title.lower():
            continue
        abstract_snip = " ".join(lines[1:])[:1000]
        items.append(
            {
                "type": "Literature",
                "source": "PubMed (live)",
                "citation": f"PMID {pmid} – {title}",
                "summary": abstract_snip,
            }
        )
    return items


if __name__ == "__main__":
    build_signal_embeddings()
    extra = enrich_evidence_with_similar_signals(
        "CardioFlow (Atenolol)", "Bradycardia"
    )
    print("Sample related evidence items:")
    for e in extra:
        print(" -", e["summary"])
