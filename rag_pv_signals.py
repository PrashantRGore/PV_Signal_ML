"""
PV Signal Detection + RAG (Audit-style)
Uses ranked_signals_2020Q2.csv with evidence columns and answers in a regulator-friendly format.
"""

import pandas as pd
from pathlib import Path
import numpy as np
from sentence_transformers import SentenceTransformer
from sklearn.metrics.pairwise import cosine_similarity
import ollama

BASE_DIR = Path(r"C:\Users\koreo\Downloads\pv-signal-ml")
RESULTS_DIR = BASE_DIR / "results"

print(f"Working in: {BASE_DIR}")

class PVSignalRAG:
    def __init__(self):
        print("[INFO] Loading signal data...")
        ranked_path = RESULTS_DIR / "ranked_signals_2020Q2.csv"
        if not ranked_path.exists():
            print(f"[ERROR] Missing: {ranked_path}")
            print("Run: python pv_signal_ml_pipeline.py")
            exit(1)

        self.signals_df = pd.read_csv(ranked_path)
        self.model = SentenceTransformer("all-MiniLM-L6-v2")

        self._build_embeddings()
        print(f"✅ RAG ready: {len(self.signals_df):,} signals")

    def _build_embeddings(self):
        self.signals_df["text"] = (
            self.signals_df["DRUGNAME_NORM"].astype(str)
            + " "
            + self.signals_df["pt"].astype(str)
            + " "
            + self.signals_df["signal_period"].astype(str)
        )
        self.embeddings = self.model.encode(self.signals_df["text"].tolist())

    def retrieve(self, query, top_k=5):
        query_emb = self.model.encode([query])
        similarities = cosine_similarity(query_emb, self.embeddings)[0]
        top_idx = np.argsort(similarities)[::-1][:top_k]
        cols = [
            "DRUGNAME_NORM",
            "pt",
            "signal_period",
            "score",
            "label_prr",
            "A_6m",
            "N_6m",
            "prr_6m",
            "A_12m",
            "N_12m",
            "model_version",
        ]
        return self.signals_df.iloc[top_idx][cols].copy()

    def generate(self, query, context):
        context_text = context.to_string(index=False)
        prompt = f"""You are a pharmacovigilance ML system being reviewed by regulators (FDA/EMA).

Query: "{query}"

The table below shows drug–event pairs with model outputs:
{context_text}

Explain in an audit-ready way:
1. Identify which pairs look like stronger signals based on score and prr_6m.
2. For each highlighted pair, report A_6m, N_6m, prr_6m, label_prr, signal_period, and model_version.
3. Briefly describe how the model works: XGBoost on FAERS data, features include A_6m, N_6m, prr_6m; labels: label_prr = 1 if prr_6m >= 0.1 and A_6m >= 3.
4. Mention global validation performance (AP ≈ 0.022, AUC ≈ 0.9999) and that the model supports human review, not replaces it.
5. Keep the answer concise: 1–2 short paragraphs and a bullet list of key evidence.
"""
        try:
            response = ollama.generate(model="llama3.2:3b", prompt=prompt)
            return response["response"]
        except Exception as e:
            return f"LLM error: {e}\nTop signals are shown in the table above."

def main():
    rag = PVSignalRAG()

    print("\n🚀 RAG READY! Try queries like:")
    print("   PROACTIV skin")
    print("   OPSUMIT cardiac")
    print("   ACETAMINOPHEN blood pressure")
    print("Type 'quit' to exit.")
    print("=" * 60)

    while True:
        query = input("\n🔍 Query: ").strip()
        if query.lower() in ["quit", "q", "exit"]:
            break

        print("\n📡 Top matches:")
        context = rag.retrieve(query, top_k=5)
        print(context.round(4).to_string(index=False))

        print("\n🤖 Audit-style analysis:")
        analysis = rag.generate(query, context)
        print(analysis)
        print("-" * 60)

if __name__ == "__main__":
    main()
