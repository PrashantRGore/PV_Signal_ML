import pandas as pd
from pathlib import Path
import json
from datetime import datetime
from langchain_community.embeddings import HuggingFaceEmbeddings
from langchain_ollama import OllamaLLM
from langchain_community.vectorstores import Chroma
from langchain_text_splitters import RecursiveCharacterTextSplitter
from langchain_core.documents import Document

class PVSignalRAGLangChain:
    def __init__(self):
        self.base_dir = Path.cwd()
        self.embeddings = HuggingFaceEmbeddings(model_name="all-MiniLM-L6-v2")
        self.llm = OllamaLLM(model="llama3.2:3b", temperature=0.1)
        self.vectorstore = self._build_vectorstore()
    
    def _load_context_docs(self):
        docs = []
        # SAR reports
        for report in (self.base_dir / 'sar_reports' / 'reports').glob('*.md'):
            with open(report, 'r', encoding='utf-8') as f:
                docs.append(Document(page_content=f.read()[:4000], metadata={"source": "SAR"}))
        # Regulatory
        regulatory = ["EMA PRR ≥ 2, Chi² ≥ 4 requires review", "FDA FAERS signals need clinical validation", "CIOMS VIII: Biological plausibility required"]
        for reg in regulatory: docs.append(Document(page_content=reg, metadata={"source": "REGULATORY"}))
        return docs
    
    def _build_vectorstore(self):
        docs = self._load_context_docs()
        splitter = RecursiveCharacterTextSplitter(chunk_size=1000, chunk_overlap=200)
        splits = splitter.split_documents(docs)
        return Chroma.from_documents(splits, self.embeddings, persist_directory="./chroma_db_pv")
    
    def generate_sar_report(self, signal_row):
        drug, event = signal_row['DRUG'], signal_row['EVENT']
        prr, chi2, cases = signal_row['PRR'], signal_row['CHISQ'], signal_row['CASES']
        query = f"{drug} {event} pharmacovigilance"
        docs = self.vectorstore.similarity_search(query, k=3)
        context = "\n".join([doc.page_content for doc in docs])
        
        response = self.llm.invoke(f"""CIOMS XIV SAR for {drug}-{event}:
PRR={prr:.1f}, Chi²={chi2:.1f}, Cases={cases}
Context: {context[:2000]}

Generate: Executive Summary + Regulatory Assessment + Recommendation""")
        
        return {
            'drug': drug, 'event': event, 'prr': prr, 'chi2': chi2, 'cases': cases,
            'sar_content': response, 'timestamp': datetime.now().isoformat()
        }
    
    def batch_generate_sars(self, top_n=20):
        signals = pd.read_csv('sar_reports/candidate_signals_2025Q1_stats_engine.csv').head(top_n)
        reports = []
        for _, signal in signals.iterrows():
            report = self.generate_sar_report(signal)
            Path('sar_reports').mkdir(exist_ok=True)
            # Sanitize drug and event names for safe filenames
            safe_drug = str(signal.DRUG).replace('/', '_').replace('\\', '_').replace(':', '_').replace('*', '_').replace('?', '_').replace('"', '_').replace('<', '_').replace('>', '_').replace('|', '_')
            safe_event = str(signal.EVENT).replace('/', '_').replace('\\', '_').replace(':', '_').replace('*', '_').replace('?', '_').replace('"', '_').replace('<', '_').replace('>', '_').replace('|', '_')
            report_file = Path('sar_reports') / f"{safe_drug}_{safe_event}_SAR.json"
            with open(report_file, 'w', encoding='utf-8') as f:
                json.dump(report, f, indent=2)
            reports.append(report)
        bundle = {'sars': reports, 'total': len(reports)}
        with open('sar_reports/sar_bundle_langchain.json', 'w', encoding='utf-8') as f:
            json.dump(bundle, f, indent=2)
        return bundle
