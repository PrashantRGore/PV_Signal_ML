"""
RAG Knowledge Base for Pharmacovigilance Literature
Ingests drug safety literature for evidence-based SAR generation
"""
import os
from pathlib import Path
from langchain_text_splitters import RecursiveCharacterTextSplitter
from langchain_community.vectorstores import Chroma
from langchain_community.embeddings import HuggingFaceEmbeddings
from langchain_core.documents import Document

class PVKnowledgeBase:
    def __init__(self, persist_dir="data/rag_knowledge_base"):
        self.persist_dir = persist_dir
        # Use medical domain embeddings
        self.embeddings = HuggingFaceEmbeddings(
            model_name="pritamdeka/S-PubMedBert-MS-MARCO",
            model_kwargs={'device': 'cpu'}
        )
        self.vectorstore = None
        
    def load_or_create(self):
        """Load existing knowledge base or create new"""
        if os.path.exists(self.persist_dir):
            try:
                self.vectorstore = Chroma(
                    persist_directory=self.persist_dir,
                    embedding_function=self.embeddings
                )
                return True
            except:
                return False
        return False
    
    def add_documents(self, documents, sources):
        """
        Add documents to knowledge base
        documents: List of text strings
        sources: List of source identifiers (PubMed IDs, WHO-ART codes, etc.)
        """
        text_splitter = RecursiveCharacterTextSplitter(
            chunk_size=500,
            chunk_overlap=50,
            separators=["\n\n", "\n", ". ", " "]
        )
        
        docs = []
        for text, source in zip(documents, sources):
            chunks = text_splitter.split_text(text)
            docs.extend([
                Document(page_content=chunk, metadata={"source": source})
                for chunk in chunks
            ])
        
        if self.vectorstore is None:
            self.vectorstore = Chroma.from_documents(
                documents=docs,
                embedding=self.embeddings,
                persist_directory=self.persist_dir
            )
        else:
            self.vectorstore.add_documents(docs)
        
        self.vectorstore.persist()
        return len(docs)
    
    def search(self, query, k=5):
        """Retrieve relevant documents"""
        if self.vectorstore is None:
            raise ValueError("Knowledge base not loaded")
        return self.vectorstore.similarity_search(query, k=k)
