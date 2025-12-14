"""
PubMed Literature Miner
"""
from Bio import Entrez
from datetime import datetime, timedelta
import re

class PubMedMiner:
    def __init__(self, email='pharmacovigilance@example.com'):
        Entrez.email = email
    
    def search_drug_event(self, drug_name, event_name, max_results=5, date_range_years=10):
        try:
            drug_clean = self._clean_term(drug_name)
            event_clean = self._clean_term(event_name)
            
            end_date = datetime.now()
            start_date = end_date - timedelta(days=365 * date_range_years)
            
            query = f'({drug_clean}[Title/Abstract]) AND ({event_clean}[Title/Abstract]) AND ("adverse effects"[Subheading] OR "adverse event"[Title/Abstract])'
            
            handle = Entrez.esearch(
                db='pubmed',
                term=query,
                retmax=max_results,
                datetype='pdat',
                mindate=start_date.strftime('%Y/%m/%d'),
                maxdate=end_date.strftime('%Y/%m/%d'),
                sort='relevance'
            )
            
            record = Entrez.read(handle)
            handle.close()
            
            pmids = record['IdList']
            
            if not pmids:
                return {'found': 0, 'papers': [], 'query': query}
            
            papers = self._fetch_paper_details(pmids)
            
            return {
                'found': len(papers),
                'papers': papers,
                'query': query,
                'search_date': datetime.now().isoformat()
            }
        except Exception as e:
            return {'found': 0, 'papers': [], 'error': str(e)}
    
    def _fetch_paper_details(self, pmids):
        try:
            ids = ','.join(pmids)
            handle = Entrez.efetch(db='pubmed', id=ids, rettype='xml')
            records = Entrez.read(handle)
            handle.close()
            
            papers = []
            for record in records['PubmedArticle']:
                try:
                    article = record['MedlineCitation']['Article']
                    
                    authors = []
                    if 'AuthorList' in article:
                        for author in article['AuthorList'][:3]:
                            if 'LastName' in author:
                                authors.append(f"{author.get('LastName', '')} {author.get('Initials', '')}")
                    
                    year = 'N/A'
                    if 'Journal' in article and 'JournalIssue' in article['Journal']:
                        pub_date = article['Journal']['JournalIssue'].get('PubDate', {})
                        year = pub_date.get('Year', 'N/A')
                    
                    paper = {
                        'pmid': str(record['MedlineCitation']['PMID']),
                        'title': article.get('ArticleTitle', 'No title'),
                        'authors': ', '.join(authors) if authors else 'Unknown',
                        'journal': article.get('Journal', {}).get('Title', 'Unknown'),
                        'year': year,
                        'url': f'https://pubmed.ncbi.nlm.nih.gov/{record["MedlineCitation"]["PMID"]}/'
                    }
                    
                    papers.append(paper)
                except:
                    continue
            
            return papers
        except:
            return []
    
    def _clean_term(self, term):
        term = re.sub(r'[^\w\s-]', '', term)
        return term.strip()

class LiteratureEvidenceGenerator:
    def __init__(self):
        self.miner = PubMedMiner()
    
    def generate_evidence_section(self, drug_name, event_name, max_papers=5):
        lit_results = self.miner.search_drug_event(drug_name, event_name, max_results=max_papers)
        
        if lit_results['found'] == 0:
            return {
                'evidence_text': f'**Literature Review:** No published literature found in PubMed for {drug_name} and {event_name} in the past 10 years.',
                'references': [],
                'n_papers': 0
            }
        
        evidence_lines = [f'**Literature Review:** {lit_results["found"]} relevant publications found:', '']
        references = []
        
        for i, paper in enumerate(lit_results['papers'], 1):
            evidence_lines.append(f'{i}. {paper["authors"]} ({paper["year"]}). {paper["title"]}. *{paper["journal"]}*. PMID: {paper["pmid"]}')
            evidence_lines.append(f'   {paper["url"]}')
            evidence_lines.append('')
            
            references.append({
                'pmid': paper['pmid'],
                'citation': f'{paper["authors"]} ({paper["year"]}). {paper["title"]}. {paper["journal"]}.',
                'url': paper['url']
            })
        
        return {
            'evidence_text': '\n'.join(evidence_lines),
            'references': references,
            'n_papers': lit_results['found']
        }
