"""
PubMed Literature Miner with RAG
Fetches relevant publications for evidence generation
"""
from Bio import Entrez
import requests
import time
from datetime import datetime, timedelta
import re

class PubMedMiner:
    '''
    Mine PubMed for drug-event literature
    Uses NCBI E-utilities API
    '''
    
    def __init__(self, email='pharmacovigilance@example.com'):
        Entrez.email = email
        self.base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
    
    def search_drug_event(self, drug_name, event_name, max_results=10, date_range_years=10):
        '''
        Search PubMed for drug-event literature
        
        Parameters:
        - drug_name: Drug name
        - event_name: Adverse event term
        - max_results: Maximum papers to retrieve
        - date_range_years: Look back period
        '''
        try:
            # Clean terms
            drug_clean = self._clean_term(drug_name)
            event_clean = self._clean_term(event_name)
            
            # Calculate date range
            end_date = datetime.now()
            start_date = end_date - timedelta(days=365 * date_range_years)
            
            # Build query
            query = f'({drug_clean}[Title/Abstract]) AND ({event_clean}[Title/Abstract]) AND ("adverse effects"[Subheading] OR "adverse event"[Title/Abstract])'
            
            # Search PubMed
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
                return {
                    'found': 0,
                    'papers': [],
                    'query': query,
                    'message': 'No literature found'
                }
            
            # Fetch details
            papers = self._fetch_paper_details(pmids)
            
            return {
                'found': len(papers),
                'papers': papers,
                'query': query,
                'search_date': datetime.now().isoformat()
            }
        
        except Exception as e:
            return {
                'found': 0,
                'papers': [],
                'error': str(e)
            }
    
    def _fetch_paper_details(self, pmids):
        '''Fetch full paper details'''
        try:
            ids = ','.join(pmids)
            handle = Entrez.efetch(db='pubmed', id=ids, rettype='xml')
            records = Entrez.read(handle)
            handle.close()
            
            papers = []
            for record in records['PubmedArticle']:
                try:
                    article = record['MedlineCitation']['Article']
                    
                    # Extract authors
                    authors = []
                    if 'AuthorList' in article:
                        for author in article['AuthorList'][:3]:  # First 3 authors
                            if 'LastName' in author:
                                authors.append(f"{author.get('LastName', '')} {author.get('Initials', '')}")
                    
                    # Extract year
                    year = 'N/A'
                    if 'Journal' in article and 'JournalIssue' in article['Journal']:
                        pub_date = article['Journal']['JournalIssue'].get('PubDate', {})
                        year = pub_date.get('Year', 'N/A')
                    
                    paper = {
                        'pmid': str(record['MedlineCitation']['PMID']),
                        'title': article.get('ArticleTitle', 'No title'),
                        'abstract': article.get('Abstract', {}).get('AbstractText', ['No abstract'])[0] if 'Abstract' in article else 'No abstract',
                        'authors': ', '.join(authors) if authors else 'Unknown',
                        'journal': article.get('Journal', {}).get('Title', 'Unknown'),
                        'year': year,
                        'url': f'https://pubmed.ncbi.nlm.nih.gov/{record["MedlineCitation"]["PMID"]}/'
                    }
                    
                    papers.append(paper)
                
                except Exception as e:
                    continue
            
            return papers
        
        except Exception as e:
            return []
    
    def _clean_term(self, term):
        '''Clean search term'''
        # Remove special characters
        term = re.sub(r'[^\w\s-]', '', term)
        return term.strip()

class LiteratureEvidenceGenerator:
    '''Generate evidence section for SAR from literature'''
    
    def __init__(self):
        self.miner = PubMedMiner()
    
    def generate_evidence_section(self, drug_name, event_name, max_papers=5):
        '''Generate literature evidence section'''
        
        # Search literature
        lit_results = self.miner.search_drug_event(drug_name, event_name, max_results=max_papers)
        
        if lit_results['found'] == 0:
            return {
                'evidence_text': f'**Literature Review:** No published literature found in PubMed for the association between {drug_name} and {event_name} in the past 10 years. This represents a novel signal requiring further investigation.',
                'references': [],
                'n_papers': 0
            }
        
        # Build evidence text
        evidence_lines = [
            f'**Literature Review:** {lit_results["found"]} relevant publications identified in PubMed:',
            ''
        ]
        
        references = []
        for i, paper in enumerate(lit_results['papers'], 1):
            # Add citation
            evidence_lines.append(f'{i}. {paper["authors"]} ({paper["year"]}). {paper["title"]}. *{paper["journal"]}*. PMID: {paper["pmid"]}')
            evidence_lines.append(f'   URL: {paper["url"]}')
            evidence_lines.append('')
            
            references.append({
                'pmid': paper['pmid'],
                'citation': f'{paper["authors"]} ({paper["year"]}). {paper["title"]}. {paper["journal"]}.',
                'url': paper['url']
            })
        
        evidence_text = '\n'.join(evidence_lines)
        
        return {
            'evidence_text': evidence_text,
            'references': references,
            'n_papers': lit_results['found'],
            'search_query': lit_results.get('query', ''),
            'search_date': lit_results.get('search_date', '')
        }
