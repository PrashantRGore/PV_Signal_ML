"""
MINIMAL Enhanced Literature Search - TESTED
Only fixes PubMed search to find those 2 ibuprofen papers
"""
from Bio import Entrez
from datetime import datetime, timedelta
import time
from typing import List, Dict

class EnhancedPubMedMiner:
    def __init__(self, email='pharmacovigilance@example.com'):
        Entrez.email = email
        Entrez.tool = "PV_Signal_ML"
    
    def search_drug_event(self, drug_name: str, event_name: str, max_results: int = 10, date_range_years: int = 10) -> Dict:
        '''Enhanced search with multiple strategies'''
        try:
            queries = self._build_queries(drug_name, event_name)
            
            all_papers = []
            for query in queries:
                papers = self._execute_search(query, max_results, date_range_years)
                if papers:
                    all_papers.extend(papers)
                    if len(all_papers) >= max_results:
                        break
                time.sleep(0.4)  # PubMed rate limiting
            
            unique_papers = self._deduplicate_papers(all_papers)
            
            return {
                'found': len(unique_papers),
                'papers': unique_papers[:max_results]
            }
        except Exception as e:
            print(f"Literature search error: {e}")
            return {'found': 0, 'papers': []}
    
    def _build_queries(self, drug_name: str, event_name: str) -> List[str]:
        '''Build multiple search strategies'''
        queries = []
        
        # Query 1: Title/Abstract search
        queries.append(f'("{drug_name}"[Title/Abstract]) AND ("{event_name}"[Title/Abstract])')
        
        # Query 2: Case reports specifically (critical for rare events)
        queries.append(f'("{drug_name}"[Title/Abstract]) AND ("{event_name}"[Title/Abstract]) AND ("Case Reports"[Publication Type])')
        
        # Query 3: Very broad (catches everything)
        queries.append(f'({drug_name}) AND ({event_name})')
        
        return queries
    
    def _execute_search(self, query: str, max_results: int, date_range_years: int) -> List[Dict]:
        '''Execute single search'''
        try:
            end_date = datetime.now()
            start_date = end_date - timedelta(days=365 * date_range_years)
            
            # Search
            handle = Entrez.esearch(
                db='pubmed',
                term=query,
                retmax=max_results,
                datetype='pdat',
                mindate=start_date.strftime('%Y/%m/%d'),
                maxdate=end_date.strftime('%Y/%m/%d')
            )
            
            record = Entrez.read(handle)
            handle.close()
            
            pmids = record.get('IdList', [])
            if not pmids:
                return []
            
            return self._fetch_details(pmids)
            
        except Exception as e:
            print(f"Search error: {e}")
            return []
    
    def _fetch_details(self, pmids: List[str]) -> List[Dict]:
        '''Fetch paper details'''
        try:
            ids = ','.join(pmids)
            handle = Entrez.efetch(db='pubmed', id=ids, rettype='xml')
            records = Entrez.read(handle)
            handle.close()
            
            papers = []
            for record in records.get('PubmedArticle', []):
                try:
                    article = record['MedlineCitation']['Article']
                    pmid = str(record['MedlineCitation']['PMID'])
                    
                    # Authors
                    authors = []
                    if 'AuthorList' in article:
                        for author in article['AuthorList'][:3]:
                            if 'LastName' in author:
                                authors.append(f"{author.get('LastName', '')} {author.get('Initials', '')}")
                    
                    # Year
                    year = 'N/A'
                    if 'Journal' in article and 'JournalIssue' in article['Journal']:
                        pub_date = article['Journal']['JournalIssue'].get('PubDate', {})
                        year = pub_date.get('Year', 'N/A')
                    
                    # Abstract
                    abstract = ''
                    if 'Abstract' in article:
                        abstract_parts = article['Abstract'].get('AbstractText', [])
                        abstract = ' '.join([str(part) for part in abstract_parts])
                    
                    paper = {
                        'pmid': pmid,
                        'title': str(article.get('ArticleTitle', 'No title')),
                        'authors': ', '.join(authors) if authors else 'Unknown',
                        'journal': str(article.get('Journal', {}).get('Title', 'Unknown')),
                        'year': str(year),
                        'abstract': abstract,
                        'url': f'https://pubmed.ncbi.nlm.nih.gov/{pmid}/'
                    }
                    
                    papers.append(paper)
                    
                except Exception as e:
                    print(f"Parse error: {e}")
                    continue
            
            return papers
            
        except Exception as e:
            print(f"Fetch error: {e}")
            return []
    
    def _deduplicate_papers(self, papers: List[Dict]) -> List[Dict]:
        '''Remove duplicates'''
        seen = set()
        unique = []
        for paper in papers:
            if paper['pmid'] not in seen:
                seen.add(paper['pmid'])
                unique.append(paper)
        return unique

class LiteratureEvidenceGenerator:
    def __init__(self):
        self.miner = EnhancedPubMedMiner()
    
    def generate_evidence_section(self, drug_name: str, event_name: str, max_papers: int = 5) -> Dict:
        '''Generate literature evidence with citations'''
        
        lit_results = self.miner.search_drug_event(drug_name, event_name, max_results=max_papers)
        
        if lit_results['found'] == 0:
            evidence_text = f'''**Literature Review (PubMed):**

No published literature found in PubMed for {drug_name} and {event_name} in the past 10 years.

**Note:** Absence of literature does not exclude causality. FAERS data may detect signals before publication.'''
            
            return {
                'evidence_text': evidence_text,
                'references': [],
                'n_papers': 0
            }
        
        # Format with full citations
        evidence_lines = [
            f'**Literature Review ({lit_results["found"]} publication(s) identified):**',
            '',
            '**Publications:**',
            ''
        ]
        
        references = []
        
        for i, paper in enumerate(lit_results['papers'], 1):
            # Full citation format
            citation_line = f"{i}. **{paper['authors']}** ({paper['year']}). {paper['title']}. *{paper['journal']}*."
            evidence_lines.append(citation_line)
            evidence_lines.append(f"   PMID: {paper['pmid']}")
            evidence_lines.append(f"   URL: {paper['url']}")
            evidence_lines.append('')
            
            references.append({
                'pmid': paper['pmid'],
                'citation': f"{paper['authors']} ({paper['year']}). {paper['title']}. {paper['journal']}.",
                'url': paper['url']
            })
        
        return {
            'evidence_text': '\n'.join(evidence_lines),
            'references': references,
            'n_papers': lit_results['found']
        }
