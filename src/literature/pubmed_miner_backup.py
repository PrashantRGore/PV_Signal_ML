"""
Enhanced Literature Mining for Pharmacovigilance
Professional-grade PubMed search with synonyms and MeSH
"""
from Bio import Entrez, Medline
from datetime import datetime, timedelta
import time
from typing import List, Dict

class DrugSynonymExpander:
    DRUG_SYNONYMS = {
        'ibuprofen': ['ibuprofen', 'advil', 'motrin', 'brufen', 'nurofen'],
        'aspirin': ['aspirin', 'acetylsalicylic acid', 'ASA'],
        'meloxicam': ['meloxicam', 'mobic'],
    }
    
    @staticmethod
    def expand(drug_name: str) -> List[str]:
        drug_lower = drug_name.lower().strip()
        for key, synonyms in DrugSynonymExpander.DRUG_SYNONYMS.items():
            if drug_lower in [s.lower() for s in synonyms]:
                return synonyms
        return [drug_name]

class EnhancedPubMedMiner:
    def __init__(self, email='pharmacovigilance@example.com'):
        Entrez.email = email
        Entrez.tool = "PV_Signal_ML"
        self.drug_expander = DrugSynonymExpander()
    
    def search_drug_event(self, drug_name: str, event_name: str, max_results: int = 10, date_range_years: int = 10) -> Dict:
        try:
            drug_synonyms = self.drug_expander.expand(drug_name)
            queries = self._build_queries(drug_synonyms, event_name)
            
            all_papers = []
            for query in queries[:3]:  # Try top 3 strategies
                papers = self._execute_search(query, max_results, date_range_years)
                if papers:
                    all_papers.extend(papers)
                    if len(all_papers) >= max_results:
                        break
                time.sleep(0.4)
            
            unique_papers = self._deduplicate_papers(all_papers)
            ranked_papers = self._rank_by_relevance(unique_papers, drug_name, event_name)
            
            return {
                'found': len(ranked_papers),
                'papers': ranked_papers[:max_results],
                'queries_used': queries[:3]
            }
        except Exception as e:
            return {'found': 0, 'papers': [], 'error': str(e)}
    
    def _build_queries(self, drug_synonyms: List[str], event_term: str) -> List[str]:
        queries = []
        drug_query = ' OR '.join([f'"{d}"[Title/Abstract]' for d in drug_synonyms])
        event_query = f'"{event_term}"[Title/Abstract]'
        
        # Query 1: Comprehensive
        queries.append(f'({drug_query}) AND ({event_query}) AND ("adverse effects"[Subheading] OR "adverse event"[Title/Abstract])')
        
        # Query 2: Case reports
        queries.append(f'({drug_query}) AND ({event_query}) AND ("Case Reports"[Publication Type] OR "case report"[Title/Abstract])')
        
        # Query 3: Broad
        queries.append(f'({drug_synonyms[0]}[Title/Abstract]) AND ({event_term}[Title/Abstract])')
        
        # Query 4: Ultra-broad
        queries.append(f'({drug_synonyms[0]}) AND ({event_term})')
        
        return queries
    
    def _execute_search(self, query: str, max_results: int, date_range_years: int) -> List[Dict]:
        try:
            end_date = datetime.now()
            start_date = end_date - timedelta(days=365 * date_range_years)
            
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
                return []
            
            return self._fetch_paper_details(pmids)
        except Exception as e:
            print(f"Search error: {e}")
            return []
    
    def _fetch_paper_details(self, pmids: List[str]) -> List[Dict]:
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
                    
                    abstract = ''
                    if 'Abstract' in article:
                        abstract_parts = article['Abstract'].get('AbstractText', [])
                        abstract = ' '.join([str(part) for part in abstract_parts])
                    
                    pub_types = []
                    if 'PublicationTypeList' in article:
                        pub_types = [str(pt) for pt in article['PublicationTypeList']]
                    
                    paper = {
                        'pmid': str(record['MedlineCitation']['PMID']),
                        'title': article.get('ArticleTitle', 'No title'),
                        'authors': ', '.join(authors) if authors else 'Unknown',
                        'journal': article.get('Journal', {}).get('Title', 'Unknown'),
                        'year': year,
                        'abstract': abstract,
                        'publication_types': pub_types,
                        'url': f'https://pubmed.ncbi.nlm.nih.gov/{record["MedlineCitation"]["PMID"]}/'
                    }
                    papers.append(paper)
                except:
                    continue
            
            return papers
        except Exception as e:
            print(f"Fetch error: {e}")
            return []
    
    def _deduplicate_papers(self, papers: List[Dict]) -> List[Dict]:
        seen_pmids = set()
        unique_papers = []
        for paper in papers:
            if paper['pmid'] not in seen_pmids:
                seen_pmids.add(paper['pmid'])
                unique_papers.append(paper)
        return unique_papers
    
    def _rank_by_relevance(self, papers: List[Dict], drug_name: str, event_name: str) -> List[Dict]:
        drug_lower = drug_name.lower()
        event_lower = event_name.lower()
        
        for paper in papers:
            score = 0
            title_lower = paper['title'].lower()
            abstract_lower = paper['abstract'].lower()
            
            if drug_lower in title_lower:
                score += 10
            if event_lower in title_lower:
                score += 10
            if drug_lower in abstract_lower:
                score += 5
            if event_lower in abstract_lower:
                score += 5
            
            pub_types = [pt.lower() for pt in paper['publication_types']]
            if 'case reports' in pub_types:
                score += 8
            
            try:
                year = int(paper['year'])
                if year >= datetime.now().year - 5:
                    score += 5
            except:
                pass
            
            paper['relevance_score'] = score
        
        return sorted(papers, key=lambda x: x.get('relevance_score', 0), reverse=True)

class LiteratureEvidenceGenerator:
    def __init__(self):
        self.miner = EnhancedPubMedMiner()
    
    def generate_evidence_section(self, drug_name: str, event_name: str, max_papers: int = 5) -> Dict:
        lit_results = self.miner.search_drug_event(drug_name, event_name, max_results=max_papers)
        
        if lit_results['found'] == 0:
            evidence_text = f'''**Literature Review (PubMed):**

**Search Strategy:** Comprehensive multi-query search
- Drug synonyms expansion
- Case report focus
- Date range: Past 10 years

**Results:** No published literature found for {drug_name} and {event_name}.

**Clinical Context:** Absence of literature may indicate rare/novel association. FAERS data may detect signals before publication.

**Queries Attempted:** {len(lit_results.get('queries_used', []))} strategies'''
            
            return {'evidence_text': evidence_text, 'references': [], 'n_papers': 0}
        
        evidence_lines = [
            f'**Literature Review ({lit_results["found"]} publications):**',
            '',
            '**Key Publications:**',
            ''
        ]
        
        references = []
        for i, paper in enumerate(lit_results['papers'], 1):
            pub_type = ' [Case Report]' if 'Case Reports' in str(paper.get('publication_types', [])) else ''
            evidence_lines.append(
                f'{i}. {paper["authors"]} ({paper["year"]}). {paper["title"]}{pub_type}. '
                f'*{paper["journal"]}*. PMID: {paper["pmid"]}'
            )
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
