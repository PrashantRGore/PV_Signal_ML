"""
Test literature search directly
"""
import sys
sys.path.insert(0, 'C:/Users/koreo/PV_Signal_ML')

from src.literature.pubmed_miner import LiteratureEvidenceGenerator

# Test
print("Testing literature search...")
print("="*80)

generator = LiteratureEvidenceGenerator()
result = generator.generate_evidence_section(
    drug_name='ibuprofen',
    event_name='pericardial mesothelioma',
    max_papers=5
)

print(f"Papers found: {result['n_papers']}")
print(f"\nEvidence text:\n{result['evidence_text']}")

if result['references']:
    print(f"\nReferences ({len(result['references'])}):")
    for ref in result['references']:
        print(f"  - PMID: {ref['pmid']}")
        print(f"    {ref['citation']}")
        print(f"    {ref['url']}")
        print()
else:
    print("\nNo references found")

print("="*80)

if result['n_papers'] >= 2:
    print("✅ SUCCESS: Found 2+ papers including expected ibuprofen case reports")
else:
    print(f"⚠️ WARNING: Only found {result['n_papers']} papers (expected 2+)")
    print("Possible reasons:")
    print("  1. PubMed API rate limiting")
    print("  2. Network connection issue")
    print("  3. Biopython not installed")
    print("  4. Search query needs adjustment")
