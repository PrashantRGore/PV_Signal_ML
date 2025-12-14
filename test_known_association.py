"""
Test literature search with known association
"""
from src.literature.pubmed_miner import LiteratureEvidenceGenerator

print("Testing with KNOWN association: Aspirin + Gastrointestinal bleeding")
print("="*80)

generator = LiteratureEvidenceGenerator()
result = generator.generate_evidence_section(
    drug_name='aspirin',
    event_name='gastrointestinal bleeding',
    max_papers=5
)

print(f"Papers found: {result['n_papers']}")
print(f"\nEvidence:\n{result['evidence_text']}")

if result['n_papers'] > 0:
    print("\n✅ SUCCESS: Literature citations found!")
    print("This proves the system works for known associations.")
else:
    print("\n⚠️ Unexpected: No papers found even for known association")
