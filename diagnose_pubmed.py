"""
Diagnose PubMed search issues
"""
from Bio import Entrez
from datetime import datetime, timedelta

Entrez.email = "test@example.com"
Entrez.tool = "PV_Signal_ML_Test"

print("="*80)
print("PUBMED SEARCH DIAGNOSTICS")
print("="*80)

# Test 1: Simple connection test
print("\n1️⃣ Testing PubMed connection...")
try:
    handle = Entrez.esearch(db='pubmed', term='aspirin', retmax=1)
    record = Entrez.read(handle)
    handle.close()
    print(f"   ✅ Connected! Found {record['Count']} total aspirin papers")
except Exception as e:
    print(f"   ❌ Connection failed: {e}")
    print("   → Check internet connection")
    exit(1)

# Test 2: Search for ibuprofen (should find many)
print("\n2️⃣ Testing ibuprofen search...")
try:
    handle = Entrez.esearch(db='pubmed', term='ibuprofen', retmax=1)
    record = Entrez.read(handle)
    handle.close()
    print(f"   ✅ Found {record['Count']} ibuprofen papers")
except Exception as e:
    print(f"   ❌ Search failed: {e}")

# Test 3: Search for pericardial mesothelioma
print("\n3️⃣ Testing 'pericardial mesothelioma' search...")
try:
    handle = Entrez.esearch(db='pubmed', term='pericardial mesothelioma', retmax=10)
    record = Entrez.read(handle)
    handle.close()
    count = record['Count']
    pmids = record.get('IdList', [])
    print(f"   ✅ Found {count} pericardial mesothelioma papers")
    if pmids:
        print(f"   📋 Sample PMIDs: {pmids[:3]}")
except Exception as e:
    print(f"   ❌ Search failed: {e}")

# Test 4: Combined search (exact query from our code)
print("\n4️⃣ Testing combined search (ibuprofen + pericardial mesothelioma)...")
try:
    query = '("ibuprofen"[Title/Abstract]) AND ("pericardial mesothelioma"[Title/Abstract])'
    print(f"   Query: {query}")
    
    handle = Entrez.esearch(db='pubmed', term=query, retmax=10)
    record = Entrez.read(handle)
    handle.close()
    count = record['Count']
    pmids = record.get('IdList', [])
    print(f"   Result: {count} papers")
    if pmids:
        print(f"   📋 PMIDs: {pmids}")
    else:
        print(f"   ⚠️ No papers found with this exact query")
except Exception as e:
    print(f"   ❌ Search failed: {e}")

# Test 5: Try broader search
print("\n5️⃣ Testing broader search (any field)...")
try:
    query = '(ibuprofen) AND (pericardial mesothelioma)'
    print(f"   Query: {query}")
    
    handle = Entrez.esearch(db='pubmed', term=query, retmax=10)
    record = Entrez.read(handle)
    handle.close()
    count = record['Count']
    pmids = record.get('IdList', [])
    print(f"   Result: {count} papers")
    if pmids:
        print(f"   📋 PMIDs: {pmids}")
        
        # Fetch details of first paper
        if pmids:
            print(f"\n   📄 Fetching details of first paper...")
            handle2 = Entrez.efetch(db='pubmed', id=pmids[0], rettype='xml')
            records = Entrez.read(handle2)
            handle2.close()
            
            article = records['PubmedArticle'][0]['MedlineCitation']['Article']
            title = article.get('ArticleTitle', 'No title')
            print(f"   Title: {title}")
except Exception as e:
    print(f"   ❌ Search failed: {e}")

# Test 6: Try with date range
print("\n6️⃣ Testing with 10-year date range...")
try:
    end_date = datetime.now()
    start_date = end_date - timedelta(days=365 * 10)
    
    query = '(ibuprofen) AND (pericardial mesothelioma)'
    print(f"   Query: {query}")
    print(f"   Date: {start_date.strftime('%Y/%m/%d')} to {end_date.strftime('%Y/%m/%d')}")
    
    handle = Entrez.esearch(
        db='pubmed',
        term=query,
        retmax=10,
        datetype='pdat',
        mindate=start_date.strftime('%Y/%m/%d'),
        maxdate=end_date.strftime('%Y/%m/%d')
    )
    record = Entrez.read(handle)
    handle.close()
    count = record['Count']
    pmids = record.get('IdList', [])
    print(f"   Result: {count} papers in past 10 years")
    if pmids:
        print(f"   📋 PMIDs: {pmids}")
except Exception as e:
    print(f"   ❌ Search failed: {e}")

# Test 7: Try even broader (just mesothelioma + ibuprofen)
print("\n7️⃣ Testing very broad search (mesothelioma + ibuprofen)...")
try:
    query = '(ibuprofen) AND (mesothelioma)'
    print(f"   Query: {query}")
    
    handle = Entrez.esearch(db='pubmed', term=query, retmax=10)
    record = Entrez.read(handle)
    handle.close()
    count = record['Count']
    pmids = record.get('IdList', [])
    print(f"   Result: {count} papers")
    if pmids:
        print(f"   📋 PMIDs: {pmids[:5]}")
except Exception as e:
    print(f"   ❌ Search failed: {e}")

print("\n" + "="*80)
print("DIAGNOSIS COMPLETE")
print("="*80)

print("\n💡 RECOMMENDATIONS:")
print("   - If Tests 1-3 pass but 4-7 fail: The drug-event combo is very rare")
print("   - If Test 7 finds papers: Use broader search terms")
print("   - If all tests fail: Network/API issue")
print("\n📝 Copy all output above and share for analysis")
