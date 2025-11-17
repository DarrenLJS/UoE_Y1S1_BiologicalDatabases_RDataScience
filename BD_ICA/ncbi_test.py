#!/usr/bin/python3

import requests
import time
import json

# Configuration
NCBI_API_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EMAIL = "s2906787@bioinfmsc9.bio.ed.ac.uk"  # Replace with your email
API_KEY = "37ee5dc16d4c1a46c163da440cc333dca209"  # Replace with your NCBI API key

def test_search_queries():
    """Test different search query formats to find what works"""
    print("=" * 60)
    print("NCBI Gene Search Diagnostic Test")
    print("=" * 60)

    # Different query formats to test
    test_queries = [
        ("BRCA1[Gene Name] AND human[Organism]", "Standard format with filters"),
        ("BRCA1 AND human", "Simple format without filters"),
        ("BRCA1", "Gene name only"),
        ("672", "Direct Gene ID"),
        ("TP53[Gene] AND Homo sapiens[Organism]", "Full species name"),
        ("BRCA1[sym] AND Homo sapiens[orgn]", "Short field codes"),
    ]

    url = f"{NCBI_API_BASE}/esearch.fcgi"

    for query, description in test_queries:
        print(f"\n[Testing] {description}")
        print(f"Query: {query}")

        params = {
            "db": "gene",
            "term": query,
            "retmax": 3,
            "retmode": "json",
            "email": EMAIL
        }
        if API_KEY:
            params["api_key"] = API_KEY

        try:
            response = requests.get(url, params=params, timeout=10)
            print(f"Status: {response.status_code}")

            if response.status_code == 200:
                data = response.json()

                # Print full response for debugging
                print(f"Full response: {json.dumps(data, indent=2)}")

                result = data.get("esearchresult", {})
                count = result.get("count", "0")
                id_list = result.get("idlist", [])

                if int(count) > 0:
                    print(f"✓ SUCCESS - Found {count} results")
                    print(f"  Gene IDs: {id_list}")
                    return True, id_list[0] if id_list else None
                else:
                    print(f"✗ No results (count: {count})")
                    # Check for error messages
                    if "errormessagelist" in data:
                        print(f"  Errors: {data['errormessagelist']}")
                    if "warninglist" in result:
                        print(f"  Warnings: {result['warninglist']}")
            else:
                print(f"✗ HTTP Error: {response.status_code}")
                print(f"Response: {response.text[:300]}")

        except Exception as e:
            print(f"✗ Exception: {e}")

        time.sleep(0.4)  # Rate limit

    return False, None

def test_direct_gene_id():
    """Test if we can fetch a gene by known ID"""
    print("\n" + "=" * 60)
    print("Testing Direct Gene ID Fetch")
    print("=" * 60)

    gene_id = "672"  # BRCA1 gene ID
    print(f"\nFetching gene ID {gene_id} (BRCA1)...")

    url = f"{NCBI_API_BASE}/esummary.fcgi"
    params = {
        "db": "gene",
        "id": gene_id,
        "retmode": "json",
        "email": EMAIL
    }
    if API_KEY:
        params["api_key"] = API_KEY

    try:
        response = requests.get(url, params=params, timeout=10)
        print(f"Status: {response.status_code}")

        if response.status_code == 200:
            data = response.json()
            print(f"Response preview: {json.dumps(data, indent=2)[:500]}...")

            result = data.get("result", {}).get(gene_id, {})
            if result:
                symbol = result.get("nomenclaturesymbol", "Unknown")
                organism = result.get("organism", {}).get("scientificname", "Unknown")
                print(f"✓ Gene found: {symbol} ({organism})")
                return True
            else:
                print(f"✗ No data in response")
        else:
            print(f"✗ HTTP Error: {response.status_code}")

    except Exception as e:
        print(f"✗ Exception: {e}")

    return False

def check_gene_database_status():
    """Check if gene database is accessible"""
    print("\n" + "=" * 60)
    print("Checking Gene Database Status")
    print("=" * 60)

    url = f"{NCBI_API_BASE}/einfo.fcgi"
    params = {
        "db": "gene",
        "retmode": "json",
        "email": EMAIL
    }
    if API_KEY:
        params["api_key"] = API_KEY

    try:
        response = requests.get(url, params=params, timeout=10)
        print(f"Status: {response.status_code}")

        if response.status_code == 200:
            data = response.json()
            db_info = data.get("einforesult", {}).get("dbinfo", [{}])[0]

            print(f"Database: {db_info.get('dbname', 'Unknown')}")
            print(f"Description: {db_info.get('description', 'Unknown')}")
            print(f"Record Count: {db_info.get('count', 'Unknown')}")
            print(f"Last Update: {db_info.get('lastupdate', 'Unknown')}")
            print(f"✓ Gene database is accessible")
            return True
        else:
            print(f"✗ Cannot access database info")

    except Exception as e:
        print(f"✗ Exception: {e}")

    return False

def run_full_diagnostic():
    """Run complete diagnostic"""
    print("\nStarting NCBI API diagnostic...\n")

    if API_KEY:
        print(f"✓ Using API key: {API_KEY[:8]}...\n")
    else:
        print("⚠ No API key provided\n")

    # Check database status
    check_gene_database_status()
    time.sleep(0.4)

    # Test direct fetch
    test_direct_gene_id()
    time.sleep(0.4)

    # Test various search queries
    success, gene_id = test_search_queries()

    print("\n" + "=" * 60)
    print("Diagnostic Complete")
    print("=" * 60)

    if success:
        print("\n✓ NCBI API is working correctly")
    else:
        print("\n✗ Search queries are not returning results")
        print("Possible issues:")
        print("  1. Gene database might be temporarily unavailable")
        print("  2. Search syntax may have changed")
        print("  3. Network/firewall blocking specific queries")
        print("  4. NCBI server issues")

# Run the diagnostic
if __name__ == "__main__":
    run_full_diagnostic()
