#!/usr/bin/python3

ens_id = [
    "ENSMUSG00000075394", 
    "ENSMUSG00000001655", 
    "ENSMUSG00000001657", 
    "ENSMUSG00000001661", 
    "ENSMUSG00000076010", 
    "ENSMUSG00000086903", 
    "ENSMUSG00000102297", 
    "ENSMUSG00000116066", 
    "ENSMUSG00000036061", 
    "ENSMUSG00000009575" 
]

import mysql.connector
import json
import pandas as pd
import requests
from Bio import Entrez
from pybiomart import Server
from pybiomart import Dataset

def main():
    print("\nAccessing Ensembl...\n")
    query_df = access_ensembl_rest(ens_id)
    print("\nDownloaded Ensembl results into SQL!\n")

    print("\nAccessing NCBI...\n")
    access_ncbi_entrezpy(query_df)
    print("\nDownloaded NCBI results into SQL!\n")
    
    print("\nAccessing InterPro...\n")
    access_interpro_api(query_df)
    print("\nDownloaded InterPro results into SQL!\n")
    
    print("\nMapping Ensembl IDs with BioMart...\n")
    xrefs_df = map_ensembl_pybiomart(query_df)
    print("\nMapped Ensembl IDs to GO IDs!\n")

    print("\nAccessing Gene Ontology (GO)...\n")
    access_go_api(xrefs_df)
    print("\nDownloaded GO results into SQL!\n")

    print("\nIntegrating SQL tables into one...\n")
    integration_df = integrate_sql()
    print("\nFinal integration completed!\n")
    print(integration_df.head())
    print(integration_df.shape)


def integrate_sql():
    connection = mysql.connector.connect(
        host = "bioinfmsc9.bio.ed.ac.uk", 
        user = "s2906787", 
        password = "r9UqKvEH", 
        database = "s2906787"
    )
    cursor = connection.cursor()
    
    cursor.execute("DROP TABLE IF EXISTS go_interpro")
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS go_interpro AS 
        SELECT 
        COALESCE(g.ensembl_id, i.ensembl_id) AS ensembl_id, 
        i.uniprot_id, 
        i.interpro_id, 
        i.interpro_name, 
        i.interpro_type, 
        g.go_id, 
        g.go_name, 
        g.go_aspect, 
        g.go_definition 
        FROM go_table g 
        LEFT JOIN interpro_table i 
        ON g.go_id = i.go_id 
        AND g.ensembl_id = i.ensembl_id 
        UNION
        SELECT 
        COALESCE(g.ensembl_id, i.ensembl_id) AS ensembl_id, 
        i.uniprot_id, 
        i.interpro_id, 
        i.interpro_name, 
        i.interpro_type, 
        g.go_id, 
        g.go_name, 
        g.go_aspect, 
        g.go_definition 
        FROM go_table g 
        RIGHT JOIN interpro_table i 
        ON g.go_id = i.go_id 
        AND g.ensembl_id = i.ensembl_id
    """)
    
    cursor.execute("DROP TABLE IF EXISTS ncbi_ensembl")
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS ncbi_ensembl AS 
        SELECT 
        COALESCE(n.ensembl_id, e.ensembl_id) AS ensembl_id, 
        n.gene_symbol, 
        n.organism, 
        n.taxon_id, 
        e.description, 
        n.chromosome, 
        e.strand, 
        n.start, 
        n.end, 
        n.entrez_id, 
        n.weight, 
        n.summary 
        FROM ncbi_gene_table n 
        LEFT JOIN ensembl_table e 
        ON n.ensembl_id = e.ensembl_id 
        UNION 
        SELECT 
        COALESCE(n.ensembl_id, e.ensembl_id) AS ensembl_id, 
        n.gene_symbol, 
        n.organism, 
        n.taxon_id, 
        e.description, 
        n.chromosome, 
        e.strand, 
        n.start, 
        n.end, 
        n.entrez_id, 
        n.weight, 
        n.summary 
        FROM ncbi_gene_table n 
        RIGHT JOIN ensembl_table e 
        ON n.ensembl_id = e.ensembl_id
    """)
    
    cursor.execute("DROP TABLE IF EXISTS ncbi_ensembl_pubmed")
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS ncbi_ensembl_pubmed AS 
        SELECT 
        ne.*, 
        p.pubmed_id, 
        p.doi, 
        p.publication_date, 
        p.author_list 
        FROM ncbi_ensembl ne 
        LEFT JOIN pubmed_table p 
        ON ne.ensembl_id = p.ensembl_id 
        UNION 
        SELECT 
        ne.*, 
        p.pubmed_id, 
        p.doi, 
        p.publication_date, 
        p.author_list 
        FROM ncbi_ensembl ne 
        RIGHT JOIN pubmed_table p 
        ON ne.ensembl_id = p.ensembl_id
    """)
    
    cursor.execute("DROP TABLE IF EXISTS final_integration")
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS final_integration AS 
        SELECT 
        COALESCE(gi.ensembl_id, nep.ensembl_id) AS ensembl_id, 
        nep.gene_symbol, 
        nep.organism, 
        nep.taxon_id, 
        nep.description, 
        nep.chromosome, 
        nep.strand, 
        nep.start, 
        nep.end, 
        nep.entrez_id, 
        nep.weight, 
        nep.summary, 
        nep.pubmed_id, 
        nep.doi, 
        nep.publication_date, 
        nep.author_list, 
        gi.uniprot_id, 
        gi.interpro_id, 
        gi.interpro_name, 
        gi.interpro_type, 
        gi.go_id, 
        gi.go_name, 
        gi.go_aspect, 
        gi.go_definition 
        FROM go_interpro gi 
        LEFT JOIN ncbi_ensembl_pubmed nep 
        ON gi.ensembl_id = nep.ensembl_id 
        UNION 
        SELECT 
        COALESCE(gi.ensembl_id, nep.ensembl_id) AS ensembl_id, 
        nep.gene_symbol, 
        nep.organism, 
        nep.taxon_id, 
        nep.description, 
        nep.chromosome, 
        nep.strand, 
        nep.start, 
        nep.end, 
        nep.entrez_id, 
        nep.weight, 
        nep.summary, 
        nep.pubmed_id, 
        nep.doi, 
        nep.publication_date, 
        nep.author_list, 
        gi.uniprot_id, 
        gi.interpro_id, 
        gi.interpro_name, 
        gi.interpro_type, 
        gi.go_id, 
        gi.go_name, 
        gi.go_aspect, 
        gi.go_definition 
        FROM go_interpro gi 
        RIGHT JOIN ncbi_ensembl_pubmed nep 
        ON gi.ensembl_id = nep.ensembl_id
    """)
    
    cursor.execute("DROP TABLE IF EXISTS go_interpro")
    cursor.execute("DROP TABLE IF EXISTS ncbi_ensembl")
    cursor.execute("DROP TABLE IF EXISTS ncbi_ensembl_pubmed")
    connection.commit()
    
    df_final = pd.read_sql("SELECT * FROM final_integration", connection)
    
    cursor.close()
    connection.close()
    return df_final

# access_ensembl_rest() takes in an arg: list of Ensembl IDs.
# Function: Query Ensembl using REST API, by looking up the provided list of Ensembl IDs.
# CREATE TABLE ensembl_table in SQL.
# If response is ok and there is at least a result, extract and insert relevant fields into ensembl_table.
# Print rows of ensembl_table to verify changes.
# Output: A temporary pandas DataFrame with cols "ensembl_id", "gene_symbol", "organism" for further queries.
def access_ensembl_rest(ens_id_list):
    connection = mysql.connector.connect(
        host = "bioinfmsc9.bio.ed.ac.uk", 
        user = "s2906787", 
        password = "r9UqKvEH", 
        database = "s2906787"
    )
    cursor = connection.cursor()
    
    cursor.execute("DROP TABLE IF EXISTS ensembl_table")
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS ensembl_table (
            ensembl_id VARCHAR(500), 
            gene_symbol VARCHAR(500), 
            organism VARCHAR(500), 
            description VARCHAR(500), 
            chromosome INTEGER(100), 
            strand INTEGER(100), 
            start INTEGER(100), 
            end INTEGER(100)
        )
    """)

    insert_ens = """
        INSERT INTO ensembl_table (
            ensembl_id, 
            gene_symbol, 
            organism, 
            description, 
            chromosome, 
            strand, 
            start, 
            end
        ) 
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
    """

    base_url = "https://rest.ensembl.org/"
    headers = {"Accept" : "application/json"}
    
    search_terms = []
    for ens_id in ens_id_list:
        print(f"\nSearching Ensembl for {ens_id}...")
        response_lookup = requests.get(f"{base_url}lookup/id/{ens_id}", headers = headers)
        
        if response_lookup.ok:
            data_ens = response_lookup.json()
            
            insert_data = (
                data_ens.get("id", None), 
                data_ens.get("display_name", None), 
                data_ens.get("species", None), 
                data_ens.get("description", None), 
                data_ens.get("seq_region_name", None), 
                data_ens.get("strand", None), 
                data_ens.get("start", None), 
                data_ens.get("end", None)
            )
            
            cursor.execute(insert_ens, insert_data)

            search_data = {
                "ensembl_id" : data_ens.get("id"), 
                "gene_symbol" : data_ens.get("display_name"), 
                "organism" : data_ens.get("species")
            }
            search_terms.append(search_data)
            print(f"Ensembl data successfully retrieved for {ens_id}!\n")

        else:
            print(f"Failed to retrieve Ensembl data for {ens_id}. Error: {response_lookup.status_code}\n")
            
            insert_data = (ens_id, None, None, None, None, None, None, None)
            cursor.execute(insert_ens, insert_data)
            connection.commit

            search_terms.append({
                "ensembl_id" : ens_id, 
                "gene_symbol" : None, 
                "organism" : None
            })

    connection.commit()
    cursor = connection.cursor(dictionary = True)
    cursor.execute("SELECT * FROM ensembl_table")
    for row in cursor:
        print(row)
    cursor.close()
    connection.close()

    search_df = pd.DataFrame(search_terms)

    def format_species_name(x):
        if pd.isna(x):
            return x
        x_list = x.split("_")
        x_list[0] = x_list[0].title()
        return " ".join(x_list)

    search_df["organism"] = search_df.apply(lambda x: format_species_name(x["organism"]), axis = 1)
    return search_df

# access_ncbi_entrezpy() takes in an arg: pandas DataFrame containing search terms.
# Function: Query NCBI gene and PubMed databases using entrezpy E-utilities.
# Entrez.esearch to get Entrez IDs for each search term.
# Entrez.elink to get PubMed IDs for each Entrez ID.
# Entrez.esummary to get summarised NCBI gene and PubMed results for each Entrez ID, 
# queried by entrez_id and pubmed_id respectively.
# CREATE TABLE ncbi_gene_table and pubmed_table in SQL.
# If response is ok and there is at least a result, extract and insert relevant fields into both tables.
# Print rows of ncbi_gene_table and pubmed_table to verify changes.
def access_ncbi_entrezpy(search_df):
    connection = mysql.connector.connect(
        host = "bioinfmsc9.bio.ed.ac.uk", 
        user = "s2906787", 
        password = "r9UqKvEH", 
        database = "s2906787"
    )
    cursor = connection.cursor()
    
    cursor.execute("DROP TABLE IF EXISTS ncbi_gene_table")
    cursor.execute("DROP TABLE IF EXISTS pubmed_table")
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS ncbi_gene_table (
            ensembl_id VARCHAR(500), 
            gene_symbol VARCHAR(500), 
            organism VARCHAR(500), 
            taxon_id VARCHAR(500), 
            chromosome INTEGER(100), 
            start INTEGER(100), 
            end INTEGER(100), 
            entrez_id VARCHAR(500), 
            weight INTEGER(100), 
            summary VARCHAR(3000)
        )
    """)
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS pubmed_table (
            ensembl_id VARCHAR(500), 
            pubmed_id JSON, 
            doi JSON, 
            publication_date JSON, 
            author_list JSON
        )
    """)

    insert_ncbi = """
        INSERT INTO ncbi_gene_table (
            ensembl_id, 
            gene_symbol, 
            organism, 
            taxon_id, 
            chromosome, 
            start, 
            end, 
            entrez_id, 
            weight, 
            summary
        ) 
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
    """
    insert_pubmed = """
        INSERT INTO pubmed_table (
            ensembl_id, 
            pubmed_id, 
            doi, 
            publication_date, 
            author_list
        ) 
        VALUES (%s, CAST(%s AS JSON), CAST(%s AS JSON), CAST(%s AS JSON), CAST(%s AS JSON))
    """
    
    Entrez.email = "s2906787@bioinfmsc9.bio.ed.ac.uk"
    Entrez.api_key = "37ee5dc16d4c1a46c163da440cc333dca209"

    for i, row in search_df.iterrows():
        terms_gene = f"{row['gene_symbol']}[Gene Name] AND {row['organism']}[Organism]"
        print(f"\nSearching NCBI gene db for {row['ensembl_id']}...")
        
        try:
            search_handle = Entrez.esearch(db = "gene", term = terms_gene, retmax = 20)
            search_record = Entrez.read(search_handle)
            search_handle.close()
            
            entrez_ids = []
            if search_record["IdList"]:
                for gene_id in search_record["IdList"]:
                    entrez_ids.append(gene_id)
                    summary_handle = Entrez.esummary(db = "gene", id = gene_id)
                    summary = Entrez.read(summary_handle)
                    summary_handle.close()
                    
                    result_ncbi = summary["DocumentSummarySet"]["DocumentSummary"][0]
                    genomic_info = result_ncbi.get("GenomicInfo", [])
                    insert_data_ncbi = (
                        row["ensembl_id"], 
                        result_ncbi.get("NomenclatureSymbol", None), 
                        result_ncbi.get("Organism", {}).get("ScientificName", None), 
                        result_ncbi.get("Organism", {}).get("TaxID", None), 
                        genomic_info[0].get("ChrLoc", None) if genomic_info and len(genomic_info) > 0 else None, 
                        genomic_info[0].get("ChrStart", None) if genomic_info and len(genomic_info) > 0 else None, 
                        genomic_info[0].get("ChrStop", None) if genomic_info and len(genomic_info) > 0 else None, 
                        gene_id, 
                        result_ncbi.get("GeneWeight", None), 
                        result_ncbi.get("Summary", None)
                    )
                    cursor.execute(insert_ncbi, insert_data_ncbi)

                print(f"NCBI gene data successfully retrieved for {row['ensembl_id']}!")
            
            else:
                print(f"Entrezpy returns no NCBI gene record for {row['ensembl_id']}.\n")
        
        except Exception as e:
            print(f"Failed to search NCBI gene db for {row['ensembl_id']}. Error: {e}\n")
        
        if entrez_ids:
            print(f"Searching NCBI PubMed db for {row['ensembl_id']}...")
            try:
                for gene_id in entrez_ids:
                    link_handle = Entrez.elink(dbfrom = "gene", db = "pubmed", id = gene_id)
                    link_record = Entrez.read(link_handle)
                    link_handle.close()

                    if link_record[0]["LinkSetDb"]:
                        pubmed_ids = [link["Id"] for link in link_record[0]["LinkSetDb"][0]["Link"]]
                        pubmed_ids = pubmed_ids[:20]

                        doi_dict = {}
                        pubdate_dict = {}
                        author_dict = {}
                        for pubmed_id in pubmed_ids:
                            summary_handle = Entrez.esummary(db = "pubmed", id = pubmed_id)
                            summary = Entrez.read(summary_handle)
                            summary_handle.close()

                            result_pm = summary[0]
                            doi_dict[pubmed_id] = result_pm.get("DOI", None)
                            pubdate_dict[pubmed_id] = result_pm.get("PubDate", None)
                            author_dict[pubmed_id] = result_pm.get("AuthorList", [])

                        print(f"PubMed data successfully retrieved for {gene_id}!")

                    else:
                        print(f"Entrezpy returns no PubMed record for {gene_id}.")

                insert_data_pubmed = (
                    row["ensembl_id"], 
                    json.dumps(pubmed_ids), 
                    json.dumps(doi_dict), 
                    json.dumps(pubdate_dict), 
                    json.dumps(author_dict)
                )
                cursor.execute(insert_pubmed, insert_data_pubmed)

            except Exception as e:
                print(f"Failed to search PubMed db for {row['ensembl_id']}. Error: {e}\n")
            finally:
                print(f"Finished searching NCBI for {row['ensembl_id']}!\n")

    connection.commit()
    cursor = connection.cursor(dictionary = True)
    cursor.execute("SELECT * FROM ncbi_gene_table")
    for row in cursor:
        print(row)
    cursor.execute("SELECT * FROM pubmed_table")
    for row in cursor:
        print(row)
    cursor.close()
    connection.close()

# access_interpro_api() takes in an arg: pandas DataFrame containing search terms.
# Function: Query UniProt database using REST API, then InterPro database using InterPro API.
# Query UniProt database using gene_symbol and organism search terms to get UniProt IDs, which are used to query InterPro database.
# CREATE TABLE interpro_table in SQL.
# If response is ok and there is at least a result, extract and insert relevant fields into interpro_table.
# Print rows of interpro_table to verify changes.
def access_interpro_api(search_df):
    connection = mysql.connector.connect(
        host = "bioinfmsc9.bio.ed.ac.uk", 
        user = "s2906787", 
        password = "r9UqKvEH", 
        database = "s2906787"
    )
    cursor = connection.cursor()
    
    cursor.execute("DROP TABLE IF EXISTS interpro_table")
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS interpro_table (
            ensembl_id VARCHAR(500), 
            uniprot_id VARCHAR(500), 
            interpro_id VARCHAR(500), 
            interpro_name VARCHAR(500), 
            interpro_type VARCHAR(500), 
            go_id VARCHAR(500) 
        )
    """)

    insert_interpro = """
        INSERT INTO interpro_table (
            ensembl_id, 
            uniprot_id, 
            interpro_id, 
            interpro_name, 
            interpro_type, 
            go_id
        ) 
        VALUES (%s, %s, %s, %s, %s, %s)
    """

    uniprot_url = "https://rest.uniprot.org/uniprotkb/search"
    
    for i, row in search_df.iterrows():
        print(f"\nSearching UniProt for {row['ensembl_id']}...")
        params = {
            "query" : f"gene:{row['gene_symbol']} AND organism_name:{row['organism']}", 
            "format" : "json"
        }
        response_uniprot = requests.get(uniprot_url, params = params)

        if response_uniprot.ok:
            uniprot_data = response_uniprot.json()
            
            uniprot_ids = {}
            for entry in uniprot_data["results"]:
                accession = entry.get("primaryAccession")
                uniprot_ids[row['ensembl_id']] = accession

            if len(uniprot_ids) < 1:
                print(f"No UniProt ID can be found for {row['ensembl_id']}.\n")
            
            else:
                interpro_url = "https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot"
                headers = {"Accept" : "application/json"}
                for ens_id, uniprot_id in uniprot_ids.items():
                    response_interpro = requests.get(f"{interpro_url}/{uniprot_id}", headers = headers)

                    if response_interpro.ok:
                        if response_interpro.text.strip():
                            interpro_data = response_interpro.json()
                            
                            for result in interpro_data.get("results", []):
                                metadata = result.get("metadata", {})
                                go_terms = metadata.get("go_terms")

                                if go_terms and isinstance(go_terms, list):
                                    for go_term in go_terms:
                                        insert_data_interpro = (
                                            ens_id, 
                                            uniprot_id, 
                                            metadata.get("accession", None), 
                                            metadata.get("name", None), 
                                            metadata.get("type", None), 
                                            go_term.get("identifier", None)
                                        )
                                        cursor.execute(insert_interpro, insert_data_interpro)

                            print(f"InterPro data successfully retrieved for {uniprot_id}!")

                        else:
                            print(f"Failed to retrieve InterPro data for {uniprot_id}.")

                print(f"Finished searching InterPro for {row['ensembl_id']}\n")

    connection.commit()
    cursor = connection.cursor(dictionary = True)
    cursor.execute("SELECT * FROM interpro_table")
    for row in cursor:
        print(row)
    cursor.close()
    connection.close()

# map_ensembl_pybiomart() takes in an arg: pandas DataFrame containing search terms.
# Function: Query Ensembl database using pybiomart, to get Gene Ontology(GO) IDs for each Ensembl ID.
# Output results to a pandas DataFrame.
# Output: A temporary pandas DataFrame with cols "ensembl_id" and "go_id" for further queries.
def map_ensembl_pybiomart(search_df):
    all_results = []
    for i, row in search_df.iterrows():
        print(f"\nMapping {row['ensembl_id']}...")
        if pd.isna(row["organism"]):
            print(f"Dataset cannot be found for {row['ensembl_id']}.\n")
            continue
        
        prefix_list = row["organism"].split(" ")
        organism_db = f"{prefix_list[0][0].lower()}{prefix_list[-1].lower()}_gene_ensembl"
    
        server = Server(host = "http://www.ensembl.org")
        mart = server["ENSEMBL_MART_ENSEMBL"]
    
        print(f"Searching BioMart for {row['ensembl_id']}...")
        dataset = Dataset(name = organism_db, host = "http://www.ensembl.org")
        result_str = dataset.query(
            attributes = [
                "ensembl_gene_id", 
                "go_id"
            ], 
            filters = {"link_ensembl_gene_id" : row["ensembl_id"]}
        ).to_json(orient = "records")
        
        result_json = json.loads(result_str)
        print(f"Successfully retrieved BioMart results for {row['ensembl_id']}!\n")

        all_results.append(pd.json_normalize(result_json))
    results_df = pd.concat(all_results, ignore_index = True)
    results_df = results_df.rename(columns = {
        "Gene stable ID" : "ensembl_id", 
        "GO term accession" : "go_id"
    })
    return results_df

# access_go_api() takes in an arg: pandas DataFrame containing search terms.
# Function: Query Gene Ontology(GO) database using QuickGO API.
# Query by go_id.
# CREATE TABLE go_table in SQL.
# If response is ok and there is at least a result, extract and insert relevant fields into go_table.
# Print rows of go_table to verify changes.
def access_go_api(search_df):
    connection = mysql.connector.connect(
        host = "bioinfmsc9.bio.ed.ac.uk", 
        user = "s2906787", 
        password = "r9UqKvEH", 
        database = "s2906787"
    )
    cursor = connection.cursor()
    
    cursor.execute("DROP TABLE IF EXISTS go_table")
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS go_table (
            ensembl_id VARCHAR(500), 
            go_id VARCHAR(500), 
            go_name VARCHAR(500), 
            go_aspect VARCHAR(500), 
            go_definition VARCHAR(1000) 
        )
    """)

    insert_go = """
        INSERT INTO go_table (
            ensembl_id, 
            go_id, 
            go_name, 
            go_aspect, 
            go_definition
        ) 
        VALUES (%s, %s, %s, %s, %s)
    """

    base_url = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/search"
    headers = {"Accept" : "application/json"}

    tracker = ""
    for i, row in search_df.iterrows():
        if tracker != row["ensembl_id"]:
            print(f"\nSearching GO for {row['ensembl_id']}...")
            tracker = row["ensembl_id"]
        
        params = {
            "query" : row["go_id"], 
            "limit" : 25, 
            "page" : 1
        }
        response_go = requests.get(base_url, params = params, headers = headers)

        if response_go.ok:
            data_go = response_go.json()
            
            results = data_go.get("results", [])
            if results and len(results) > 0:
                for result in results:
                    if result.get("id") == row["go_id"]:
                        definition = result.get("definition", {})
                        insert_data_go = (
                            row["ensembl_id"], 
                            result.get("id", None), 
                            result.get("name", None), 
                            result.get("aspect", None), 
                            definition.get("text", None) if definition else None
                        )
                        cursor.execute(insert_go, insert_data_go)
                        print(f"Successfully retrieved GO data for {row['go_id']}!")

        else:
            print(f"Failed to retrieve GO data for {row['go_id']}. Error: {response_go.status_code}")

    connection.commit()
    cursor = connection.cursor(dictionary = True)
    cursor.execute("SELECT * FROM go_table")
    for row in cursor:
        print(row)
    cursor.close()
    connection.close()

if __name__ == "__main__":
    main()


