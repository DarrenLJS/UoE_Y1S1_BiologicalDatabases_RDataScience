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
from pybiomart import Server
from pybiomart import Dataset

def main():
    print("\nAccessing Ensembl...\n")
    query_df = access_ensembl_rest(ens_id)
    print("\nDownloaded Ensembl results into SQL!\n")

    print("\nAccessing NCBI...\n")
    access_ncbi_entrez(query_df)
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
    
    cursor.execute("DROP TABLE IF EXISTS all_results")
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS all_results AS 
        SELECT 
        COALESCE(gi.ensembl_id, ne.ensembl_id) AS ensembl_id, 
        ne.gene_symbol, 
        ne.organism, 
        ne.taxon_id, 
        ne.description, 
        ne.chromosome, 
        ne.strand, 
        ne.start, 
        ne.end, 
        ne.entrez_id, 
        ne.weight, 
        ne.summary, 
        gi.uniprot_id, 
        gi.interpro_id, 
        gi.interpro_name, 
        gi.interpro_type, 
        gi.go_id, 
        gi.go_name, 
        gi.go_aspect, 
        gi.go_definition 
        FROM go_interpro gi 
        LEFT JOIN ncbi_ensembl ne 
        ON gi.ensembl_id = ne.ensembl_id 
        UNION 
        SELECT 
        COALESCE(gi.ensembl_id, ne.ensembl_id) AS ensembl_id, 
        ne.gene_symbol, 
        ne.organism, 
        ne.taxon_id, 
        ne.description, 
        ne.chromosome, 
        ne.strand, 
        ne.start, 
        ne.end, 
        ne.entrez_id, 
        ne.weight, 
        ne.summary, 
        gi.uniprot_id, 
        gi.interpro_id, 
        gi.interpro_name, 
        gi.interpro_type, 
        gi.go_id, 
        gi.go_name, 
        gi.go_aspect, 
        gi.go_definition 
        FROM go_interpro gi 
        RIGHT JOIN ncbi_ensembl ne 
        ON gi.ensembl_id = ne.ensembl_id
    """)
    
    cursor.execute("DROP TABLE IF EXISTS go_interpro")
    cursor.execute("DROP TABLE IF EXISTS ncbi_ensembl")
    connection.commit()
    
    df_final = pd.read_sql("SELECT * FROM all_results", connection)
    
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

# access_ncbi_entrez() takes in an arg: pandas DataFrame containing search terms.
# Function: Query NCBI gene database using entrez E-utilities.
# Esearch to get Entrez IDs for each search term.
# Esummary to get summarised NCBI gene for each Entrez ID, 
# CREATE TABLE ncbi_gene_table in SQL.
# If response is ok and there is at least a result, extract and insert relevant fields into ncbi_gene_table.
# Print rows of ncbi_gene_table to verify changes.
def access_ncbi_entrez(search_df):
    connection = mysql.connector.connect(
        host = "bioinfmsc9.bio.ed.ac.uk", 
        user = "s2906787", 
        password = "r9UqKvEH", 
        database = "s2906787"
    )
    cursor = connection.cursor()
    
    cursor.execute("DROP TABLE IF EXISTS ncbi_gene_table")
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
    
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    email = "s2906787@bioinfmsc9@bio.ed.ac.uk"
    api_key = "37ee5dc16d4c1a46c163da440cc333dca209"

    for i, row in search_df.iterrows():
        terms_gene = f"{row['gene_symbol']}[Gene Name] AND {row['organism']}[Organism]"
        params_esearch = {
            "db" : "gene", 
            "term" : terms_gene, 
            "retmax" : 20, 
            "retmode" : "json", 
            "email" : email, 
            "api_key" : api_key
        }
        print(f"\nSearching NCBI gene db for {row['ensembl_id']}...")
        
        try:
            response_esearch = requests.get(f"{base_url}esearch.fcgi", params = params_esearch)
            search_result = response_esearch.json()
            id_list = search_result.get("esearchresult", {}).get("idlist", [])
            
            if id_list:
                for gene_id in id_list:
                    params_esummary = {
                        "db" : "gene", 
                        "id" : gene_id, 
                        "retmode" : "json", 
                        "email" : email, 
                        "api_key" : api_key
                    }
                    response_esummary = requests.get(f"{base_url}esummary.fcgi", params = params_esummary)
                    summary = response_esummary.json()
                    
                    result_ncbi = summary_result.get("result", {}).get(gene_id, {})
                    genomic_info = result_ncbi.get("genomicinfo", [])
                    insert_data_ncbi = (
                        row["ensembl_id"], 
                        result_ncbi.get("nomenclaturesymbol", None), 
                        result_ncbi.get("organism", {}).get("scientificname", None), 
                        result_ncbi.get("organism", {}).get("taxid", None), 
                        genomic_info[0].get("chrloc", None) if genomic_info else None, 
                        genomic_info[0].get("chrstart", None) if genomic_info else None, 
                        genomic_info[0].get("chrstop", None) if genomic_info else None, 
                        gene_id, 
                        result_ncbi.get("geneweight", None), 
                        result_ncbi.get("summary", None)
                    )
                    cursor.execute(insert_ncbi, insert_data_ncbi)

                print(f"NCBI gene data successfully retrieved for {row['ensembl_id']}!")
            
            else:
                print(f"Entrez returns no NCBI gene record for {row['ensembl_id']}.\n")
        
        except Exception as e:
            print(f"Failed to search NCBI gene db for {row['ensembl_id']}. Error: {e}\n")

    connection.commit()
    cursor = connection.cursor(dictionary = True)
    cursor.execute("SELECT * FROM ncbi_gene_table")
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


