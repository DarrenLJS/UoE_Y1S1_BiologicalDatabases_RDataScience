#!/usr/bin/python3

ens_ids = [
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
import time
import subprocess
from pybiomart import Server
from pybiomart import Dataset

# main() takes in an arg: a list of Ensembl IDs.
# Function: Creates MySQL database and runs the pipeline.
# Runs access_ensembl_rest, access_interpro_api, map_ensembl_pybiomart, access_go_api functions.
# Once final_integration and final_summary tables are produced in MySQL, parse the tables from MySQL to Pandas DataFrame.
# Output 1: BD_ICA_integration_table.tsv file, containing final_integration table in tab-separated values parsed from Pandas DataFrame.
# Output 2: BD_ICA_summary_table.tsv file, containing final_summary table in tab-separated values parsed from Pandas DataFrame.
def main(ens_id_list):
    
    # Establish MySQL connection and cursor.
    # Create a new database if not exists.
    # Close cursor and connection.
    host = "bioinfmsc9.bio.ed.ac.uk"
    user = "s2906787"
    pwd = "r9UqKvEH"
    database_name = "s2906787"
    print(f"Creating MySQL database {database_name}...")
    connection = mysql.connector.connect(
        host = host, 
        user = user, 
        password = pwd
    )
    cursor = connection.cursor()
    cursor.execute(f"CREATE DATABASE IF NOT EXISTS {database_name}")
    cursor.execute(f"USE {database_name}")
    print(f"Database {database_name} successfully created or already exists!")
    cursor.close()
    connection.close()
    
    # Pipeline
    print("\nAccessing Ensembl...\n")
    query_df = access_ensembl_rest(ens_id_list, host, user, pwd, database_name)
    print("\nDownloaded Ensembl results into SQL!\n")
    
    print("\nAccessing InterPro...\n")
    access_interpro_api(query_df, host, user, pwd, database_name)
    print("\nDownloaded InterPro results into SQL!\n")
    
    print("\nMapping Ensembl IDs with BioMart...\n")
    xrefs_df = map_ensembl_pybiomart(query_df)
    print("\nMapped Ensembl IDs to GO IDs!\n")

    print("\nAccessing Gene Ontology (GO)...\n")
    access_go_api(xrefs_df, host, user, pwd, database_name)
    print("\nDownloaded GO results into SQL!\n")

    print("\nIntegrating SQL tables into one...\n")
    integration_df, summary_df = integrate_sql(ens_id_list, host, user, pwd, database_name)
    
    # After parsing MySQL final_integration table to Pandas DataFrame, reformat columns with multiple values per cell into lists.
    cols_format = [
        "uniprot_id", 
        "interpro_id", 
        "go_id", 
        "molecular_function", 
        "biological_process", 
        "cellular_component", 
        "protein_domain_family", 
        "all_interpro_features"
    ]
    for col in cols_format:
        if col in integration_df.columns:
            integration_df[col] = integration_df[col].apply(lambda x: x.strip().split("; ") if isinstance(x, str) else x)
    
    # Save final_integration and final_summary tables to .tsv files.
    # Print the first few lines of each table to the terminal.
    integration_df.to_csv("BD_ICA_integration_table.tsv", sep = "\t", header = True, index = False)
    summary_df.to_csv("BD_ICA_summary_table.tsv", sep = "\t", header = True, index = False)
    subprocess.call("head BD_ICA_integration_table.tsv", shell = True)
    subprocess.call("head BD_ICA_summary_table.tsv", shell = True)
    print("\nFinal integration completed! Output: BD_ICA_integration_table.tsv, BD_ICA_summary_table.tsv\n")

# integrate_sql() takes in args: a list of Ensembl IDs, and login details to MySQL database.
# Function: Integrates all the downloaded data tables into one integrated table in MySQL.
# Firstly, FULL OUTER JOIN all data tables into one all_results table.
# Next, concat and group data by ensembl_id, into final_integration table.
# Lastly, create final_summary table of ensembl_id, uniprot_id, interpro_id, and go_id, 
# excluding those protein_coding with NULL uniprot_id anf interpro_id.
# Output: two Pandas DataFrames parsed from MySQL final_integration and final_summary tables.
def integrate_sql(ens_id_list, server, username, password, database):
    
    # Establish MySQL connection and cursor.
    connection = mysql.connector.connect(
        host = server, 
        user = username, 
        password = password, 
        database = database
    )
    cursor = connection.cursor()
    
    # FULL OUTER JOIN go_table and interpro_table into go_interpro table.
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
    
    # FULL OUTER JOIN go_interpro and ensembl_table into all_results table.
    cursor.execute("DROP TABLE IF EXISTS all_results")
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS all_results AS 
        SELECT 
        COALESCE(gi.ensembl_id, e.ensembl_id) AS ensembl_id, 
        e.gene_symbol, 
        e.organism, 
        e.description, 
        e.gene_type, 
        e.chromosome, 
        e.strand, 
        e.start, 
        e.end, 
        gi.uniprot_id, 
        gi.interpro_id, 
        gi.interpro_name, 
        gi.interpro_type, 
        gi.go_id, 
        gi.go_name, 
        gi.go_aspect, 
        gi.go_definition 
        FROM go_interpro gi 
        LEFT JOIN ensembl_table e 
        ON gi.ensembl_id = e.ensembl_id 
        UNION 
        SELECT 
        COALESCE(gi.ensembl_id, e.ensembl_id) AS ensembl_id, 
        e.gene_symbol, 
        e.organism, 
        e.description, 
        e.gene_type, 
        e.chromosome, 
        e.strand, 
        e.start, 
        e.end, 
        gi.uniprot_id, 
        gi.interpro_id, 
        gi.interpro_name, 
        gi.interpro_type, 
        gi.go_id, 
        gi.go_name, 
        gi.go_aspect, 
        gi.go_definition 
        FROM go_interpro gi 
        RIGHT JOIN ensembl_table e 
        ON gi.ensembl_id = e.ensembl_id
    """)
    
    # Concat and group data by ensembl_id from all_results to final_integration table.
    # Sort by the original order of ens_id_list.
    sorter = ", ".join([f"'{ens_id}'" for ens_id in ens_id_list])
    cursor.execute("DROP TABLE IF EXISTS final_integration")
    cursor.execute(f"""
        CREATE TABLE IF NOT EXISTS final_integration AS 
        SELECT 
        ensembl_id, 
        gene_symbol, 
        organism, 
        description, 
        gene_type, 
        GROUP_CONCAT(
            DISTINCT uniprot_id ORDER BY uniprot_id SEPARATOR "; "
        ) AS uniprot_id, 
        GROUP_CONCAT(
            DISTINCT interpro_id ORDER BY interpro_id SEPARATOR "; "
        ) AS interpro_id, 
        GROUP_CONCAT(
            DISTINCT go_id ORDER BY go_id SEPARATOR "; "
        ) AS go_id, 
        GROUP_CONCAT(
            DISTINCT CASE 
                WHEN go_aspect = 'molecular_function' THEN go_name 
            END 
            ORDER BY go_name 
            SEPARATOR '; '
        ) AS molecular_function, 
        GROUP_CONCAT(
            DISTINCT CASE 
                WHEN go_aspect = 'biological_process' THEN go_name 
            END 
            ORDER BY go_name 
            SEPARATOR '; '
        ) AS biological_process, 
        GROUP_CONCAT(
            DISTINCT CASE 
                WHEN go_aspect = 'cellular_component' THEN go_name 
            END 
            ORDER BY go_name 
            SEPARATOR '; '
        ) AS cellular_component, 
        GROUP_CONCAT(
            DISTINCT CASE 
                WHEN interpro_type IN ('domain', 'family') THEN interpro_name 
            END 
            ORDER BY interpro_name 
            SEPARATOR '; '
        ) AS protein_domain_family, 
        GROUP_CONCAT(
            DISTINCT CONCAT(interpro_name, '(', interpro_type, ')') 
            ORDER BY interpro_name 
            SEPARATOR '; '
        ) AS all_interpro_features, 
        chromosome, 
        strand, 
        start, 
        end 
        FROM all_results 
        GROUP BY ensembl_id, gene_symbol, organism, description, gene_type, chromosome, strand, start, end 
        ORDER BY FIELD(ensembl_id, {sorter})
    """)
    
    # Create a final_summary table of ensembl_id, uniprot_id, interpro_id, and go_id, 
    # excluding those protein_coding with NULL uniprot_id and interpro_id.
    sorter = ", ".join([f"'{ens_id}'" for ens_id in ens_id_list])
    cursor.execute("DROP TABLE IF EXISTS final_summary")
    cursor.execute(f"""
        CREATE TABLE IF NOT EXISTS final_summary AS 
        SELECT DISTINCT 
        ensembl_id, 
        gene_symbol, 
        organism, 
        gene_type, 
        uniprot_id, 
        interpro_id, 
        go_id 
        FROM all_results 
        WHERE NOT (
            gene_type = 'protein_coding' 
            AND uniprot_id IS NULL 
            AND interpro_id IS NULL
        ) 
        ORDER BY FIELD(ensembl_id, {sorter})
    """)

    # Drop go_interpro table. Commit changes.
    # Parse final_integration and final_summary tables from MySQL to Pandas DataFrame.
    # Close cursor and connection.
    # Return parsed Pandas DataFrames.
    cursor.execute("DROP TABLE IF EXISTS go_interpro")
    connection.commit()
    
    df_final = pd.read_sql("SELECT * FROM final_integration", connection)
    df_summary = pd.read_sql("SELECT * FROM final_summary", connection)

    cursor.close()
    connection.close()
    return df_final, df_summary

# access_ensembl_rest() takes in args: a list of Ensembl IDs, and login details to MySQL database.
# Function: Queries Ensembl using REST API, by looking up the provided list of Ensembl IDs.
# Extracts and downloads results into MySQL, creating ensembl_table table.
# If response is ok and there is at least a result, extract and insert relevant fields into ensembl_table.
# Print rows of ensembl_table to verify changes.
# Output: A temporary Pandas DataFrame with cols "ensembl_id", "gene_symbol", "organism" for further queries.
def access_ensembl_rest(ens_id_list, server, username, password, database):
    
    # Establish MySQL connection and cursor.
    # Create new table ensembl_table.
    # Prepare MySQL insert query, insert_ens.
    connection = mysql.connector.connect(
        host = server, 
        user = username, 
        password = password, 
        database = database
    )
    cursor = connection.cursor()
    
    cursor.execute("DROP TABLE IF EXISTS ensembl_table")
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS ensembl_table (
            ensembl_id VARCHAR(500), 
            gene_symbol VARCHAR(500), 
            organism VARCHAR(500), 
            description VARCHAR(500), 
            gene_type VARCHAR(500), 
            chromosome VARCHAR(10), 
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
            gene_type, 
            chromosome, 
            strand, 
            start, 
            end
        ) 
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)
    """
    
    # Query Ensembl by looking up the provided list of Ensembl IDs.
    base_url = "https://rest.ensembl.org/"
    headers = {"Accept" : "application/json"}
    search_terms = []
    for ens_id in ens_id_list:
        print(f"\nSearching Ensembl for {ens_id}...")
        time.sleep(0.5)
        response_lookup = requests.get(f"{base_url}lookup/id/{ens_id}", headers = headers)
        
        # If response is ok, extract and download data into MySQL, by inserting into ensembl_table.
        # Also, extract search terms for Pandas DataFrame.
        # Else, print error and insert None values.
        if response_lookup.ok:
            data_ens = response_lookup.json()
            insert_data = (
                data_ens.get("id", None), 
                data_ens.get("display_name", None), 
                data_ens.get("species", None), 
                data_ens.get("description", None), 
                data_ens.get("biotype", None), 
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
            insert_data = (ens_id, None, None, None, None, None, None, None, None)
            cursor.execute(insert_ens, insert_data)

            search_terms.append({
                "ensembl_id" : ens_id, 
                "gene_symbol" : None, 
                "organism" : None
            })
    
    # Commit changes to MySQL.
    # Print rows of ensembl_table to verify changes.
    # Close cursor and connection.
    connection.commit()
    cursor = connection.cursor(dictionary = True)
    cursor.execute("SELECT * FROM ensembl_table")
    for row in cursor:
        print(row)
    cursor.close()
    connection.close()
    
    # Create a Pandas DataFrame of search terms for each Ensembl ID.
    search_df = pd.DataFrame(search_terms)
    
    # Reformat "organism" search term from e.g. mus_musculus to Mus musculus
    def format_species_name(x):
        if pd.isna(x):
            return x
        x_list = x.split("_")
        x_list[0] = x_list[0].title()
        return " ".join(x_list)
    search_df["organism"] = search_df.apply(lambda x: format_species_name(x["organism"]), axis = 1)
    
    return search_df

# access_interpro_api() takes in args: a Pandas DataFrame containing search terms, and login details to MySQL database.
# Function: Queries UniProt database using REST API, then InterPro database using InterPro API.
# Extracts and downloads results in MySQL, creating interpro_table table.
# Query UniProt database using gene_symbol and organism search terms to get UniProt IDs, which are used to query InterPro database.
# If response is ok and there is at least a result, extract and insert relevant fields into interpro_table.
# Print rows of interpro_table to verify changes.
def access_interpro_api(search_df, server, username, password, database):
    
    # Establish MySQL connection and cursor.
    # Create new table interpro_table.
    # Prepare MySQL insert query, insert_interpro.
    connection = mysql.connector.connect(
        host = server, 
        user = username, 
        password = password, 
        database = database
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
    
    # Query UniProt by looking up gene_symbol and organism search terms from provided Pandas DataFrame
    uniprot_url = "https://rest.uniprot.org/uniprotkb/search"
    for i, row in search_df.iterrows():
        print(f"\nSearching UniProt for {row['ensembl_id']}...")
        params = {
            "query" : f"gene:{row['gene_symbol']} AND organism_name:{row['organism']}", 
            "format" : "json"
        }
        time.sleep(0.5)
        response_uniprot = requests.get(uniprot_url, params = params)
        
        # If response from UniProt is ok, retrieve and check for UniProt ID.
        if response_uniprot.ok:
            uniprot_data = response_uniprot.json()
            uniprot_ids = {}
            for entry in uniprot_data["results"]:
                accession = entry.get("primaryAccession")
                uniprot_ids[row['ensembl_id']] = accession

            # If no UniProt ID, print error.
            if len(uniprot_ids) < 1:
                print(f"No UniProt ID can be found for {row['ensembl_id']}.\n")
            
            # If have UniProt IDs, query InterPro by looking up the UniProt IDs.
            else:
                interpro_url = "https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot"
                headers = {"Accept" : "application/json"}
                for ens_id, uniprot_id in uniprot_ids.items():
                    time.sleep(0.5)
                    response_interpro = requests.get(f"{interpro_url}/{uniprot_id}", headers = headers)
                    
                    # If response from InterPro is ok, retrieve and check for GO IDs.
                    # Else, insert ens_id and uniprot_id only.
                    if response_interpro.ok:
                        if response_interpro.text.strip():
                            interpro_data = response_interpro.json()
                            
                            for result in interpro_data.get("results", []):
                                metadata = result.get("metadata", {})
                                go_terms = metadata.get("go_terms")
                                
                                # If have GO IDs, extract and download data into MySQL, inserting into interpro_table.
                                # Else, extract and download relevant data. Insert go_id as None.
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

                                else:
                                    insert_data_interpro = (
                                        ens_id,
                                        uniprot_id,
                                        metadata.get("accession", None),
                                        metadata.get("name", None),
                                        metadata.get("type", None),
                                        None
                                    )
                                    cursor.execute(insert_interpro, insert_data_interpro)

                            print(f"InterPro data successfully retrieved for {uniprot_id}!")

                        else:
                            print(f"Failed to retrieve InterPro data for {uniprot_id}.")
                            insert_data_interpro = (
                                ens_id, 
                                uniprot_id, 
                                None, 
                                None, 
                                None, 
                                None
                            )
                            cursor.execute(insert_interpro, insert_data_interpro)

                print(f"Finished searching InterPro for {row['ensembl_id']}!\n")
    
    # Commit changes.
    # Print rows of interpro_table to verify changes.
    # Close cursor and connection.
    connection.commit()
    cursor = connection.cursor(dictionary = True)
    cursor.execute("SELECT * FROM interpro_table")
    for row in cursor:
        print(row)
    cursor.close()
    connection.close()

# map_ensembl_pybiomart() takes in an arg: a Pandas DataFrame containing search terms.
# Function: Queries Ensembl using pybiomart, to get Gene Ontology(GO) IDs for each Ensembl ID.
# Save results to a Pandas DataFrame.
# Output: A temporary Pandas DataFrame with cols "ensembl_id" and "go_id" for further queries.
def map_ensembl_pybiomart(search_df):
    
    # In BioMart, have to select dataset based on organism first.
    all_results = []
    for i, row in search_df.iterrows():
        print(f"\nMapping {row['ensembl_id']}...")

        # If organism is NaN, print error and continue.
        # Else, convert organism to database search term e.g. Mus musculus to mmusculus.
        if pd.isna(row["organism"]):
            print(f"Dataset cannot be found for {row['ensembl_id']}.\n")
            continue
        prefix_list = row["organism"].split(" ")
        join_letters = "".join(list(map(lambda x: x[0], prefix_list[:-1])))
        organism_db = f"{join_letters.lower()}{prefix_list[-1].lower()}_gene_ensembl"
        
        # Query Ensembl through BioMart using organism dataset and ensembl_id.
        # Extract associated GO IDs.
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
    
    # Create a Pandas DataFrame with ensembl_id and go_id, and return it.
    results_df = pd.concat(all_results, ignore_index = True)
    results_df = results_df.rename(columns = {
        "Gene stable ID" : "ensembl_id", 
        "GO term accession" : "go_id"
    })
    return results_df

# access_go_api() takes in args: a Pandas DataFrame containing search terms, and login details to MySQL database.
# Function: Queries Gene Ontology(GO) database using QuickGO API, by looking up GO IDs.
# Extracts and downloads data into MySQL, creating go_table table.
# If response is ok and there is at least a result, extract and insert relevant fields into go_table.
# Print rows of go_table to verify changes.
def access_go_api(search_df, server, username, password, database):
    
    # Establish MySQL connection and cursor.
    # Create new table go_table.
    # Prepare MySQL insert query, insert_go.
    connection = mysql.connector.connect(
        host = server, 
        user = username, 
        password = password, 
        database = database
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
    
    # Query QuickGO by looking up go_id from provided Pandas DataFrame.
    base_url = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/search"
    headers = {"Accept" : "application/json"}
    
    # Note corresponding ensembl_id
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
        
        # If response is ok, check if there are results.
        # If there are results, extract and download data into MySQL, inserting into go_table.
        # Else, print error.
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
                print(f"Failed to retrieve GO data for {row['go_id']}.")

        else:
            print(f"Failed to retrieve GO data for {row['go_id']}. Error: {response_go.status_code}")

    # Commit changes.
    # Print rows of go_table to verify changes.
    # Close cursor and connection.
    connection.commit()
    cursor = connection.cursor(dictionary = True)
    cursor.execute("SELECT * FROM go_table")
    for row in cursor:
        print(row)
    cursor.close()
    connection.close()

if __name__ == "__main__":
    main(ens_ids)


