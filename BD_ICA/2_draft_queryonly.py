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
    "ENSMUSG00000009575", 
    "test_error", 
    "ENSMUSG000000095"
]

import subprocess
import glob
import json
import pandas as pd
import requests
from Bio import Entrez
from pybiomart import Server
from pybiomart import Dataset

def main():
    print("\nAccessing Ensembl...\n")
    query_df = access_ensembl_rest(ens_id)
    print("\nDownloaded JSON data for Ensembl results!\n")

    print("\nAccessing NCBI gene...\n")
    access_ncbi_entrezpy(query_df)
    print("\nDownloaded JSON data for NCBI results!\n")
    
    print("\nAccessing InterPro...\n")
    access_interpro_api(query_df)
    print("\nDownloaded JSON data for InterPro results!\n")
    
    print("\nMapping Ensembl IDs with BioMart...\n")
    xrefs_df = map_ensembl_pybiomart(query_df)
    print("\nDownloaded JSON data for BioMart results!\n")

    print("\nAccessing Gene Ontology (GO)...\n")
    access_go_api(xrefs_df)
    print("\nDownloaded JSON data for GO results!\n")

# access_ensembl_rest() takes in an arg: list of Ensembl IDs.
# Function: Query Ensembl using REST API, by looking up the provided list of Ensembl IDs.
# If response is ok and there is at least a result, download the raw JSON-formatted result into ./ensembl_json/ dir.
# Output: ./ensembl_json/{ensembl_id}.json files, each containing raw JSON result from querying each ensembl_id.
# Output: A temporary pandas DataFrame with cols "ensembl_id", "gene_symbol", "organism" for further queries.
def access_ensembl_rest(ens_id_list):
    subprocess.call("rm -rf ensembl_json", shell = True)
    subprocess.call("mkdir -p ensembl_json", shell = True)
    base_url = "https://rest.ensembl.org/"
    headers = {"Accept" : "application/json"}
    
    all_results = []
    for ens_id in ens_id_list:
        print(f"\nSearching Ensembl for {ens_id}...")
        response_lookup = requests.get(f"{base_url}lookup/id/{ens_id}", headers = headers)
        
        if response_lookup.ok:
            data = response_lookup.json()
            with open(f"./ensembl_json/{ens_id}.json", "w") as file:
                json.dump(data, file, indent = 2)
            
            filter_data = {
                "ensembl_id" : data.get("id"), 
                "gene_symbol" : data.get("display_name"), 
                "organism" : data.get("species")
            }
            all_results.append(filter_data)
            print(f"Ensembl data successfully retrieved for {ens_id}!\n")

        else:
            print(f"Failed to retrieve Ensembl data for {ens_id}. Error: {response_lookup.status_code}\n")
            all_results.append({
                "ensembl_id" : ens_id, 
                "gene_symbol" : None, 
                "organism" : None
            })

    results_df = pd.DataFrame(all_results)

    def format_species_name(x):
        if pd.isna(x):
            return x
        x_list = x.split("_")
        x_list[0] = x_list[0].title()
        return " ".join(x_list)

    results_df["organism"] = results_df.apply(lambda x: format_species_name(x["organism"]), axis = 1)
    return results_df

# access_ncbi_entrezpy() takes in an arg: pandas DataFrame containing search terms.
# Function: Query NCBI gene and PubMed databases using entrezpy E-utilities.
# Entrez.esearch to get Entrez IDs for each search term.
# Entrez.elink to get PubMed IDs for each Entrez ID.
# Entrez.esummary to get summarised NCBI gene and PubMed results for each Entrez ID.
# If response is ok and there is at least a result, download the raw JSON-formatted result into ./ncbi_json/ dir.
# Separate the results of NCBI gene and PubMed queries into ./ncbi_json/ncbi_gene/ and ./ncbi_json/pubmed/ dirs respectively.
# Results are grouped by ensembl_id, and are queried by entrez_id and pubmed_id.
# Output: ./ncbi_json/ncbi_gene/{ensembl_id}/{entrez_id}.json and ./ncbi_json/pubmed/{ensembl_id}/{pubmed_id}.json files, 
# each containing raw JSON result from querying each entrez_id or pubmed_id.
def access_ncbi_entrezpy(search_df):
    subprocess.call("rm -rf ncbi_json", shell = True)
    Entrez.email = "s2906787@bioinfmsc9.bio.ed.ac.uk"
    
    for i, row in search_df.iterrows():
        subprocess.call(f"mkdir -p ncbi_json/ncbi_gene/{row['ensembl_id']}", shell = True)
        subprocess.call(f"mkdir -p ncbi_json/pubmed/{row['ensembl_id']}", shell = True)
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
                    
                    result = summary["DocumentSummarySet"]["DocumentSummary"][0]
                    with open(f"./ncbi_json/ncbi_gene/{row['ensembl_id']}/{gene_id}.json", "w") as file:
                        json.dump(result, file, indent = 2)

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

                        for pubmed_id in pubmed_ids:
                            summary_handle = Entrez.esummary(db = "pubmed", id = pubmed_id)
                            summary = Entrez.read(summary_handle)
                            summary_handle.close()

                            result = summary[0]
                            with open(f"./ncbi_json/pubmed/{row['ensembl_id']}/{pubmed_id}.json", "w") as file:
                                json.dump(result, file, indent = 2)
                            
                        print(f"PubMed data successfully retrieved for {gene_id}!")

                    else:
                        print(f"Entrezpy returns no PubMed record for {gene_id}.")

            except Exception as e:
                print(f"Failed to search PubMed db for {row['ensembl_id']}. Error: {e}\n")
            finally:
                print(f"Finished searching NCBI PubMed db for {row['ensembl_id']}!\n")

# access_interpro_api() takes in an arg: pandas DataFrame containing search terms.
# Function: Query UniProt database using REST API, then InterPro database using InterPro API.
# Query UniProt database to get UniProt IDs, which are used to query InterPro database.
# If response is ok and there is at least a result, download the raw JSON-formatted result into ./interpro_json/ dir.
# Results are grouped by ensembl_id, and are queried by uniprot_id.
# Output: ./interpro_json/{ensembl_id}/{uniprot_id}.json files, each containing raw JSON result from querying each uniprot_id.
def access_interpro_api(search_df):
    subprocess.call("rm -rf interpro_json", shell = True)
    uniprot_url = "https://rest.uniprot.org/uniprotkb/search"
    
    for i, row in search_df.iterrows():
        subprocess.call(f"mkdir -p interpro_json/{row['ensembl_id']}", shell = True)
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
                            with open(f"./interpro_json/{ens_id}/{uniprot_id}.json", "w") as file:
                                json.dump(interpro_data, file, indent = 2)
                            print(f"InterPro data successfully retrieved for {uniprot_id}!")

                        else:
                            print(f"Failed to retrieve InterPro data for {uniprot_id}.")

                print(f"Finished searching InterPro for {row['ensembl_id']}\n")

# map_ensembl_pybiomart() takes in an arg: pandas DataFrame containing search terms.
# Function: Query Ensembl database using pybiomart, to get Gene Ontology(GO) IDs for each Ensembl ID.
# If response is ok and there is at least a result, download the raw JSON-formatted result into ./xrefs_json/ dir.
# Output: ./xrefs_json/{ensembl_id}.json files, each containing raw JSON result from querying each ensembl_id.
# Output: A temporary pandas DataFrame with cols "ensembl_id" and "go_id" for further queries.
def map_ensembl_pybiomart(search_df):
    subprocess.call("rm -rf xrefs_json", shell = True)
    subprocess.call("mkdir -p xrefs_json", shell = True)
    
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
        with open(f"./xrefs_json/{row['ensembl_id']}.json", "w") as file:
            json.dump(result_json, file, indent = 2)
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
# If response is ok and there is at least a result, download the raw JSON-formatted result into ./go_json/ dir.
# Results are grouped by ensembl_id, and are queried by go_id.
# Output: ./go_json/{ensembl_id}/{go_id}.json files, each containing raw JSON result from querying each go_id.
def access_go_api(search_df):
    subprocess.call("rm -rf go_json", shell = True)
    subprocess.call("mkdir -p go_json", shell = True)
    base_url = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/search"
    headers = {"Accept" : "application/json"}

    tracker = ""
    for i, row in search_df.iterrows():
        if tracker != row["ensembl_id"]:
            print(f"\nSearching GO for {row['ensembl_id']}...")
            tracker = row["ensembl_id"]
            subprocess.call(f"mkdir -p ./go_json/{tracker}", shell = True)
        
        params = {
            "query" : row["go_id"], 
            "limit" : 25, 
            "page" : 1
        }
        response_go = requests.get(base_url, params = params, headers = headers)

        if response_go.ok:
            data_go = response_go.json()
            file_name = f"{row['go_id']}".replace(":", "_")
            with open(f"./go_json/{tracker}/{file_name}.json", "w") as file:
                json.dump(data_go, file, indent = 2)
            print(f"Successfully retrieved GO data for {row['go_id']}!")

        else:
            print(f"Failed to retrieve GO data for {row['go_id']}. Error: {response_go.status_code}")

if __name__ == "__main__":
    main()

