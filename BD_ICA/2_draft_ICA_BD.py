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

import os
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
    
    print("\nMerging data into pandas DataFrame...\n")
    ensembl_df = ensembl_cleaner()
    ncbi_gene_df = ncbi_gene_cleaner()
    pubmed_df = pubmed_cleaner()
    interpro_df = interpro_cleaner()
    go_df = go_cleaner()
    final_df = final_merge(ensembl_df, ncbi_gene_df, pubmed_df, interpro_df, go_df)
    print("\nData merging completed! Pandas DataFrame ready for parsing to SQL!\n")
    print(final_df)

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

# ensembl_cleaner() takes in no arg.
# Function: Extracts Ensembl gene data from ./ensembl_json/*.json, skipping empty files.
# Extracts relevant fields, and reformats "organism" values to proper format.
# Output: pandas DataFrame containing extracted Ensembl data.
def ensembl_cleaner():
    data = []
    json_files = glob.glob("./ensembl_json/*.json")
    
    for file_path in json_files:
        try:
            with open(file_path, "r") as f:
                content = json.load(f)
            if not content:
                continue
                
            row = {
                "ensembl_id" : content.get("id", pd.NA),
                "gene_symbol" : content.get("display_name", pd.NA),
                "organism" : content.get("species", pd.NA),
                "description" : content.get("description", pd.NA),
                "chromosome" : content.get("seq_region_name", pd.NA),
                "strand" : content.get("strand", pd.NA),
                "start" : content.get("start", pd.NA),
                "end" : content.get("end", pd.NA)
            }
            data.append(row)
            
        except Exception  as e:
            print(f"Error processing {file_path}: {e}")
            continue
    
    df = pd.DataFrame(data)
    def format_species_name(x):
        if pd.isna(x):
            return pd.NA
        x_list = x.split("_")
        x_list[0] = x_list[0].title()
        return " ".join(x_list)

    df["organism"] = df.apply(lambda x: format_species_name(x["organism"]), axis = 1)
    return df

# ncbi_gene_cleaner() takes in no arg.
# Function: Extracts NCBI gene data from ./ncbi_json/ncbi_gene/*/*.json, skipping empty files.
# Extracts relevant fields, and extracts entrez_id from file name and ensembl_id from dir name.
# Output: pandas DataFrame containing extracted NCBI gene data.
def ncbi_gene_cleaner():
    data = []
    ensembl_dirs = glob.glob("./ncbi_json/ncbi_gene/*")
    
    for ensembl_dir in ensembl_dirs:
        if not os.path.isdir(ensembl_dir):
            continue
            
        ensembl_id = os.path.basename(ensembl_dir)
        json_files = glob.glob(os.path.join(ensembl_dir, "*.json"))
        
        for file_path in json_files:
            try:
                with open(file_path, "r") as f:
                    content = json.load(f)
                if not content:
                    continue
                
                entrez_id = os.path.splitext(os.path.basename(file_path))[0]
                organism = content.get("Organism", {})
                scientific_name = organism.get("ScientificName", pd.NA)
                taxon_id = organism.get("TaxID", pd.NA)
                genomic_info = content.get("GenomicInfo", [])
                
                chr_loc = pd.NA
                chr_start = pd.NA
                chr_stop = pd.NA
                if genomic_info and len(genomic_info) > 0:
                    chr_loc = genomic_info[0].get("ChrLoc", pd.NA)
                    chr_start = genomic_info[0].get("ChrStart", pd.NA)
                    chr_stop = genomic_info[0].get("ChrStop", pd.NA)
                
                row = {
                    "ensembl_id" : ensembl_id,
                    "gene_symbol" : content.get("NomenclatureSymbol", pd.NA),
                    "organism" : scientific_name,
                    "taxon_id" : taxon_id,
                    "chromosome" : chr_loc,
                    "start" : chr_start,
                    "end" : chr_stop,
                    "entrez_id" : entrez_id,
                    "weight" : content.get("GeneWeight", pd.NA),
                    "summary" : content.get("Summary", pd.NA)
                }
                data.append(row)
                
            except Exception as e:
                print(f"Error processing {file_path}: {e}")
                continue
    
    return pd.DataFrame(data)

# pubmed_cleaner() takes in no arg.
# Function: Extracts PubMed data from ./ncbi_json/pubmed/*/*.json, skipping empty files.
# Extracts relevant fields, and extracts pubmed_id from file name and ensembl_id from dir name.
# Aggregate data by ensembl_id.
# Output: pandas DataFrame containing extracted PubMed data.
def pubmed_cleaner():
    data_dict = {}
    ensembl_dirs = glob.glob("./ncbi_json/pubmed/*")
    
    for ensembl_dir in ensembl_dirs:
        if not os.path.isdir(ensembl_dir):
            continue
            
        ensembl_id = os.path.basename(ensembl_dir)
        json_files = glob.glob(os.path.join(ensembl_dir, "*.json"))
        
        if ensembl_id not in data_dict:
            data_dict[ensembl_id] = {
                "ensembl_id" : ensembl_id,
                "pubmed_id" : [],
                "doi" : {},
                "publication_date" : {},
                "author_list" : {}
            }
        
        for file_path in json_files:
            try:
                with open(file_path, "r") as f:
                    content = json.load(f)
                if not content:
                    continue
                
                pubmed_id = os.path.splitext(os.path.basename(file_path))[0]
                data_dict[ensembl_id]["pubmed_id"].append(pubmed_id)
                data_dict[ensembl_id]["doi"][pubmed_id] = content.get("DOI", pd.NA)
                data_dict[ensembl_id]["publication_date"][pubmed_id] = content.get("PubDate", pd.NA)
                data_dict[ensembl_id]["author_list"][pubmed_id] = content.get("AuthorList", [])
                
            except Exception as e:
                print(f"Error processing {file_path}: {e}")
                continue
    
    data = list(data_dict.values())
    return pd.DataFrame(data)

# interpro_cleaner() takes in no arg.
# Function: Extracts InterPro data from ./interpro_json/*/*.json, skipping empty files.
# Extracts relevant fields, and extracts uniprot_id from file name and ensembl_id from dir name.
# go_id is the primary key here.
# Output: pandas DataFrame containing extracted InterPro data.
def interpro_cleaner():
    data = []
    ensembl_dirs = glob.glob("./interpro_json/*")
    
    for ensembl_dir in ensembl_dirs:
        if not os.path.isdir(ensembl_dir):
            continue
            
        ensembl_id = os.path.basename(ensembl_dir)
        json_files = glob.glob(os.path.join(ensembl_dir, "*.json"))
        
        for file_path in json_files:
            try:
                with open(file_path, "r") as f:
                    content = json.load(f)
                if not content or "results" not in content:
                    continue
                
                uniprot_id = os.path.splitext(os.path.basename(file_path))[0]
                
                for result in content.get("results", []):
                    metadata = result.get("metadata", {})
                    go_terms = metadata.get("go_terms")
                    
                    if go_terms and isinstance(go_terms, list):
                        for go_term in go_terms:
                            row = {
                                "ensembl_id" : ensembl_id,
                                "uniprot_id" : uniprot_id,
                                "interpro_id" : metadata.get("accession", pd.NA),
                                "interpro_name" : metadata.get("name", pd.NA),
                                "interpro_type" : metadata.get("type", pd.NA),
                                "go_id" : go_term.get("identifier", pd.NA)
                            }
                            data.append(row)
                    else:
                        row = {
                            "ensembl_id" : ensembl_id,
                            "uniprot_id" : uniprot_id,
                            "interpro_id" : metadata.get("accession", pd.NA),
                            "interpro_name" : metadata.get("name", pd.NA),
                            "interpro_type" : metadata.get("type", pd.NA),
                            "go_id" : pd.NA
                        }
                        data.append(row)
                        
            except Exception as e:
                print(f"Error processing {file_path}: {e}")
                continue
    
    return pd.DataFrame(data)

# go_cleaner() takes in no arg.
# Function: Extracts Gene Ontology(GO) data from ./go_json/*/*.json, skipping empty files.
# Filters and extracts the result matching the go_id in file name.
# Extracts relevant fields, and extracts go_id from file name and ensembl_id from dir name.
# Output: pandas DataFrame containing extracted GO data.
def go_cleaner():
    data = []
    ensembl_dirs = glob.glob("./go_json/*")
    
    for ensembl_dir in ensembl_dirs:
        if not os.path.isdir(ensembl_dir):
            continue
            
        ensembl_id = os.path.basename(ensembl_dir)
        json_files = glob.glob(os.path.join(ensembl_dir, "*.json"))
        
        for file_path in json_files:
            try:
                with open(file_path, "r") as f:
                    content = json.load(f)
                if not content or "results" not in content:
                    continue
                
                go_id_from_file = os.path.splitext(os.path.basename(file_path))[0]
                go_id_formatted = go_id_from_file.replace("_", ":", 1)
                
                for result in content.get("results", []):
                    if result.get("id") == go_id_formatted:
                        definition = result.get("definition", {})
                        row = {
                            "ensembl_id" : ensembl_id,
                            "go_id" : result.get("id", pd.NA),
                            "go_name" : result.get("name", pd.NA),
                            "go_aspect" : result.get("aspect", pd.NA),
                            "go_definition" : definition.get("text", pd.NA) if definition else pd.NA
                        }
                        data.append(row)
                        break
                        
            except Exception as e:
                print(f"Error processing {file_path}: {e}")
                continue
    
    return pd.DataFrame(data)

# final_merge() takes in 5 args: all 5 pd DataFrames
# Function: Merges all 5 dataframes into one.
# Merges on go_id as the primary key.
# Prefers NCBI gene data over Ensembl, as NCBI is secondary database while Ensembl is integrated database.
# Output: final_df which is the final pandas DataFrame ready for parsing to SQL.
def final_merge(ensembl_df, ncbi_gene_df, pubmed_df, interpro_df, go_df):
    
    final_df = interpro_df.merge(go_df, on = ["ensembl_id", "go_id"], how = "outer")
    final_df = final_df.merge(ensembl_df, on = "ensembl_id", how = "outer", suffixes = ("", "_ens"))
    final_df = final_df.merge(ncbi_gene_df, on = "ensembl_id", how = "outer", suffixes = ("", "_ncbi"))
    final_df = final_df.merge(pubmed_df, on = "ensembl_id", how = "outer")
    
    if "gene_symbol_ncbi" in final_df.columns:
        final_df["gene_symbol"] = final_df["gene_symbol_ncbi"].fillna(final_df["gene_symbol"])
        final_df.drop("gene_symbol_ncbi", axis = 1, inplace = True)
    
    if "organism_ncbi" in final_df.columns:
        final_df["organism"] = final_df["organism_ncbi"].fillna(final_df["organism"])
        final_df.drop("organism_ncbi", axis = 1, inplace = True)
    
    if "chromosome_ncbi" in final_df.columns:
        final_df["chromosome"] = final_df["chromosome_ncbi"].fillna(final_df["chromosome"])
        final_df.drop("chromosome_ncbi", axis = 1, inplace = True)
    
    if "start_ncbi" in final_df.columns:
        final_df["start"] = final_df["start_ncbi"].fillna(final_df["start"])
        final_df.drop("start_ncbi", axis = 1, inplace = True)
    
    if "end_ncbi" in final_df.columns:
        final_df["end"] = final_df["end_ncbi"].fillna(final_df["end"])
        final_df.drop("end_ncbi", axis = 1, inplace = True)
    
    column_order = [
        "ensembl_id", "gene_symbol", "organism", "taxon_id", "description", 
        "chromosome", "strand", "start", "end", 
        "entrez_id", "weight", "summary", 
        "pubmed_id", "doi", "publication_date", "author_list", 
        "uniprot_id", "interpro_id", "interpro_name", "interpro_type", 
        "go_id", "go_name", "go_aspect", "go_definition"
    ]
    final_df = final_df[column_order]
    return final_df

if __name__ == "__main__":
    main()

