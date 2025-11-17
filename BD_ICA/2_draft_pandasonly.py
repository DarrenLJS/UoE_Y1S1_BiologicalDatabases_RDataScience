#!/usr/bin/python3

import glob
import json
import pandas as pd
import os

def main():
    print("Merging data into pandas DataFrame...")
    ensembl_df = ensembl_cleaner()
    ncbi_gene_df = ncbi_gene_cleaner()
    pubmed_df = pubmed_cleaner()
    interpro_df = interpro_cleaner()
    go_df = go_cleaner()
    final_df = final_merge(ensembl_df, ncbi_gene_df, pubmed_df, interpro_df, go_df)
    print("Data merging completed! Pandas DataFrame ready for parsing to SQL!")
    print(final_df)

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
