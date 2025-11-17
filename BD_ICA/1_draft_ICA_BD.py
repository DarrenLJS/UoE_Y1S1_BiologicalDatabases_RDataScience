#!/usr/bin/bash

from pybiomart import Server
from pybiomart import Dataset
import entrezpy.esearch.esearcher
import entrezpy.efetch.efetcher
import entrezpy.log.logger

from Bio import Entrez
import time
import pandas as pd

def search_and_get_gene_info(ens_df):
    Entrez.email = "s2906787@bioinfmsc9.bio.ed.ac.uk"
    results = []
    
    for i, row in ens_df.iterrows():
        term = f"{row['gene_name']}[Gene Name] AND {row['species']}[Organism]"
        
        try:
            # Step 1: Search for IDs
            search_handle = Entrez.esearch(db="gene", term=term, retmax=1)
            search_record = Entrez.read(search_handle)
            search_handle.close()
            
            if search_record['IdList']:
                gene_id = search_record['IdList'][0]  # Get first result
                
                time.sleep(0.34)
                
                # Step 2: Get summary info
                summary_handle = Entrez.esummary(db="gene", id=gene_id)
                summary = Entrez.read(summary_handle)
                summary_handle.close()
                
                doc = summary['DocumentSummarySet']['DocumentSummary'][0]
                
                results.append({
                    'query_gene': row['gene_name'],
                    'query_species': row['species'],
                    'gene_id': gene_id,
                    'symbol': doc['Name'],
                    'description': doc['Description'],
                    'chromosome': doc.get('Chromosome', 'N/A'),
                    'map_location': doc.get('MapLocation', 'N/A'),
                    'organism': doc['Organism']['ScientificName']
                })
            else:
                results.append({
                    'query_gene': row['gene_name'],
                    'query_species': row['species'],
                    'gene_id': 'Not found',
                    'symbol': 'N/A',
                    'description': 'N/A',
                    'chromosome': 'N/A',
                    'map_location': 'N/A',
                    'organism': 'N/A'
                })
                
        except Exception as e:
            print(f"Error on row {i}: {e}")
            
        time.sleep(0.34)
    
    return pd.DataFrame(results)

# Use it
results_df = search_and_get_gene_info(ens_df)
print(results_df)

def map_ensembl_pybiomart(ens_df):
    species_df = ens_df.loc[:, ["ensembl_id", "species"]]

    dataset_name = []
    for i, row in species_df.iterrows():
        if pd.isna(row["species"]):
            print(f"Dataset cannot be found for {row['ensembl_id']}.")
            dataset_name.append(None)
            continue
        prefix_list = row["species"].split("_")
        dataset = f"{prefix_list[0][0].lower()}{prefix_list[-1].lower()}_gene_ensembl"
        dataset_name.append(dataset)
    
    species_df["biomart_dataset"] = dataset_name
    species_df = species_df.loc[:, ["ensembl_id", "biomart_dataset"]]

    server = Server(host = "http://www.ensembl.org")
    mart = server["ENSEMBL_MART_ENSEMBL"]
    
    all_results = []
    for i, row in species_df.iterrows():
        print(f"Searching BioMart for {row['ensembl_id']}...")
        if pd.isna(row["biomart_dataset"]):
            print(f"No result found for {row['ensembl_id']}.")
            continue
        dataset = Dataset(name = row["biomart_dataset"], host = "http://www.ensembl.org")
        result = dataset.query(
            attributes = [
                "ensembl_gene_id",
                "entrezgene_id",
                "uniprotswissprot"
            ],
            filters = {"link_ensembl_gene_id" : row["ensembl_id"]}
        )
        clean_result = result.dropna()
        clean_result = clean_result.rename(columns = {
            "Gene stable ID": "ensembl_id",
            "NCBI gene (formerly Entrezgene) ID": "entrez_id",
            "UniProtKB/Swiss-Prot ID": "uniprot_id"
        })
        print(clean_result)
        all_results.append(clean_result)
        print(f"BioMart result retrieved for {row['ensembl_id']}!")
    
    results_df = pd.concat(all_results, ignore_index = True)
    print(results_df)
    print(species_df)
    xrefs_df = results_df.merge(species_df, on = "ensembl_id", how = "outer", sort = False)
    return xrefs_df

def access_ensembl_pybiomart(id_list):
    server = Server(host = 'http://www.ensembl.org')
    mart = server['ENSEMBL_MART_ENSEMBL']
    mart.list_datasets()
    
    def get_dataset_organism(ids):
        prefix_dataset = {
            "ENSG" : ["hsapiens_gene_ensembl", "Homo sapiens"], 
            "ENSMUSG" : ["mmusculus_gene_ensembl", "Mus musculus"], 
            "ENSDARG" : ["drerio_gene_ensembl", "Danio rerio"], 
            "ENSXETG" : ["xtropicalis_gene_ensembl", "Xenopus tropicalis"], 
            "ENSGALG" : ["ggallus_gene_ensembl", "Gallus gallus"], 
            "ENSBTAG" : ["btaurus_gene_ensembl", "Bos taurus"], 
            "ENSRNOG" : ["rnorvegicus_gene_ensembl", "Rattus norvegicus"], 
            "ENSOARG" : ["oaries_gene_ensembl" "Ovis aries"], 
            "ENSCSAG" : ["cfamiliaris_gene_ensembl", "Canis familiaris"]
        }
        for prefix, dataset in prefix_dataset.items():
            if ids.startswith(prefix):
                return dataset[0], dataset[1]
        return None
    
    all_results = []
    for i in id_list:
        dt_name, organism_name = get_dataset_organism(i)
        if not dt_name:
            print(f"Dataset cannot be found for {i}!")
            continue
        
        print(f"Searching for {i}...")
        dataset = Dataset(name = dt_name, host = 'http://www.ensembl.org')
        result = dataset.query(
            attributes = [
                'ensembl_gene_id', 
                'external_gene_name', 
                'chromosome_name', 
                'description', 
                'start_position', 
                'end_position'
            ],
            filters = {'link_ensembl_gene_id' : i}
        )
        result["Organism"] = organism_name
        print(result)

        if result.empty:
            print("No result found for {i}!")
            continue
        all_results.append(result)

    if len(all_results) < 1:
        print("No results found!")
        return pd.DataFrame()
    merged_df = pd.concat(all_results, ignore_index = True)
    print(merged_df.head(10))
    return merged_df

def access_ncbi_entrezpy(name_list):
    
    entrezpy.log.logger.set_level('WARN')
    
    for i in name_list:
        e = entrezpy.esearch.esearcher.Esearcher(
            "entrezpy", 
            "s2906787@bioinfmsc9.bio.ed.ac.uk", 
            apikey=None, 
            apikey_var=None, 
            threads=None, 
            qid=None
        )
        searcher = e.inquire({
            'db' : 'nucleotide',
            'term' : i,
            'retmax' : '20',
            'rettype' : 'uilist'
        })
        print(searcher.result.count, searcher.result.uids)

def access_ensembl_rest(ens_id_list):
    base_url = "https://rest.ensembl.org/"
    headers = {"Content-Type" : "application/json"}
    results = []

    for ens_id in ens_id_list:
        print(f"Searching Ensembl for {ens_id}...")
        response_lookup = requests.get(f"{base_url}lookup/id/{ens_id}", headers = headers)

        if response_lookup.ok:
            data = response_lookup.json()
            print(f"\n{json.dumps(data, indent = 2)}\n")
            filter_data = {
                "ensembl_id" : data.get("id"),
                "gene_name" : data.get("display_name"),
                "species" : data.get("species"),
                "description" : data.get("description"),
                "chromosome" : data.get("seq_region_name"),
                "strand" : data.get("strand"),
                "start" : data.get("start"),
                "end" : data.get("end")
            }

            response_xids = requests.get(f"{base_url}xrefs/id/{ens_id}", headers = headers)
            xids = response_xids.json()
            print(f"\n{json.dumps(xids, indent = 2)}\n")
            results.append(filter_data)
            print(f"Ensembl data retrieved for {ens_id}!")

        else:
            print(f"Failed to request Ensembl data for {ens_id}. Error: {response_lookup.status_code}")
            results.append({"ensembl_id" : ens_id})

    return pd.DataFrame(results)

def access_interpro_api(ens_id_list):
    base_url = "https://www.ebi.ac.uk/interpro/api/"
    query = "entry/interpro/Ensembl:"
    results = []

    for ens_id in ens_id_list:
        print(f"Searching InterPro for {ens_id}...")
        response = requests.get(f"{base_url}{query}{ens_id}")

        if response.ok:
            data = response.json()
            #filter_data =

def mapping_ensembl_uniprot(ens_id_list):
    base_url = "https://rest.ensembl.org/"
    response_post = requests.post(f"{base_url}idmapping/run", data = {
        "from" : "Ensembl", 
        "to" : "UniProtKB", 
        "ids" : ",".join(ens_id_list)
    })

    if response_post.ok:
        job_id = response_post.json()["jobId"]
        job = f"{base_url}idmapping/details/{job_id}"
        link = session.get(job)
        
        if link.ok:
            search = link.json()["redirectURL"]
            print(json.dumps(search, indent = 2))

    else:
        print(f"Failed to map Ensembl ids to UniProtKB ids. Error: {response_post.status_code}")

def map_ensembl_uniprot(ens_df):
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    uniprot_ids = []

    for idx, row in ens_df.iterrows():
        transcript_id = row["canonical_transcript"]
        if not transcript_id:
            uniprot_ids.append(None)
            continue

        print(f"Mapping canonical transcript {transcript_id} → UniProt...")
        params = {
            "query": f"xref:Ensembl-{transcript_id}",
            "format": "json"
        }
        response = requests.get(base_url, params=params)

        if response.ok:
            data = response.json()
            # Prefer canonical or reviewed entry
            canonical_entry = None
            reviewed_entry = None

            for result in data.get("results", []):
                if result.get("isCanonical"):
                    canonical_entry = result.get("primaryAccession")
                    break
                if result.get("entryType") == "Swiss-Prot":
                    reviewed_entry = result.get("primaryAccession")

            # Priority: canonical → reviewed → first entry
            if canonical_entry:
                uniprot_ids.append(canonical_entry)
            elif reviewed_entry:
                uniprot_ids.append(reviewed_entry)
            elif data.get("results"):
                uniprot_ids.append(data["results"][0].get("primaryAccession"))
            else:
                uniprot_ids.append(None)
        else:
            print(f"Failed UniProt request for {transcript_id}. Error: {response.status_code}")
            uniprot_ids.append(None)

    ens_df["uniprot_id"] = uniprot_ids
    return ens_df


