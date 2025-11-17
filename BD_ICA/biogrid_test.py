#!/usr/bin/python3

import pandas as pd
import requests
import json

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
base_url = "https://webservice.thebiogrid.org/interactions"
params = {
    "accessKey" : "cabf28e59b02c44454c83ce11423f27c", 
    "format" : "json", 
    "geneList" : "|".join(ens_id), 
    "searchIds" : "true", 
    "includeInteractors" : "true", 
    "includeHeader" : "true", 
    "additionalIdentifierTypes" : "ENSEMBL"
}
headers = {"Accept" : "application/json"}

response = requests.get(f"{base_url}", params = params)
lookup_data = response.json()

with open("biogrid.json", "w") as file:
    json.dump(lookup_data, file, indent = 2)

