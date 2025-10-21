#!/usr/bin/python3

from pybiomart import Dataset
from pybiomart import Server

#Connect to the Ensembl server & print Marts
server = Server(host='http://www.ensembl.org')
print(server.list_marts())

#list the datasets  from the server
mart = server['ENSEMBL_MART_ENSEMBL']
mart.list_datasets()

#grab the required dataset
dataset = Dataset(name='mmusculus_gene_ensembl',host='http://www.ensembl.org')

#print attributes
print(dataset.list_attributes())

#print filters
print(dataset.list_filters())

#generate query results- from chosen filters and attributes
result =dataset.query(attributes=['ensembl_gene_id','external_gene_name','chromosome_name','description'],filters={'chromosome_name': ['1','2']})

#print query
print(result)
