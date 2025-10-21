library(rentrez)
# wrapper for https://www.ncbi.nlm.nih.gov/books/NBK25500/
#See https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html

#list databases
entrez_dbs()

#explore pubmed database
entrez_db_summary("pubmed")
entrez_db_searchable("pubmed")

#perform a query
result <-entrez_search(db   = "pubmed",
              term = "(Oct4[TITL]) AND (mouse[TITL])",retmode="json")

#buried in the return file is a list of IDS of the hits
result$file$esearchresult$idlist

#grab the one we want-just picked one ID from the list as an example
entrez_fetch(db="pubmed",id="35229007",rettype="text")

#or how about getting  a sequence record?
entrez_db_summary("nuccore")
entrez_db_searchable("nuccore")

query <- "Mouse[Organism] AND Pou5f1[Gene]"
nuc_search <- entrez_search(db="nuccore", query, use_history=TRUE,retmode="json")

#grap one or more record(here we just grab the first hit)
nuc_seqs <- entrez_fetch(db = "nuccore", WebEnv = NULL, id=nuc_search$file$esearchresult$idlist,rettype = "fasta", retmax = 1)

