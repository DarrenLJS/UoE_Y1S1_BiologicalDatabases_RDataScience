#########################################################################
#                      Accessing BioMart                                #
#See   https://bioconductor.org/packages/release/bioc/vignettes/biomaRt #
#       /inst/doc/accessing_ensembl.html                                # 
#########################################################################

#BiomaRt is already installed on the server- just load the library

library(biomaRt)

#Note: sometimes marts are unavailable but we can use mirrrors
#https://www.ensembl.org - main site
#https://useast.ensembl.org -US East coast mirror
#https://uswest.ensembl.org -US West coast mirror

#Try UK mirror first
ensembl_host <-"https://www.ensembl.org"

#Ensembl genes to example ensembl query
test_ensembl_genes <-c("ENSMUSG00000059552","ENSMUSG00000017146")

listEnsembl()

ensembl <- useEnsembl(biomart = "genes",host=ensembl_host)

datasets <- listDatasets(ensembl)

#Here I use the grep to find a data set of interest (eg the mouse)
datasets[grep("musculus",datasets[,1]),]

ensembl <- useDataset(dataset = "mmusculus_gene_ensembl", mart = ensembl)

#Filters are the restrictions you apply to results eg only Ensembl genes
filters = listFilters(ensembl)
head(filters)

#Attributes are the information you want back from the query
attributes = listAttributes(ensembl)
head(attributes)

resultAnnot <- biomaRt::getBM(values=test_ensembl_genes,attributes = c("ensembl_gene_id","external_gene_name","chromosome_name","description"),filters="ensembl_gene_id",mart=ensembl)

#Deal with  the results
class(resultAnnot)
head(resultAnnot)




