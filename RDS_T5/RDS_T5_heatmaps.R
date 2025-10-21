data(airway, package="airway")

library(SummarizedExperiment)
assays(airway)
expression <- assays(airway)[[1]]
head(expression)
nrow(expression)

variances <- apply(expression, 1, var)
expression <- expression[order(variances, decreasing=TRUE)[1:50],] 
head(expression)

library(org.Hs.eg.db) 
annotation <- merge(toTable(org.Hs.egSYMBOL), toTable(org.Hs.egENSEMBL2EG)) 
gene_names <- annotation$symbol[match(rownames(expression), annotation$ensembl_id)] 
rownames(expression)[! is.na(gene_names)] <- paste(gene_names[!is.na(gene_names)], rownames(expression)[! is.na(gene_names)], sep = ' / ')

experiment <- data.frame(colData(airway)) 
head(experiment) 

library(pheatmap)
expression <- log2(expression + 1)

pheatmap(
  expression,
  fontsize_row = 6,
  fontsize_col=12
)

pheatmap(
  expression,
  scale='row',
  fontsize_row = 6,
  fontsize_col=12
)

pheatmap(
  expression,
  scale='row',
  annotation_col=experiment,
  fontsize_row = 6,
  fontsize_col=12
)

pheatmap(
  expression,
  scale='row',
  annotation_col=experiment[, c('cell', 'dex')],
  fontsize_row = 6,
  fontsize_col=12
)

pheatmap(
  expression,
  scale='row',
  annotation_col=experiment[, c('cell', 'dex')],
  cluster_col=FALSE,
  fontsize_row = 6,
  fontsize_col=12
)

experiment <- experiment[with(experiment, order(dex, cell)), ] 
expression <- expression[,rownames(experiment)] 

pheatmap(
  expression,
  scale='row',
  annotation_col=experiment[, c('cell', 'dex')],
  cluster_col=FALSE,
  fontsize_row = 6,
  fontsize_col=12
)

pheatmap( 
  expression, 
  scale='row', 
  annotation_col=experiment[,c('cell', 'dex')], 
  cluster_col=FALSE, 
  annotation_colors = list('cell' = c('N052611' = 'red', 'N061011' = 'blue', 'N080611' = 'green', 'N61311' = 'yellow')), 
  fontsize_row = 6, 
  fontsize_col = 12 
)

pheatmap( 
  expression, 
  scale='row', 
  annotation_col=experiment[,c('cell', 'dex')], 
  cluster_col=FALSE, 
  annotation_colors = list('cell' = c('N052611' = 'red', 'N061011' = 'blue', 'N080611' = 'green', 'N61311' = 'yellow')), 
  main = 'My awesome plot', 
  legend = FALSE, 
  fontsize_row = 6, 
  fontsize_col = 12 
)