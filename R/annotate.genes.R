#' @title annotate genes with GO terms
#'
#' @description Generate a GO term binary matrix from a paired DGE dataset
#' @param OrgDb object returned by AnnotationHub.
#' @param keyType type of genes. e.g. 'ENSEMBL', 'ENTREZID', 'SYMBOL'.
#' @param genes a vector of genes.
#' @param GOEnrichment if TRUE, annotate genes with GO enrichment. Default FALSE.
#' @param BgGenes background genes if use GO enrichment.
#' @param ont One of 'ALL', 'BP', 'BF', or 'CC'.
#' @export
#' @import AnnotationHub
#' @import clusterProfiler
#' @return ann a dataframe of the GO term binary matrix
#' @examples  \dontrun{}
annotate.genes <- function(OrgDb, keyType, genes, GOEnrichment=FALSE, BgGenes=NULL, ont='ALL'){
  ds.keys <- sapply(genes, as.character)
  unv.keys <- sapply(BgGenes, as.character)

  if (GOEnrichment == FALSE){
    # assign GO terms directly
    GO.res <- select(OrgDb, ds.keys, c(keyType, 'GOALL', 'ONTOLOGYALL'), keyType)
    GO.res <- GO.res[complete.cases(GO.res['GOALL']),]
    if (ont != 'ALL'){
      GO.res[GO.res['ONTOLOGYALL'] == ont,]
    }
  }
  else {
    # GO enrichment
    GO.res <- enrichGO(gene=ds.keys, OrgDb=OrgDb, keyType=keyType, ont=ont, universe=unv.keys)
  }

  # construct binary GO term matrix
  ann <- data.frame(genes=ds.keys)
  rownames(ann) <- ds.keys
  if (GOEnrichment == FALSE){
    for (row in 1:nrow(GO.res)){
      GO.term <- GO.res[row, 'GO']
      if (!GO.term %in% colnames(ann))
        ann[GO.term] <- 0
      ann[rownames(ann)==GO.res[row, keyType], GO.term] <- 1
    }
  }
  else {
    for (row in 1:nrow(GO.res)){
      GO.term <- as.character(GO.res[row, 'Description'])
      if (is.na(GO.term))
	next
      if (!GO.term %in% colnames(ann))
        ann[GO.term] <- 0
      geneIDs <- unlist(strsplit(GO.res[row, 'geneID'], '/'))
      ann[ds.keys %in% geneIDs, GO.term] <- 1
    }
  }
  ann <- ann[,-1]
  return(ann)
}
