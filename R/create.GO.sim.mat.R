#' @export
#' @import GOSemSim
#' @import AnnotationDbi
#' @import AnnotationHub
create.GO.sim.mat <- function(gtf.file.path=NULL, genes=NULL, OrgDb, ont='BP', keyType, computeIC=FALSE, measure='Wang', out.file.path=NULL) {
  if (!is.null(gtf.file.path)){
    temp.folder <- '/tmp/dgeclustering'
    system(paste('mkdir -p', temp.folder))
    geneids.filepath <- file.path(temp.folder, 'geneids.tsv')
    path <- paste(system.file(package='DGEclustering'), 'gene_ids_from_gtf.py', sep='/')
    system(paste(path,
                 '<', gtf.file.path,
                 '>', geneids.filepath))
    genes <- unique(read.table(geneids.filepath, header=FALSE, stringsAsFactors=FALSE)$V1)
  }
  semData <- godata(OrgDb=orgdb, ont=ont, computeIC=computeIC)
  res <- mgeneSim(genes, semData, measure=measure)
  if(!is.null(out.file.path)){
    write.table(res, out.file.path, sep='\t', header=TRUE)
  }
  else {
    return(res)
  }
}
