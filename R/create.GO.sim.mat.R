#' @export
#' @import GOSemSim
#' @import AnnotationDbi
#' @import AnnotationHub
create.GO.sim.mat <- function(gtf.file.path, OrgDb, ont='BP', computeIC=FALSE, keyType, measure='Wang', out.file.path) {
  temp.folder <- '/tmp/dgeclustering'
  system(paste('mkdir -p', temp.folder))
  geneids.filepath <- file.path(temp.folder, 'geneids.tsv')
  path <- paste(system.file(package='DGEclustering'), 'gene_ids_from_gtf.py', sep='/')
  system(paste(path,
	       '<', gtf.file.path,
	       '>', geneids.filepath))
  genes <- unique(read.table(geneids.filepath, header=FALSE, stringsAsFactors=FALSE)$V1)

}
