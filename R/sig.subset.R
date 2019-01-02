#' @title
#' @description
#' @param datasets
#' @param geneCol the column name of gene IDs
#' @param x.dsNumber
#' @param y.dsNumber
#' @param x.threshold (adj)p-value threshold of the first dataset
#' @param y.threshold (adj)p-value threshold of the second dataset
#' @param adjPvalue if TRUE, use adjusted p-value for the threshold; if FALSE, use p-value
#' @export
#' @import png
#' @import ggplot2
sig.subset <- function(datasets, geneCol, x.dsNumber=1, y.dsNumber=2, x.threshold=0.05, y.threshold=0.05, adjPvalue=TRUE, dotSize=20) {
  # export datasets to temp folder
  temp.folder <- '/tmp/dgeclustering'
  system(paste('mkdir -p', temp.folder))
  filepaths <- file.path(temp.folder, paste0(names(datasets), '.tsv'))
  for (i in 1:length(datasets)){ 
    write.table(datasets[[i]], file=filepaths[i], sep='\t', row.names=FALSE)
  }
  if (is.null(geneCol)) {
    stop('gene column not specified.')
  }
  python.boolean.convert <- function(bool) {
    if (bool == TRUE)
      return('1')
    else
      return('0')
  }
  path <- paste(system.file(package='DGEclustering'), 'sig_subset.py', sep='/')
  system(paste(path,
             '-f', paste(filepaths, collapse=' '),
             '-g', geneCol,
	     '-n1', x.dsNumber - 1,
             '-n2', y.dsNumber - 1,
             '-d', temp.folder,
             '-o', 'temp',
             '-x', x.threshold,
             '-y', y.threshold,
             '-a', python.boolean.convert(adjPvalue),
             '-z', dotSize))
  # plot the scatter plot
  img <- readPNG(file.path(temp.folder, 'temp_sig_plot.png'))
  g <- grid::rasterGrob(img, interpolate=TRUE) 
  p <- ggplot() +
    annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)
  if (length(datasets) == 2) {
    dis <- read.table(file.path(temp.folder, 'temp_disagreeing_genes.tsv'), header=TRUE, 
      check.names=FALSE, sep='\t', stringsAsFactors=FALSE, quote='')
    con <- read.table(file.path(temp.folder, 'temp_agreeing_genes.tsv'), header=TRUE, 
      check.names=FALSE, sep='\t', stringsAsFactors=FALSE, quote='')
    res <- list(dis, con, p)
    names(res) <- c('dis', 'con', 'p')
    return(res)
  } 
  else {
    dat <- read.table(file.path(temp.folder, 'temp_all_sig_genes.tsv'), header=TRUE, check.names=FALSE, 
      sep='\t', stringsAsFactors=FALSE, quote='')
    res <- list(dat, p)
    names(res) <- c('dat', 'p')
    return(res)
  }
}
