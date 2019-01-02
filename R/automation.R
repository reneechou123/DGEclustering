#' @title automated plotting of diagnostic plots
#' @description From a starting directory, recursively search for DE output files and their corresponding R markdown files, 
#'   store them in the SQLite database, then plot the diagnostic plots for the paired datasets
#' @param rootDir path to the starting directory
#' @param geneCol the column name of gene IDs
#' @param x.threshold (adj) p-value threshold of the first dataset
#' @param y.threshold (adj) p-value threshold of the second dataset
#' @param adjPvalue if TRUE, use adjusted p-value for the threshold; if FALSE, use p-value
#' @param qqPlot if TRUE, plot Q-Q plots
#' @param fishPlot if TRUE, plot fish plots
#' @param scatterPlot if TRUE, plot scatter plots
#' @export
#' @import RSQLite
automation <- function(rootDir, geneCol, x.threshold=0.05, y.threshold=0.05, adjPvalue=TRUE, qqPlot=TRUE, fishPlot=TRUE, scatterPlot=TRUE, dotSize=20){
  
  python.boolean.convert <- function(bool) {
    if (bool == TRUE)
      return('1')
    else
      return('0')
  }
  path <- paste(system.file(package='DGEclustering'), 'automation.py', sep='/')
  system(paste(path,
               '-d', rootDir,
               '-g', geneCol,
               '-x', x.threshold,
               '-y', y.threshold,
               '-a', python.boolean.convert(adjPvalue),
               '-q', python.boolean.convert(qqPlot),
               '-f', python.boolean.convert(fishPlot),
               '-s', python.boolean.convert(scatterPlot),
               '-z', dotSize))
}

