#' @title
#' @description
#' @param expressions gene expression dataset
#' @param annotations GO term annotation dataset
#' @param clust.methods one of 'intego', 'genclust', 'ward' (for clustering without annotations).
#' @param nb.group number of clustering groups
#' @param genclust.priori If TRUE, use intego result as a priori. Default is FALSE. 
#' @param nb.generation number of generations for genclust
#' @param LIM.ASSO
#' @param LIM.COR
#' @export
#' @import InteGO
#' @import rlist
#' @import GOSemSim
#' @import AnnotationDbi
#' @import scales
#' @references \url{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-42}
#' @return
#' @examples  \dontrun{}
DGE.clust <- function(expressions, annotations=NULL, integrate.method='intego', clust.method='agnes', nb.group, OrgDb=NULL, ont='BP', keyType=NULL, genclust.priori=FALSE, nb.generation=500, LIM.ASSO=4, LIM.COR=0.5, nb.dim=NULL){
  nb.dim.ex <- ncol(expressions)
  nb.dim.an <- min((nrow(annotations) - 1), (ncol(annotations) - 1))
  expressions <- scale(expressions)

  evaluate <- function(groups, expressions, annotations){
    options(warn=-1) # turn off warnings from unused evaluate indicator
    eva <- Indicators(groups, expressions, annotations)
    options(warn=0) # turn warnings back on
    original.scores <- c()
    for(i in 1:length(groups)){
      if (is.na(eva[[i]][[1]]))
        original.scores <- c(original.scores, -1/3)
      else
        original.scores <- c(original.scores, eva[[i]][[1]])
    }
    scores <- rescale(original.scores, to=c(0,1))
    ave <- round(sum(scores) / length(groups), 2)
    eva.res <- list(ave, scores, original.scores)
    names(eva.res) <- c('average', 'scores', 'original.scores')
    return(eva.res)
  }

  if (!is.null(annotations)){
    if (integrate.method == 'intego'){
      integrated.matrix <- Integration(annotations, expressions, nb.dim.ex, LIM.ASSO, LIM.COR)
      integrated.matrix <- apply(integrated.matrix, 2, as.factor)
      if (is.null(nb.dim)){
        MCA <- MCAsimple(integrated.matrix)[, 1:nb.dim.an]
      }
      else {
        MCA <- MCAsimple(integrated.matrix)[, 1:nb.dim]
      }
      DIST <- dist(MCA, diag=TRUE, upper=TRUE)
    }
    else{
      if (is.null(OrgDb) && is.null(keyType)){
        stop('the argument OrgDb is required for new.distance integration method.')
      }
      semData <- godata(OrgDb=OrgDb, ont='BP', computeIC=FALSE)
      genes <- mapIds(OrgDb, keys=rownames(expressions), column='ENTREZID', keytype=keyType, multiVals='first')
      genes <- genes[!is.na(genes)] 
      GO.sim <- mgeneSim(genes, semData, measure='Wang')
      sub.genes <- genes[genes %in% colnames(GO.sim)]
      sub.expressions <- expressions[rownames(expressions) %in% names(sub.genes),]
      sub.annotations <- annotations[rownames(annotations) %in% names(sub.genes),]
      exp.sim <- as.matrix(dist(sub.expressions, diag=TRUE, upper=TRUE))
      integrated.matrix <- exp.sim * GO.sim
      PCA <- PCAsimple(integrated.matrix)$ind[, 1:2]
      DIST <- dist(PCA, diag=TRUE, upper=TRUE)
    }
  
    if (clust.method != 'genclust'){ # set intego as default
      groups <- clustering(DIST, mode='Classification', nb.group=nb.group)
    }

    else{
      # input for gensclust does not need to be scaled (i.e. MCA is already scaled)
      # generate MCA input for Genclust (i.e. input is GO term-integrated)
      line1 <-
        data.frame(as.list(c(
          nrow(MCA), ncol(MCA), 10, rep(NA, ncol(MCA) - 2)
        )))
      colnames(line1) <- seq(1, ncol(line1))

      line2 <- list('gene')
      for (i in seq(1, ncol(MCA))) {
        line2[i + 1] <- paste('f', i, sep = '')
      }
      line2 <- data.frame(line2)
      colnames(line2) <- colnames(line1)

      gen.input <- NULL
      gen.input <- cbind(as.data.frame(rownames(MCA)), as.data.frame(MCA))
      colnames(gen.input) <- colnames(line1)

      package.dir <- system.file(package="DGEclustering")
      gen.input <- rbind(line1, line2, gen.input)
      write.table(
        gen.input,
        file = paste(package.dir, 'genclust_sig_data.tsv', sep='/'),
        sep = '\t',
        col.names = FALSE,
        row.names = FALSE,
        na = '',
        quote = FALSE
      )

      # generate initialization file
      filepath <- paste(package.dir, 'out.tmp', sep='/')
      fileConn <- file(filepath)
      write('#Comment\n#Comment',
            filepath,
            append = FALSE,
            sep = '\n')
      if (genclust.priori==FALSE){
          set.seed(123)
          ran <- split(seq(1, nrow(MCA)), sample(1:nb.group, nrow(MCA), replace = TRUE))
          for (i in 1:length(ran)) {
            write(paste(length(ran[[i]]), paste(ran[[i]], collapse = ' '), sep = '\t'),
                  filepath,
                  append = TRUE,
                  sep = '\n')
          }
      }
      else {
        DIST <- dist(MCA, diag = TRUE, upper = TRUE)
        groups <- clustering(DIST, mode='Classification', nb.group=nb.group)
        rownames(gen.input) <- seq(1, nrow(gen.input))
        for (i in 1:nb.group){
          temp <- c()
          for (j in 1:length(groups[[i]])){
            temp <- c(temp, as.numeric(rownames(gen.input)[gen.input[, 1] == groups[[i]][j]]) - 2)
          }
          write(paste(length(groups[[i]]), paste(temp, collapse = ' '), sep = '\t'),
                filepath,
                append = TRUE,
                sep = '\n')
        }
      }
      close(fileConn)

      # run genclust
      system(
        paste('genclust',
              paste(package.dir, 'genclust_sig_data.tsv', sep='/'),
              nb.group,
              nb.generation,
              paste(package.dir, 'genclust_out.txt', sep='/'),
              0,
              0
              ), ignore.stdout = TRUE)

      # import genclust result
      gen.out <-
        read.table(paste(package.dir, 'genclust_out.txt', sep='/'),
                   header = FALSE,
                   sep = '\t')
      gen.out[] <- lapply(gen.out, as.character)
  
      x <- list()
      g <- 0
      for (row in 1:(nrow(gen.out) - 1)) {
        string.list <- unlist(strsplit(gen.out[row, 1], ' '))
        if (string.list[1] == 'CLUSTER') {
          if (g > 0) {
            x <- list.append(x, temp)
          }
          g <- g + 1
          temp <- c()
          i <- 1
        }
        else{
          temp[i] <- string.list[2]
          i <- i + 1
        }
      }
      groups <- list.append(x, temp)
      names(groups) <- paste('Group', 1:g, sep = '.')
    }
    # create log
    vignette <- strsplit(
            paste(paste('## Expressions input:', dim(expressions)[2], 'cols.', paste(colnames(expressions), collapse=' ')),
                  paste('## Annotations input:', dim(annotations)[2], 'cols.', paste(colnames(annotations), collapse=' ')),
                  paste('## Number of groups:', nb.group), sep='\n'), '\n')
  }

  else { # annotations = NULL
    DIST <- dist(expressions, diag=TRUE, upper=TRUE)
    groups <- clustering(DIST, mode='Classification', nb.group=nb.group)
    integrated.matrix <- NULL
    MCA <- NULL
    # create log
    vignette <- strsplit(
            paste(paste('## Annotations input:', dim(annotations)[2], 'cols.', paste(colnames(annotations), collapse=' ')),
                  paste('## Number of groups:', nb.group), sep='\n'), '\n')
  }
  
  if (integrate.method == 'intego'){
    evaluation <- evaluate(groups, expressions, annotations)
    res <- list(groups, integrated.matrix, MCA, vignette, evaluation)
    names(res) <- c('groups', 'integrated.matrix', 'MCA', 'vignette', 'evaluation')
  }
  else {
    evaluation <- evaluate(groups, sub.expressions, sub.annotations)
    res <- list(groups, integrated.matrix, PCA, vignette, evaluation)
    names(res) <- c('groups', 'integrated.matrix', 'PCA', 'vignette', 'evaluation')
  }
  return(res)
}
