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
#' @references \url{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-42}
#' @return
#' @examples  \dontrun{}
DGE.clust <- function(expressions, annotations=NULL, integrate.method='intego', clust.method='agnes', nb.group, OrgDb=NULL, ont='BP', keyType=NULL, alpha=1, genclust.priori=FALSE, nb.generation=500, LIM.ASSO=4, LIM.COR=0.5, nb.dim=NULL, sim.mat=NULL, random.seed=123){
  nb.dim.ex <- ncol(expressions)
  nb.dim.an <- min((nrow(annotations) - 1), (ncol(annotations) - 1))
  expressions <- scale(expressions)

  EVALUATE <- function(groups, expressions){
    #--- code modified from InteGO ---#
    INDIC = function(group.element){
    if (length(group.element) > 1){
      correlations = cor(t(expressions)[,group.element])
      res = sum(correlations[lower.tri(correlations, diag = F)]) / (length(group.element) * (length(group.element) - 1) / 2)
    } 
    else 	
      res = -1/3 # the lowest score
    names(res) = names(group.element)
    return(res)
    }   
    original.scores = c()
    for (i in 1:length(groups))
      original.scores <- c(original.scores, INDIC(groups[[i]]))
    scores <- (original.scores + 1/3) * (3/4) # rescale scores from range (-1/3 ~ 1) to range (0 ~ 1)
    ave <- round(sum(scores) / length(groups), 2)
    eva.res <- list(ave, scores, original.scores)
    names(eva.res) <- c('average', 'scores', 'original.scores')
    return(eva.res)
  }

  PVALUES <- function(groups, original.scores, expressions, nb.sim=100, abaque=NULL, threshold=0.1){
    set.seed(random.seed)
    size = function(x){return(as.numeric(lapply(1:length(x), function(g){return(length(x[[g]]))})))}
    taille.g = unique(size(groups))
    taille.g = taille.g[taille.g > 1]  
    ABAQUE = function(taille, abaque){
      simu = replicate(nb.sim, EVALUATE(groups = list(sample(row.names(expressions))[1:taille]), 
					  expressions = expressions)$original.scores)
      # simu = t(simplify2array(simu))
      abaque[[taille]] = simu ; names(abaque)[taille] = taille
      return(abaque)
    }
    for (i in 1:length(taille.g)){
      abaque = ABAQUE(taille.g[i], abaque)
    }
    abaque <- as.data.frame(do.call(cbind, abaque))
    PVAL = function(i, groups, original.scores, abaque){
      if (length(groups[[i]]) != 1){
        pval = length(which(as.numeric(abaque[[length(groups[[i]])]]) > original.scores[i])) / 
		length(abaque[[length(groups[[i]])]])
      } 
      else{ 
	pval = NA
      }
      return(pval)
    }
    pvalues = unlist(lapply(1:length(groups), PVAL, groups, original.scores, abaque))
    sig.clusters <- length(pvalues[pvalues<0.1 & !is.na(pvalues)]) / length(pvalues)
    pval.res <- list(pvalues, sig.clusters, abaque)
    names(pval.res) = names('pvalues', 'proportion sig. clusters', 'simulated sets')
    return(pval.res)
    # proportion of pvalues
  }
  #--- code modified from InteGO ---#

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
      genes <- rownames(expressions)
      if (is.null(sim.mat)){
        semData <- godata(OrgDb=OrgDb, ont=ont, keytype=keyType, computeIC=FALSE)
        GO.sim <- mgeneSim(genes, semData, measure='Wang')
      }
      else {
        GO.sim <- sim.mat[rownames(sim.mat) %in% genes, colnames(sim.mat) %in% genes]
      }
      # semantic dissimlarity may result in some genes missing
      sub.genes <- genes[genes %in% colnames(GO.sim)]
      sub.expressions <- expressions[rownames(expressions) %in% sub.genes,]
      sem.dis <- 1 - GO.sim

      # dissimilarity matrix for expression
      exp.dis <- as.matrix(dist(sub.expressions, diag=TRUE, upper=TRUE))

      # integration
      integrated.matrix <- sem.dis ^ alpha * exp.dis
      integrated.matrix <- scale(integrated.matrix)
      PCA <- PCAsimple(integrated.matrix)$ind[, 1:nb.dim]
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
    evaluation <- EVALUATE(groups, expressions)
    pvalues <- PVALUES(groups, evaluation$original.scores, expressions) 
    res <- list(groups, integrated.matrix, MCA, vignette, list(evaluation, pvalues))
    names(res) <- c('groups', 'integrated.matrix', 'MCA', 'vignette', 'evaluation')
  }
  else {
    evaluation <- EVALUATE(groups, sub.expressions)
    pvalues <- PVALUES(groups, evaluation$original.scores, sub.expressions)
    res <- list(groups, integrated.matrix, PCA, vignette, list(evaluation, pvalues))
    names(res) <- c('groups', 'integrated.matrix', 'PCA', 'vignette', 'evaluation')
  }
  return(res)
}
