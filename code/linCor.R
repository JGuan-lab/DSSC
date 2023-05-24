linCor <- function(yref,
                   iters = 100,
                   pval = .01,
                   n.types = NULL,
                   scree = 'cumvar',
                   logTransform = F){
  
  if (all(class(yref) != c('matrix'))){
    stop("matrix not supplied in yref")
  }
  
  if (logTransform){
    yref = log2(yref+1)
  }
  
  row.means = rowSums(yref)
  yref = yref[row.means != 0,]
  
  lo <- linseed::LinseedObject$new(yref)
  lo$calculatePairwiseLinearity()
  lo$calculateSpearmanCorrelation()
  lo$calculateSignificanceLevel(iters)
  lo$filterDatasetByPval(pval)
  # lo$svdPlot()
  if (all(lo$genes$pvals > pval)) {
    return(list(prop = matrix(rep(100,n.types * nrow(yref)),
                              ncol = n.types),
                sigs = NA))
  }
  
  if (is.null(n.types)){
    n.types = findNumberCells(yref,scree = scree)
  }
  
  lo$setCellTypeNumber(n.types)
  lo$project("full")
  # lo$projectionPlot(color="filtered")
  
  lo$project("filtered")
  lo$smartSearchCorners(dataset="filtered", error="norm")
  lo$deconvolveByEndpoints()
  # linseed::plotProportions(lo$proportions)
  
  # We can also use tSNE to haave an idea of how data looks like when dimensionally reduced.
  

  # lets select 100 genes closest to 
  # lo$selectGenes(100)
  # lo$tsnePlot()

  
  # To compare with actual proportions you can use `dotPlotProportions` function

  # data("proportionsLiverBrainLung")
  # dotPlotPropotions(lo$proportions, proportionsLiverBrainLung[, 10:42], guess=TRUE)
  
  return(list(prop = lo$proportions,
              sig = lo$signatures))
  
}
