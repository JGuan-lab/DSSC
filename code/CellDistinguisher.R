my_gecd_CellDistinguisher <- function (exprLinear, genesymb = NULL, numCellClasses = 2, minDistinguisherAlternatives = 1, 
          maxDistinguisherAlternatives = 5, minAlternativesLengthsNormalized = 0.7, 
          expressionQuantileForScale = 0.75, expressionQuantileForFilter = 1, 
          expressionConcentrationRatio = 0.333, probesWithGenesOnly = FALSE, 
          verbose = 0) 
{
  ptm <<- proc.time()
  expressionPreprocessing <- function(exprLinear, genesymb, 
                                      expressionQuantileForScale, expressionQuantileForFilter, 
                                      expressionConcentrationRatio, probesWithGenesOnly, verbose) {
    if (length(unique(rownames(exprLinear))) != length(rownames(exprLinear))) {
      mesg <- paste("There are only", length(unique(rownames(exprLinear))), 
                    "unique row names among the", length(rownames(exprLinear)), 
                    "rows of exprLinear; there are duplicates.\n")
      stop(mesg)
    }
    if (probesWithGenesOnly) {
      if (is.null(genesymb)) {
        mesg <- "probesWithGenesOnly=TRUE and genesymb=NULL filters out all probes.\n"
        stop(mesg)
      }
      ProbesWithoutGenes <- (genesymb == "")
      if ((verbose > 3) & (sum(ProbesWithoutGenes) > 0)) {
        cat("Eliminating probes without gene symbols:\n")
        print(rownames(exprLinear)[ProbesWithoutGenes])
      }
      exprLinear <- exprLinear[!ProbesWithoutGenes, ]
      genesymb <- genesymb[!ProbesWithoutGenes]
    }
    cNaN <- sum(is.nan(exprLinear))
    cNA <- sum(is.na(exprLinear)) - cNaN
    if (cNA + cNaN > 0) {
      mesg <- "The exprLinear matrix includes prohibited"
      if (cNA > 0) {
        mesg <- paste(mesg, "NA")
      }
      if (cNA * cNaN > 0) {
        mesg <- paste(mesg, "and")
      }
      if (cNaN > 0) {
        mesg <- paste(mesg, "NaN")
      }
      mesg <- paste(mesg, "values.\n")
      stop(mesg)
    }
    exprLinear <- t(t(exprLinear)/colSums(exprLinear))
    if (verbose > 3) {
      print(proc.time() - ptm)
      ptm <<- proc.time()
      cat("Columns (samples) normalized\n")
    }
    if (expressionQuantileForFilter < 1) {
      MyExpressionLimit <- quantile(exprLinear, p = expressionQuantileForFilter)
      ProbesWithExtremeExpression <- apply(exprLinear, 
                                           1, function(probe) {
                                             return(sum(probe > MyExpressionLimit) > 0)
                                           })
      if ((verbose > 3) & (sum(ProbesWithExtremeExpression) > 
                           0)) {
        cat("Eliminating probes with outlier expression:\n")
        print(rownames(exprLinear)[ProbesWithExtremeExpression])
      }
      exprLinear <- exprLinear[!ProbesWithExtremeExpression, 
      ]
      if (!is.null(genesymb)) {
        genesymb <- genesymb[!ProbesWithExtremeExpression]
      }
    }
    if (expressionConcentrationRatio > 0) {
      ProbesWithConcentratedExpression <- apply(exprLinear, 
                                                1, function(probe) {
                                                  top <- sort(probe, decreasing = TRUE)[1:2]
                                                  return(top[2] <= expressionConcentrationRatio * 
                                                           top[1])
                                                })
      if ((verbose > 3) & (sum(ProbesWithConcentratedExpression) > 
                           0)) {
        cat("Eliminating probes with concentrated expression:\n")
        print(rownames(exprLinear)[ProbesWithConcentratedExpression])
      }
      exprLinear <- exprLinear[!ProbesWithConcentratedExpression, 
      ]
      if (!is.null(genesymb)) {
        genesymb <- genesymb[!ProbesWithConcentratedExpression]
      }
    }
    matrixRank <- rankMatrix(exprLinear)[1]
    if (matrixRank < numCellClasses) {
      mesg <- sprintf("The expression matrix is of low rank.  Ask for numCellClasses <= %d.\n", 
                      matrixRank)
      stop(mesg)
    }
    if (verbose > 3) {
      cat("gecd_CellDistinguisher: rank test complete\n")
    }
    aDilution <- quantile(exprLinear[exprLinear > 1e-12], 
                          probs = c(0, expressionQuantileForScale), names = FALSE, 
                          type = 6)
    dilution <- aDilution[2] - aDilution[1]
    exprLinear <- exprLinear + dilution
    exprLinear <- t(t(exprLinear)/colSums(exprLinear))
    if (verbose > 2) {
      print(proc.time() - ptm)
      ptm <<- proc.time()
      cat(sprintf("Expression diluted with Expression[%f quantile] = %g\n", 
                  expressionQuantileForScale, dilution))
    }
    if (verbose > 3) {
      print(proc.time() - ptm)
      ptm <<- proc.time()
      cat("Q created, implicitly\n")
    }
    exprLinearBar <- exprLinear/(exprLinear %*% colSums(exprLinear))[, 
                                                                     1]
    if (verbose > 3) {
      print(proc.time() - ptm)
      ptm <<- proc.time()
      cat("Qbar created, implicitly; i.e., rows of Q normalized\n")
    }
    if (TRUE) {
      exprLinearBarAdj <- exprLinearBar
    }else {
      exprLinearBarAdj <- t(t(exprLinearBar) - colMeans(exprLinearBar))
      if (verbose > 3) {
        print(proc.time() - ptm)
        ptm <<- proc.time()
        cat("Qbar centered, implicitly\n")
      }
    }
    tExprLinear_exprLinear <- gecd_MatrixChainMultiplication(t(exprLinear), 
                                                             exprLinear, verbose = verbose - 3)
    return(list(exprLinear = exprLinear, tExprLinear_exprLinear = tExprLinear_exprLinear, 
                exprLinearBar = exprLinearBar, exprLinearBarAdj = exprLinearBarAdj, 
                genesymb = genesymb))
  }
  findCandidateDistinguishers <- function(exprLinear, tExprLinear_exprLinear, 
                                          exprLinearBarAdj, numCellClasses) {
    selectDistantMarker <- function(exprLinearBarAdj, exprLinear, 
                                    tExprLinear_exprLinear, passOneDistinguishers, allLengths) {
      lengths2 <- gecd_MatrixChainMultiplication(exprLinearBarAdj, 
                                                 tExprLinear_exprLinear, t(exprLinearBarAdj), 
                                                 diagonalOnly = TRUE, verbose = verbose - 3)
      distinguisher <- which.max(lengths2)
      passOneDistinguishers <- c(passOneDistinguishers, 
                                 distinguisher)
      allLengths <- c(allLengths, sqrt(lengths2[distinguisher]))
      if (verbose > 2) {
        iDistinguisher <- length(passOneDistinguishers)
        print(proc.time() - ptm)
        ptm <<- proc.time()
        cat(sprintf("First pass: 1 CellClass[%d] distinguisher found = %d \"%s\" at distance %g\n", 
                    iDistinguisher, distinguisher, rownames(exprLinear)[distinguisher], 
                    allLengths[iDistinguisher]))
      }
      return(list(passOneDistinguishers = passOneDistinguishers, 
                  allLengths = allLengths))
    }
    passOneDistinguishers <- NULL
    allLengths <- NULL
    if (1 <= numCellClasses) {
      iDistinguisher <- 1
      ans <- selectDistantMarker(exprLinearBarAdj, exprLinear, 
                                 tExprLinear_exprLinear, passOneDistinguishers, 
                                 allLengths)
      passOneDistinguishers <- ans$passOneDistinguishers
      allLengths <- ans$allLengths
      rm(ans)
    }
    if (2 <= numCellClasses) {
      iDistinguisher <- 2
      exprLinearBarAdj <- t(t(exprLinearBarAdj) - exprLinearBarAdj[passOneDistinguishers[iDistinguisher - 
                                                                                           1], ])
      ans <- selectDistantMarker(exprLinearBarAdj, exprLinear, 
                                 tExprLinear_exprLinear, passOneDistinguishers, 
                                 allLengths)
      passOneDistinguishers <- ans$passOneDistinguishers
      allLengths <- ans$allLengths
      rm(ans)
    }
    if (3 <= numCellClasses) {
      for (iDistinguisher in 3:numCellClasses) {
        project <- exprLinearBarAdj[passOneDistinguishers[iDistinguisher - 
                                                            1], , drop = FALSE]
        tExprLinear_exprLinear_tProject <- gecd_MatrixChainMultiplication(tExprLinear_exprLinear, 
                                                                          t(project), verbose = verbose - 3)
        projectionNorm <- solve(gecd_MatrixChainMultiplication(project, 
                                                               tExprLinear_exprLinear_tProject, verbose = verbose - 
                                                                 3))
        exprLinearBarAdj <- exprLinearBarAdj - gecd_MatrixChainMultiplication(exprLinearBarAdj, 
                                                                              tExprLinear_exprLinear_tProject, projectionNorm, 
                                                                              project, verbose = verbose - 3)
        ans <- selectDistantMarker(exprLinearBarAdj, 
                                   exprLinear, tExprLinear_exprLinear, passOneDistinguishers, 
                                   allLengths)
        passOneDistinguishers <- ans$passOneDistinguishers
        allLengths <- ans$allLengths
        rm(ans)
      }
    }
    return(list(passOneDistinguishers = passOneDistinguishers, 
                allLengths = allLengths, exprLinearBarAdj = exprLinearBarAdj))
  }
  adjustDistinguishers <- function(passOneDistinguishers, exprLinearBarAdj, 
                                   exprLinear, tExprLinear_exprLinear, maxDistinguisherAlternatives, 
                                   noDistinguishersYet) {
    projectOut <- function(passOneDistinguishers, exprLinearBarAdj, 
                           exprLinear, tExprLinear_exprLinear, noDistinguishersYet) {
      if (length(passOneDistinguishers) > 0 && noDistinguishersYet) {
        exprLinearBarAdj <- t(t(exprLinearBarAdj) - exprLinearBarAdj[passOneDistinguishers[1], 
        ])
        passOneDistinguishers <- passOneDistinguishers[-1]
      }
      if (length(passOneDistinguishers) > 0) {
        project <- exprLinearBarAdj[passOneDistinguishers, 
                                    , drop = FALSE]
        tExprLinear_exprLinear_tProject <- gecd_MatrixChainMultiplication(tExprLinear_exprLinear, 
                                                                          t(project), verbose = verbose - 3)
        projectionNormInv <- gecd_MatrixChainMultiplication(project, 
                                                            tExprLinear_exprLinear_tProject, verbose = verbose - 
                                                              3)
        if (kappa(projectionNormInv) > 1e+06) {
          mesg <- sprintf("The 'projectionNormInv' matrix is of low rank (kappa = %g).  Ask for fewer numCellClasses.\n", 
                          kappa(projectionNormInv))
          stop(mesg)
        }
        projectionNorm <- solve(projectionNormInv)
        exprLinearBarAdj <- exprLinearBarAdj - gecd_MatrixChainMultiplication(exprLinearBarAdj, 
                                                                              tExprLinear_exprLinear_tProject, projectionNorm, 
                                                                              project, verbose = verbose - 3)
      }
      return(list(exprLinearBarAdj = exprLinearBarAdj))
    }
    stopifnot(nrow(exprLinearBarAdj) == nrow(exprLinear))
    stopifnot(ncol(exprLinearBarAdj) == ncol(exprLinear))
    stopifnot(nrow(tExprLinear_exprLinear) == ncol(exprLinear))
    stopifnot(ncol(tExprLinear_exprLinear) == ncol(exprLinear))
    len <- length(passOneDistinguishers)
    if (len == 1) {
      project <- exprLinearBarAdj[passOneDistinguishers[1], 
                                  , drop = FALSE]
      lengths <- gecd_MatrixChainMultiplication(exprLinearBarAdj, 
                                                tExprLinear_exprLinear, t(project), verbose = verbose - 
                                                  3)
      lengths <- lengths/sqrt(lengths[passOneDistinguishers[1]])
      bestDistinguishers <- order(lengths, decreasing = TRUE)[1:maxDistinguisherAlternatives]
      bestLengths <- lengths[bestDistinguishers]
      bestLengthsNormalized <- bestLengths/bestLengths[1]
      if (verbose > 2) {
        print(proc.time() - ptm)
        ptm <<- proc.time()
        cat(sprintf("Second pass: %d distinguisher(s) found for the cell class with pass-one distinguisher %s (length = %g)\n", 
                    maxDistinguisherAlternatives, rownames(exprLinear)[passOneDistinguishers[1]], 
                    lengths[passOneDistinguishers[1]]))
        distinguishers <- rownames(exprLinear)[bestDistinguishers]
        if (verbose > 3) {
          print(rbind(distinguishers, bestLengths, bestLengthsNormalized))
        }else {
          print(rbind(distinguishers, bestLengths))
        }
      }
      return(list(bestDistinguishers = bestDistinguishers, 
                  bestLengths = bestLengths, bestLengthsNormalized = bestLengthsNormalized))
    }else     {
      mid1 <- (len + 1)%/%2
      mid2 <- mid1 + 1
      projectOutExprLinearBarAdj1 <- projectOut(passOneDistinguishers[mid2:len], 
                                                exprLinearBarAdj, exprLinear, tExprLinear_exprLinear, 
                                                noDistinguishersYet)$exprLinearBarAdj
      ans1 <- Recall(passOneDistinguishers[1:mid1], projectOutExprLinearBarAdj1, 
                     exprLinear, tExprLinear_exprLinear, maxDistinguisherAlternatives, 
                     noDistinguishersYet = FALSE)
      projectOutExprLinearBarAdj2 <- projectOut(passOneDistinguishers[1:mid1], 
                                                exprLinearBarAdj, exprLinear, tExprLinear_exprLinear, 
                                                noDistinguishersYet)$exprLinearBarAdj
      ans2 <- Recall(passOneDistinguishers[mid2:len], projectOutExprLinearBarAdj2, 
                     exprLinear, tExprLinear_exprLinear, maxDistinguisherAlternatives, 
                     noDistinguishersYet = FALSE)
      return(list(bestDistinguishers = cbind(ans1$bestDistinguishers, 
                                             ans2$bestDistinguishers), bestLengths = cbind(ans1$bestLengths, 
                                                                                           ans2$bestLengths), bestLengthsNormalized = cbind(ans1$bestLengthsNormalized, 
                                                                                                                                            ans2$bestLengthsNormalized)))
    }
  }
  goodEnoughDistinguishers <- function(distinguishers, lengths, 
                                       threshold = 0.5, minLength = 1) {
    result <- distinguishers
    result[(lengths < threshold) & (row(lengths) > minLength)] <- NA
    tab <- table(result)
    MyGoodNames <- names(tab)[tab == 1]
    result[!(result %in% MyGoodNames) & !(row(result) == 
                                            1)] <- NA
    return(result)
  }
  ans <- expressionPreprocessing(exprLinear, genesymb, expressionQuantileForScale, 
                                 expressionQuantileForFilter, expressionConcentrationRatio, 
                                 probesWithGenesOnly, verbose)
  exprLinear <- ans$exprLinear
  tExprLinear_exprLinear <- ans$tExprLinear_exprLinear
  exprLinearBar <- ans$exprLinearBar
  exprLinearBarAdj <- ans$exprLinearBarAdj
  genesymb <- ans$genesymb
  rm(ans)
  if (verbose > 1) {
    print(proc.time() - ptm)
    ptm <<- proc.time()
    cat("gecd_CellDistinguisher: preprocessing complete\n")
  }
  ans <- findCandidateDistinguishers(exprLinear, tExprLinear_exprLinear, 
                                     exprLinearBarAdj, numCellClasses)
  passOneDistinguishers <- ans$passOneDistinguishers
  allLengths <- ans$allLengths
  rm(ans)
  if (verbose > 1) {
    print(proc.time() - ptm)
    ptm <<- proc.time()
    cat("gecd_CellDistinguisher: findCandidateDistinguishers complete\n")
  }
  ans <- adjustDistinguishers(passOneDistinguishers, exprLinearBarAdj, 
                              exprLinear, tExprLinear_exprLinear, maxDistinguisherAlternatives, 
                              noDistinguishersYet = TRUE)
  bestDistinguishers <- ans$bestDistinguishers
  bestLengths <- ans$bestLengths
  bestLengthsNormalized <- ans$bestLengthsNormalized
  rm(ans)
  if (verbose > 1) {
    print(proc.time() - ptm)
    ptm <<- proc.time()
    cat("gecd_CellDistinguisher: distinguisher adjustment complete\n")
  }
  bestDistinguishers <- goodEnoughDistinguishers(bestDistinguishers, 
                                                 bestLengthsNormalized, threshold = minAlternativesLengthsNormalized, 
                                                 minLength = minDistinguisherAlternatives)
  bestLengths[is.na(bestDistinguishers)] <- NA
  bestLengthsNormalized[is.na(bestDistinguishers)] <- NA
  if (nrow(bestDistinguishers) > 1) {
    bubbleNA <- function(x) {
      y <- rep(NA, length(x))
      z <- x[!is.na(x)]
      if (length(z) > 0) {
        y[1:length(z)] <- z
      }
      return(y)
    }
    bestDistinguishers <- apply(bestDistinguishers, 2, bubbleNA)
    bestLengths <- apply(bestLengths, 2, bubbleNA)
    bestLengthsNormalized <- apply(bestLengthsNormalized, 
                                   2, bubbleNA)
  }
  lastGoodRow <- sum(apply(bestDistinguishers, 1, function(rank) {
    sum(!is.na(rank)) > 0
  }))
  bestDistinguishers <- bestDistinguishers[1:lastGoodRow, , 
                                           drop = FALSE]
  bestLengths <- bestLengths[1:lastGoodRow, , drop = FALSE]
  bestLengthsNormalized <- bestLengthsNormalized[1:lastGoodRow, 
                                                 , drop = FALSE]
  bestDistinguishersGeneNames <- NULL
  passOneDistinguishersGeneNames <- NULL
  if (!is.null(genesymb)) {
    passOneDistinguishersGeneNames <- as.character(genesymb[passOneDistinguishers])
    bestDistinguishersGeneNames <- matrix(genesymb[bestDistinguishers], 
                                          nrow(bestDistinguishers), ncol(bestDistinguishers))
  }
  if (!is.null(rownames(exprLinear))) {
    passOneDistinguishers <- rownames(exprLinear)[passOneDistinguishers]
    bestDistinguishers <- matrix(rownames(exprLinear)[bestDistinguishers], 
                                 nrow(bestDistinguishers), ncol(bestDistinguishers))
  }
  bestDistinguishersFrame <- data.frame(biolproc = as.vector(col(bestDistinguishers)), 
                                        rank = as.vector(row(bestDistinguishers)), probeset = as.vector(bestDistinguishers), 
                                        genesymb = NA, length = as.vector(bestLengths), lengthNormalized = as.vector(bestLengthsNormalized), 
                                        stringsAsFactors = FALSE)
  passOneFrame <- data.frame(probeset = passOneDistinguishers, 
                             genesymb = NA, length = allLengths, stringsAsFactors = FALSE)
  if (!is.null(genesymb)) {
    passOneFrame$genesymb <- passOneDistinguishersGeneNames
    bestDistinguishersFrame$genesymb <- as.vector(bestDistinguishersGeneNames)
  }
  bestDistinguishersFrame <- bestDistinguishersFrame[!is.na(bestDistinguishersFrame$length), 
  ]
  rownames(passOneFrame) <- 1:nrow(passOneFrame)
  rownames(bestDistinguishersFrame) <- paste(bestDistinguishersFrame$biolproc, 
                                             bestDistinguishersFrame$rank, sep = "_")
  if (verbose > 1) {
    print(proc.time() - ptm)
    ptm <<- proc.time()
    cat("gecd_CellDistinguisher: post processing complete\n")
  }
  return(list(bestDistinguishers = bestDistinguishers, bestDistinguishersGeneNames = bestDistinguishersGeneNames, 
              bestLengths = bestLengths, bestLengthsNormalized = bestLengthsNormalized, 
              passOneDistinguishers = passOneDistinguishers, passOneDistinguishersGeneNames = passOneDistinguishersGeneNames, 
              passOneLengths = allLengths, bestDistinguishersFrame = bestDistinguishersFrame, 
              passOneFrame = passOneFrame))
}

my_gecd_DeconvolutionByDistinguishers <- function (exprLinear, bestDistinguishers, nonNegativeOnly = TRUE, 
          convexSolution = TRUE, verbose = 0) 
{
  exprLinear <- t(t(exprLinear)/colSums(exprLinear))
  if (is.matrix(bestDistinguishers)) {
    topDistinguishers <- bestDistinguishers[1, ]
  }else {
    topDistinguishers <- bestDistinguishers
  }
  if ((sum(is.na(topDistinguishers)) > 0) | (length(topDistinguishers) != 
                                             length(unique(topDistinguishers)))) {
    mesg <- "The supplied bestDistinguishers parameter contains NA or duplicates in the first row.\n"
    stop(mesg)
  }
  rowNormalization <- (exprLinear %*% colSums(exprLinear))[, 
                                                           1]
  exprLinearBar <- exprLinear/rowNormalization
  exprLinearBar[which(rowNormalization == 0), ] <- 0
  tExprLinear_exprLinear <- gecd_MatrixChainMultiplication(t(exprLinear), 
                                                           exprLinear, verbose = verbose - 3)
  project <- exprLinearBar[topDistinguishers, , drop = FALSE]
  projectionNormInv <- gecd_MatrixChainMultiplication(project, 
                                                      tExprLinear_exprLinear, t(project), verbose = verbose - 
                                                        3)
  if (kappa(projectionNormInv) > 1e+06) {
    mesg <- sprintf("The 'projectionNormInv' matrix is of low rank (kappa = %g).  Ask for fewer numCellClasses.\n", 
                    kappa(projectionNormInv))
    stop(mesg)
  }
  projectionNorm <- solve(projectionNormInv)
  convexCoefficientsUnconstrained <- gecd_MatrixChainMultiplication(exprLinearBar, 
                                                                    tExprLinear_exprLinear, t(project), projectionNorm, verbose = verbose - 
                                                                      3)
  if (convexSolution) {
    constraintCoefficients <- matrix(1, 1, length(topDistinguishers))
    constraintConstants <- matrix(1, nrow(exprLinear), 1)
    constraintNormInv <- gecd_MatrixChainMultiplication(constraintCoefficients, 
                                                        projectionNorm, t(constraintCoefficients), verbose = verbose - 
                                                          3)
    if (kappa(constraintNormInv) > 1e+06) {
      mesg <- sprintf("The 'constraintNormInv' matrix is of low rank (kappa = %g).  Ask for fewer numCellClasses.\n", 
                      kappa(constraintNormInv))
      stop(mesg)
    }
    constraintNorm <- solve(constraintNormInv)
    productNorm <- gecd_MatrixChainMultiplication(constraintNorm, 
                                                  constraintCoefficients, projectionNorm, verbose = verbose - 
                                                    3)
    convexCoefficientsConstrained <- convexCoefficientsUnconstrained + 
      (gecd_MatrixChainMultiplication(constraintConstants, 
                                      productNorm, verbose = verbose - 3) - gecd_MatrixChainMultiplication(convexCoefficientsUnconstrained, 
                                                                                                           t(constraintCoefficients), productNorm, verbose = verbose - 
                                                                                                             3))
  }else {
    convexCoefficientsConstrained <- convexCoefficientsUnconstrained
  }
  if (nonNegativeOnly) {
    epsilon <- 1e-06
    for (row in which(rowSums(convexCoefficientsConstrained < 
                              -epsilon) > 0)) {
      if (convexSolution) {
        constraintCoefficients <- rep(1, ncol(convexCoefficientsConstrained))
        constraintConstants <- 1
      }else {
        constraintCoefficients <- NULL
        constraintConstants <- NULL
      }
      replacement <- convexCoefficientsConstrained[row, 
      ]
      while (sum(replacement < -epsilon) > 0) {
        for (col in which(replacement < -epsilon)) {
          basis <- rep(0, ncol(convexCoefficientsConstrained))
          basis[col] <- 1
          constraintCoefficients <- rbind(constraintCoefficients, 
                                          basis)
          constraintConstants <- cbind(constraintConstants, 
                                       0)
        }
        constraintNorm <- solve(gecd_MatrixChainMultiplication(constraintCoefficients, 
                                                               projectionNorm, t(constraintCoefficients), 
                                                               verbose = verbose - 3))
        productNorm <- gecd_MatrixChainMultiplication(constraintNorm, 
                                                      constraintCoefficients, projectionNorm, verbose = verbose - 
                                                        3)
        replacement <- convexCoefficientsUnconstrained[row, 
                                                       , drop = FALSE] + (gecd_MatrixChainMultiplication(constraintConstants, 
                                                                                                         productNorm, verbose = verbose - 3) - gecd_MatrixChainMultiplication(convexCoefficientsUnconstrained[row, 
                                                                                                                                                                                                              , drop = FALSE], t(constraintCoefficients), 
                                                                                                                                                                              productNorm, verbose = verbose - 3))
      }
      convexCoefficientsConstrained[row, ] <- replacement
    }
    convexCoefficientsConstrained[convexCoefficientsUnconstrained < 
                                    epsilon] <- 0
  }
  cellSubclassSignatures <- convexCoefficientsConstrained * 
    (exprLinear %*% colSums(exprLinear))[, 1]
  cellSubclassSignatures <- t(t(cellSubclassSignatures)/colSums(cellSubclassSignatures))
  goodProbes <- unique(bestDistinguishers[!is.na(bestDistinguishers)])
  exprLinearDistinguishers <- exprLinear[goodProbes, ]
  exprLinearDistinguishers <- t(t(exprLinearDistinguishers)/colSums(exprLinearDistinguishers))
  cellSubclassSignaturesDistinguishers <- cellSubclassSignatures[goodProbes, 
  ]
  cellSubclassSignaturesDistinguishers <- t(t(cellSubclassSignaturesDistinguishers)/colSums(cellSubclassSignaturesDistinguishers))
  projectionNormInv <- t(cellSubclassSignaturesDistinguishers) %*% 
    cellSubclassSignaturesDistinguishers
  if (kappa(projectionNormInv) > 1e+06) {
    mesg <- sprintf("The 'projectionNormInv' matrix is of low rank (kappa = %g).  Ask for fewer numCellClasses.\n", 
                    kappa(projectionNormInv))
    stop(mesg)
  }
  projectionNorm <- solve(projectionNormInv)
  sampleCompositionsUnconstrained <- gecd_MatrixChainMultiplication(projectionNorm, 
                                                                    t(cellSubclassSignaturesDistinguishers), exprLinearDistinguishers, 
                                                                    verbose = verbose - 3)
  if (convexSolution) {
    constraintCoefficients <- matrix(1, length(topDistinguishers), 
                                     1)
    constraintConstants <- matrix(1, 1, ncol(exprLinear))
    constraintNormInv <- gecd_MatrixChainMultiplication(t(constraintCoefficients), 
                                                        projectionNorm, constraintCoefficients, verbose = verbose - 
                                                          3)
    if (kappa(constraintNormInv) > 1e+06) {
      mesg <- sprintf("The 'constraintNormInv' matrix is of low rank (kappa = %g).  Ask for fewer numCellClasses.\n", 
                      kappa(constraintNormInv))
      stop(mesg)
    }
    constraintNorm <- solve(constraintNormInv)
    productNorm <- gecd_MatrixChainMultiplication(projectionNorm, 
                                                  constraintCoefficients, constraintNorm, verbose = verbose - 
                                                    3)
    sampleCompositionsConstrained <- sampleCompositionsUnconstrained + 
      (gecd_MatrixChainMultiplication(productNorm, constraintConstants, 
                                      verbose = verbose - 3) - gecd_MatrixChainMultiplication(productNorm, 
                                                                                              t(constraintCoefficients), sampleCompositionsUnconstrained, 
                                                                                              verbose = verbose - 3))
  }else {
    sampleCompositionsConstrained <- sampleCompositionsUnconstrained
  }
  if (nonNegativeOnly) {
    epsilon <- 1e-06
    for (col in which(colSums(sampleCompositionsConstrained < 
                              -epsilon) > 0)) {
      if (convexSolution) {
        constraintCoefficients <- rep(1, nrow(sampleCompositionsConstrained))
        constraintConstants <- 1
      }else {
        constraintCoefficients <- NULL
        constraintConstants <- NULL
      }
      replacement <- sampleCompositionsConstrained[, col]
      while (sum(replacement < -epsilon) > 0) {
        for (row in which(replacement < -epsilon)) {
          basis <- rep(0, nrow(sampleCompositionsConstrained))
          basis[row] <- 1
          constraintCoefficients <- cbind(constraintCoefficients, 
                                          basis)
          constraintConstants <- rbind(constraintConstants, 
                                       0)
        }
        constraintNorm <- solve(gecd_MatrixChainMultiplication(t(constraintCoefficients), 
                                                               projectionNorm, constraintCoefficients, verbose = verbose - 
                                                                 3))
        productNorm <- gecd_MatrixChainMultiplication(projectionNorm, 
                                                      constraintCoefficients, constraintNorm, verbose = verbose - 
                                                        3)
        replacement <- sampleCompositionsUnconstrained[, 
                                                       col, drop = FALSE] + (gecd_MatrixChainMultiplication(productNorm, 
                                                                                                            constraintConstants, verbose = verbose - 3) - 
                                                                               gecd_MatrixChainMultiplication(productNorm, 
                                                                                                              t(constraintCoefficients), sampleCompositionsUnconstrained[, 
                                                                                                                                                                         col, drop = FALSE], verbose = verbose - 
                                                                                                                3))
      }
      sampleCompositionsConstrained[, col] <- replacement
    }
    sampleCompositionsConstrained[sampleCompositionsConstrained < 
                                    epsilon] <- 0
  }
  return(list(cellSubclassSignatures = cellSubclassSignatures, 
              sampleCompositions = sampleCompositionsConstrained))
}
