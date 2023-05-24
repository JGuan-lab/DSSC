SCDC_prop_my <- function (bulk.eset, sc.eset, ct.varname, sample, ct.sub, iter.max = 1000, 
          nu = 1e-04, epsilon = 0.01, truep = NULL, weight.basis = T, 
          ct.cell.size = NULL, Transform_bisque = F, ...) 
{
  bulk.eset <- bulk.eset[rowSums(exprs(bulk.eset)) > 0, , drop = FALSE]
  ct.sub <- intersect(ct.sub, unique(sc.eset@phenoData@data[, 
                                                            ct.varname]))
  sc.basis <- SCDC_basis_my(x = sc.eset, ct.sub = ct.sub, ct.varname = ct.varname, 
                         sample = sample, ct.cell.size = ct.cell.size)
  
  # rownames(sc.basis$sigma) <- rownames(as.matrix(res))
  commongenes <- intersect(rownames(sc.basis$basis.mvw), rownames(bulk.eset))
  if (length(commongenes) < 0.2 * min(dim(sc.eset)[1], dim(bulk.eset)[1])) {
    stop("Too few common genes!")
  }
  message(paste("Used", length(commongenes), "common genes..."))
  if (weight.basis) {
    basis.mvw <- sc.basis$basis.mvw[commongenes, ct.sub]
  }else {
    basis.mvw <- sc.basis$basis[commongenes, ct.sub]
  }
  if (Transform_bisque) {
    GenerateSCReference <- function(sc.eset, ct.sub, ct.varname) {
      cell.labels <- base::factor(sc.eset[[ct.varname]])
      all.cell.types <- base::levels(cell.labels)
      aggr.fn <- function(ct.sub) {
        base::rowMeans(Biobase::exprs(sc.eset)[, cell.labels == 
                                                 ct.sub, drop = F])
      }
      template <- base::numeric(base::nrow(sc.eset))
      sc.ref <- base::vapply(all.cell.types, aggr.fn, template)
      return(sc.ref)
    }
    sc.ref <- GenerateSCReference(sc.eset, ct.sub, ct.varname)[commongenes, 
                                                               , drop = F]
    ncount <- table(sc.eset@phenoData@data[, sample], sc.eset@phenoData@data[, 
                                                                             ct.varname])
    true.prop <- ncount/rowSums(ncount, na.rm = T)
    sc.props <- round(true.prop[complete.cases(true.prop), 
                                ], 2)
    Y.train <- sc.ref %*% t(sc.props[, colnames(sc.ref)])
    X.pred <- exprs(bulk.eset)[commongenes, ]
    sample.names <- base::colnames(Biobase::exprs(bulk.eset))
    template <- base::numeric(base::length(sample.names))
    base::names(template) <- sample.names
    SemisupervisedTransformBulk <- function(gene, Y.train, 
                                            X.pred) {
      Y.train.scaled <- base::scale(Y.train[gene, , drop = T])
      Y.center <- base::attr(Y.train.scaled, "scaled:center")
      Y.scale <- base::attr(Y.train.scaled, "scaled:scale")
      n <- base::length(Y.train.scaled)
      shrink.scale <- base::sqrt(base::sum((Y.train[gene, 
                                                    , drop = T] - Y.center)^2)/n + 1)
      X.pred.scaled <- base::scale(X.pred[gene, , drop = T])
      Y.pred <- base::matrix((X.pred.scaled * shrink.scale) + 
                               Y.center, dimnames = base::list(base::colnames(X.pred), 
                                                               gene))
      return(Y.pred)
    }
    Y.pred <- base::matrix(base::vapply(X = commongenes, 
                                        FUN = SemisupervisedTransformBulk, FUN.VALUE = template, 
                                        Y.train, X.pred, USE.NAMES = TRUE), nrow = base::length(sample.names))
    indices <- base::apply(Y.pred, MARGIN = 2, FUN = function(column) {
      base::anyNA(column)
    })
    if (base::any(indices)) {
      if (sum(!indices) == 0) {
        base::stop("Zero genes left for decomposition.")
      }
      Y.pred <- Y.pred[, !indices, drop = F]
      sc.ref <- sc.ref[!indices, , drop = F]
    }
    results <- base::as.matrix(base::apply(Y.pred, 1, function(b) {
      sol <- lsei::pnnls(sc.ref, b, sum = 1)
      return(sol$x)
    }))
    prop.est.mvw <- t(results)
    colnames(prop.est.mvw) <- colnames(sc.ref)
    rownames(prop.est.mvw) <- colnames(bulk.eset)
    yhat <- sc.ref %*% results
    colnames(yhat) <- colnames(bulk.eset)
    yobs <- exprs(bulk.eset)
    yeval <- SCDC_yeval(y = yobs, yest = yhat, yest.names = c("SCDC"))
    peval <- NULL
    if (!is.null(truep)) {
      peval <- SCDC_peval(ptrue = truep, pest = prop.est.mvw, 
                          pest.names = c("SCDC"), select.ct = ct.sub)
    }
  }else {
    xbulk <- getCPM0(exprs(bulk.eset)[commongenes, , drop = F])
    sigma <- sc.basis$sigma[commongenes, ct.sub]
    ALS.S <- sc.basis$sum.mat[ct.sub]
    N.bulk <- ncol(bulk.eset)
    valid.ct <- (colSums(is.na(sigma)) == 0) & (colSums(is.na(basis.mvw)) == 
                                                  0) & (!is.na(ALS.S))
    if (sum(valid.ct) <= 1) {
      stop("Not enough valid cell type!")
    }
    message(paste("Used", sum(valid.ct), "cell types in deconvolution..."))
    basis.mvw <- basis.mvw[, valid.ct]
    ALS.S <- ALS.S[valid.ct]
    sigma <- sigma[, valid.ct]
    prop.est.mvw <- NULL
    yhat <- NULL
    yhatgene.temp <- rownames(basis.mvw)
    for (i in 1:N.bulk) {
      basis.mvw.temp <- basis.mvw
      xbulk.temp <- xbulk[, i] * 100
      sigma.temp <- sigma
      message(paste(colnames(xbulk)[i], "has common genes", 
                    sum(xbulk[, i] != 0), "..."))
      lm <- nnls::nnls(A = basis.mvw.temp, b = xbulk.temp)
      delta <- lm$residuals
      wt.gene <- 1/(nu + delta^2 + colSums((lm$x * ALS.S)^2 * 
                                             t(sigma.temp)))
      x.wt <- xbulk.temp * sqrt(wt.gene)
      b.wt <- sweep(basis.mvw.temp, 1, sqrt(wt.gene), "*")
      lm.wt <- nnls::nnls(A = b.wt, b = x.wt)
      prop.wt <- lm.wt$x/sum(lm.wt$x)
      delta <- lm.wt$residuals
      for (iter in 1:iter.max) {
        wt.gene <- 1/(nu + delta^2 + colSums((lm.wt$x * 
                                                ALS.S)^2 * t(sigma.temp)))
        x.wt <- xbulk.temp * sqrt(wt.gene)
        b.wt <- sweep(basis.mvw.temp, 1, sqrt(wt.gene), 
                      "*")
        lm.wt <- nnls::nnls(A = b.wt, b = x.wt)
        delta.new <- lm.wt$residuals
        prop.wt.new <- lm.wt$x/sum(lm.wt$x)
        if (sum(abs(prop.wt.new - prop.wt)) < epsilon) {
          prop.wt <- prop.wt.new
          delta <- delta.new
          R2 <- 1 - var(xbulk.temp - basis.mvw.temp %*% 
                          as.matrix(lm.wt$x))/var(xbulk.temp)
          message("WNNLS Converged at iteration ", iter)
          break
        }
        prop.wt <- prop.wt.new
        delta <- delta.new
      }
      R2 <- 1 - var(xbulk.temp - basis.mvw.temp %*% as.matrix(lm.wt$x))/var(xbulk.temp)
      prop.est.mvw <- rbind(prop.est.mvw, prop.wt)
      yhat.temp <- basis.mvw.temp %*% as.matrix(lm.wt$x)
      yhatgene.temp <- intersect(rownames(yhat.temp), yhatgene.temp)
      yhat <- cbind(yhat[yhatgene.temp, ], yhat.temp[yhatgene.temp, 
                                                     ])
    }
    colnames(prop.est.mvw) <- colnames(basis.mvw)
    rownames(prop.est.mvw) <- colnames(xbulk)
    colnames(yhat) <- colnames(xbulk)
    yobs <- exprs(bulk.eset)
    yeval <- SCDC_yeval(y = yobs, yest = yhat, yest.names = c("SCDC"))
    peval <- NULL
    if (!is.null(truep)) {
      peval <- SCDC_peval(ptrue = truep, pest = prop.est.mvw, 
                          pest.names = c("SCDC"), select.ct = ct.sub)
    }
  }
  return(list(prop.est.mvw = prop.est.mvw, basis.mvw = basis.mvw, 
              yhat = yhat, yeval = yeval, peval = peval))
}


SCDC_basis_my <- function (x, ct.sub = NULL, ct.varname, sample, ct.cell.size = NULL) 
{
  if (is.null(ct.sub)) {
    ct.sub <- unique(x@phenoData@data[, ct.varname])
  }
  ct.sub <- ct.sub[!is.na(ct.sub)]
  x.sub <- x[, x@phenoData@data[, ct.varname] %in% ct.sub]
  x.sub <- x.sub[rowSums(exprs(x.sub)) > 0, ]
  countmat <- exprs(x.sub)
  ct.id <- droplevels(as.factor(x.sub@phenoData@data[, ct.varname]))
  sample.id <- as.character(x.sub@phenoData@data[, sample])
  ct_sample.id <- paste(ct.id, sample.id, sep = "%")
  mean.mat <- sapply(unique(ct_sample.id), function(id) {
    y = as.matrix(countmat[, ct_sample.id %in% id])
    apply(y, 1, sum, na.rm = TRUE)/sum(y)
  })
  mean.id <- do.call("rbind", strsplit(unique(ct_sample.id), 
                                       split = "%"))
  sigma <- sapply(unique(mean.id[, 1]), function(id) {
    y = mean.mat[, mean.id[, 1] %in% id]
    if (is.null(dim(y))) {
      res = rep(0, length(y))
      message("Warning: the cell type [", id, "] is only available in at most 1 subject!")
    }else {
      res = apply(y, 1, var, na.rm = TRUE)
      sigm.names <- rownames(as.matrix(res))
    }
    return(res)
  })
  rownames(sigma) <- rownames(mean.mat)
  
  sum.mat2 <- sapply(unique(sample.id), function(sid) {
    sapply(unique(ct.id), function(id) {
      y = as.matrix(countmat[, ct.id %in% id & sample.id %in% 
                               sid])
      if (ncol(y) > 0) {
        out = sum(y)/ncol(y)
      }
      else {
        out = 0
      }
      return(out)
    })
  })
  rownames(sum.mat2) <- unique(ct.id)
  colnames(sum.mat2) <- unique(sample.id)
  if (is.null(ct.cell.size)) {
    sum.mat <- rowMeans(sum.mat2, na.rm = T)
  }else {
    if (is.null(names(ct.cell.size))) {
      message("Cell size factor vector requires cell type names...")
      break
    }else {
      sum.mat <- ct.cell.size
    }
  }
  basis <- sapply(unique(mean.id[, 1]), function(id) {
    z <- sum.mat[mean.id[, 1]]
    mean.mat.z <- t(t(mean.mat) * z)
    y = as.matrix(mean.mat.z[, mean.id[, 1] %in% id])
    apply(y, 1, mean, na.rm = TRUE)
  })
  my.max <- function(x, ...) {
    y <- apply(x, 1, max, na.rm = TRUE)
    if (median(y, na.rm = T) == 0) {
      outx = y
    }else {
      outx = y/median(y, na.rm = T)
    }
    return(outx)
  }
  var.adj <- sapply(unique(sample.id), function(sid) {
    my.max(sapply(unique(ct.id), function(id) {
      y = countmat[, ct.id %in% id & sample.id %in% sid, 
                   drop = FALSE]
      if (ncol(y) > 0) {
        out = apply(y, 1, var, na.rm = T)
      }else {
        out = rep(0, nrow(y))
      }
      return(out)
    }), na.rm = T)
  })
  colnames(var.adj) <- unique(sample.id)
  q15 <- apply(var.adj, 2, function(zz) {
    z1 = min(zz[zz > 0])
    z2 = quantile(zz, 0.15, na.rm = T)
    return(max(z1, z2))
  })
  q85 <- apply(var.adj, 2, quantile, probs = 0.85, na.rm = T)
  var.adj.q <- t(apply(var.adj, 1, function(y) {
    y[y < q15] <- q15[y < q15]
    y[y > q85] <- q85[y > q85]
    return(y)
  }))
  var.adj.q <- t(apply(var.adj, 1, function(y) {
    y[y < q15] <- q15[y < q15]
    y[y > q85] <- q85[y > q85]
    return(y)
  }))
  message("Creating Basis Matrix adjusted for maximal variance weight")
  mean.mat.mvw <- sapply(unique(ct_sample.id), function(id) {
    sid = unlist(strsplit(id, "%"))[2]
    y = as.matrix(countmat[, ct_sample.id %in% id])
    yy = sweep(y, 1, sqrt(var.adj.q[, sid]), "/")
    apply(yy, 1, sum, na.rm = TRUE)/sum(yy)
  })
  basis.mvw <- sapply(unique(mean.id[, 1]), function(id) {
    z <- sum.mat[mean.id[, 1]]
    mean.mat.z <- t(t(mean.mat.mvw) * z)
    y = as.matrix(mean.mat.z[, mean.id[, 1] %in% id])
    apply(y, 1, mean, na.rm = TRUE)
  })
  basis.mvw <- basis.mvw[, ct.sub]
  sigma <- sigma[, ct.sub]
  basis <- basis[, ct.sub]
  sum.mat <- sum.mat[ct.sub]
  return(list(basis = basis, sum.mat = sum.mat, sigma = sigma, 
              basis.mvw = basis.mvw, var.adj.q = var.adj.q))
}
