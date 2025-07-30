#' Runs DSSC
#'
#' @param data_bulk the input dropout bulk data.
#'
#' @param data_ref the reference GEP matrix
#'
#' @param k Number of cell types used for matrix initialization
#'
#' @param parameter the vector of parameters.
#' lambda1 is the value of lambda1 in the mathematical model;
#' lambda2 is the value of lambda2 in the mathematical model;
#' lambdaC is the value of lambda3 in the mathematical model.
#'
#' @param error the threshold of the error between the current imputed matrix and the previous one.
#' Default is 1e-5.
#'
#' @param seedd Random seed, default is 44.
#'
#' @param Ss sample similarity matrix, obtained through the SM function, eg: SM(t(data_bulk))
#'
#' @param Sg gene similarity matrix, obtained through the SM function, eg: SM (data_bulk)
#'
#' @param all_number the maximum iterations, the default is 500.
#'
#' @param ReturnL Returns the value L of the objective function for each iteration.
#' Default is F.
#'
#' @examples
#' # Set up the parameter and matrix used in DSSC
#' Ss <- SM(t(data_bulk))
#' Sg <- SM(data_bulk)
#' lambda1 <- 1e-03
#' lambda2 <- 1e+00
#' lambdaC <- 1e+01
#'
#' # Run DSSC
#' result <- DSSC(data_bulk,data_ref,k,lambda1, lambda2, lambdaC, Ss, Sg)
#'
#' @return Two matrices, the first is the GEP matrix C, and the second is the cell type proportion matrix P
#'
#' @rdname DSSC
#'
#' @export

DSSC <- function(data_bulk, data_ref, lambda1, lambda2, lambdaC, error= 10^-5,
                 seedd=44, Ss, Sg, all_number = 500, ReturnL = F){
    data_bulk <- as.matrix(data_bulk)
    if(storage.mode(data_bulk) != "double"){
        storage.mode(data_bulk) <- "double"
    }
    if (is.null(dim(data_ref))) {
        k = data_ref
        data_ref = 0
    } else {
        k = dim(data_ref)[2]
    }
    gn <- dim(data_bulk)
    g <- gn[1]
    n <- gn[2]
    # Ss <- SM(t(data_bulk))
    # Sg <- SM(data_bulk)
    set.seed(seedd)
    c <- matrix(runif(g*k) ,g,k)
    set.seed(seedd)
    p <- matrix(runif(k*n) ,k,n)
    p <- NORM(p)
    # error = 10^-5
    jump = 0
    L = array()
    while(jump < all_number){
        c_init <- c
        p_init <- p
        c <- c*((data_bulk%*%t(p)+ lambdaC*data_ref +2*lambda2*Sg%*%c))/(c%*%p%*%t(p)+ lambdaC*c+ 2*lambda2*c%*%t(c)%*%c)
        p <- p*(t(c)%*%data_bulk+2*lambda1*p%*%Ss)/(t(c)%*%c%*%p+2*lambda1*p%*%t(p)%*%p)
        jump = jump+1

        rata_change_c <- Change_rata(c_init,c)
        rata_change_p <- Change_rata(p_init,p)
        if(rata_change_p == Inf || is.na(rata_change_p))rata_change_p = 100
        if(rata_change_c == Inf || is.na(rata_change_c))rata_change_c = 100

        if(ReturnL){
            L[jump] <- norm(data_bulk-c%*%p, type = "F")^2 + lambdaC*norm(c-data_ref, type = "F")^2
            +lambda1*norm(Ss-t(p)%*%p, type = "F")^2 + +lambda2*norm(Sg-c%*%t(c), type = "F")^2
        }

        if(all(rata_change_p < error & rata_change_c < error)){
            break
        }

    }
    c[which(c < 0.000001) ] <- 0
    if(ReturnL){
        return(list(c=c,p=p,jump=jump,L=L))
    }else{
        return(list(c=c,p=p,jump=jump))
    }
}

# get  similarity matrix
SM <- function(A){

    score = 1/(2-cor(t(A)))

    return(score)
}


NORM <- function(A){
    v <- colSums(A)
    for (i in 1:length(v)) {
        if(v[i]!=0){
            A[,i] <- A[,i]/v[i]
        }
    }
    d <- diag(v)
    A <- A%*%solve(d)
    return(A)
}

# Make the sum of cell types corresponding to the proportional matrix sample one
SumEqual_1 <- function(p){
    return(apply(p,2,function(x) x/sum(x)))
}

# the threshold of the error between the current imputed matrix and the previous one.
Change_rata <- function(x_init,x){
    a = norm(x_init - x,"F")^2
    b = norm(x_init,"F")^2
    return(a/b)
}

simulation <- function(scData,leastNum=50, plotmarker = F,norm1 = "CPM",log2.threshold = 1) {
    scData$simulate1 <- scSimulateSplit(scData,leastNum=50, plotmarker = F,norm1 = "CPM",log2.threshold = 1)
    bulkData <- scSimulateShift(scData,"all",standardization=TRUE)
    return(bulkData)
}
cross_validation <- function(bulk,
                             ref,
                             n_folds = 5,
                             seedd = 44,
                             TerCondition = 10^-8,
                             lambda1 = c(0,10^-3,10^-2,10^-1,1,10),
                             lambda2 = c(0,10^-3,10^-2,10^-1,1,10),
                             lambdaC = c(0,10^-1,10^0,10^1,10^2,10^3),
                             max_num = 1500){

    # library(rBayesianOptimization)
    set.seed(seedd)
    vec <- vector()
    for (i in 1:n_folds) {
        vec <- c(vec, rep(i,nrow(bulk)*ncol(bulk)/n_folds))
    }
    sample_matrix <- matrix(sample(vec), nrow = nrow(bulk))
    mask <- list();mask_test <- list()
    for (i in 1:n_folds) {
        mask[[i]] <- sample_matrix
        mask[[i]][mask[[i]] == i] <- 0
        mask[[i]][mask[[i]] > 0] <- 1

        mask_test[[i]] <- sample_matrix
        mask_test[[i]][mask_test[[i]] != i] <- 0
        mask_test[[i]][mask_test[[i]]  > 0] <- 1
    }
    pb <- txtProgressBar(style = 3)
    star_time <- Sys.time()
    Ss <- list()
    Sg <- list()
    for (i in 1:n_folds) {
        Ss[[i]] <- SM(t(bulk*mask[[i]]))
        Sg[[i]] <- SM(bulk*mask[[i]])
    }
    num = 0
    para_lambda_44_8 = NULL

    library(dplyr)
    for (dir_i in 1:length(lambda1)){
        for (dir_j in 1:length(lambda2)){
            for (dir_k in 1:length(lambdaC)){
                result_CDSC <- list()
                result1 <- data.frame(0,0,0,0)
                for (i in 1:n_folds) {
                    result_CDSC[[i]] <- DSSC(bulk*mask[[i]],
                                             ref,
                                             lambda1[dir_i], lambda2[dir_j], lambdaC[dir_k],
                                             TerCondition,seedd,
                                             Ss[[i]],Sg[[i]],
                                             max_num)
                    result_c <- result_CDSC[[i]][[1]]
                    result_p <- result_CDSC[[i]][[2]]
                    temp <- try(Row_label(result_c, as.matrix(ref), leastnum = 3),silent = FALSE)
                    if ('try-error' %in% class(temp)) {
                        ctlabels <- colnames(ref)
                    } else {
                        ctlabels <- Row_label(result_c, as.matrix(ref), leastnum = 3)
                    }
                    colnames(result_c) <- ctlabels
                    rownames(result_p) <- ctlabels
                    result1 = result1 + data.frame(getPearsonRMSE(ref, result_c),
                                                   getPearsonRMSE(bulk*mask_test[[i]], (result_c%*%result_p)*mask_test[[i]]))
                    num = num + 1
                }
                result1 = result1 / n_folds
                para_lambda_44_8 =  rbind(para_lambda_44_8, data.frame(result1, lambda1[dir_i], lambda2[dir_j], lambdaC[dir_k]))

                setTxtProgressBar(pb, num/(length(lambda1)*length(lambda2)*length(lambdaC)*n_folds))
            }
        }
    }
    colnames(para_lambda_44_8) <- c("RMSE.C", "PCC.C", "RMSE.T", "PCC.T", "lambda1", "lambda2", "lambdaC")
    return(para_lambda_44_8)
}