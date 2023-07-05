#' Runs CDSC
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
#' # Set up the parameter and matrix used in CDSC
#' Ss <- SM(t(data_bulk))
#' Sg <- SM(data_bulk)
#' lambda1 <- 1e-03
#' lambda2 <- 1e+00
#' lambdaC <- 1e+01
#' 
#' # Run CDSC
#' result <- CDSC(data_bulk,data_ref,k,lambda1, lambda2, lambdaC, Ss, Sg)
#'
#' @return Two matrices, the first is the GEP matrix C, and the second is the cell type proportion matrix P
#'
#' @rdname CDSC
#'
#' @export

CDSC <- function(data_bulk,data_ref, k, lambda1, lambda2,lambdaC,error= 10^-5,seedd=44,Ss,Sg,all_number = 500,ReturnL = F){
  data_bulk <- as.matrix(data_bulk)
  if(storage.mode(data_bulk) != "double"){
    storage.mode(data_bulk) <- "double"
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


# NORM <- function(A){
#   v <- colSums(A)
#   for (i in 1:length(v)) {
#     if(v[i]!=0){
#       A[,i] <- A[,i]/v[i]
#     }
#   }
#   # d <- diag(v)
#   # A <- A%*%solve(d)
#   return(A)
# }

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
