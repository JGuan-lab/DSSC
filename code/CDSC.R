
CDSC_3 <- function(X,c0, k, lambda1, lambda2,lambdaC,error= 10^-5,seedd,Ss,Sg,all_number = 500,ReturnL = F){
  X <- as.matrix(X)
  if(storage.mode(X) != "double"){
    storage.mode(X) <- "double"
  }
  
  gn <- dim(X)
  g <- gn[1]
  n <- gn[2]
  # Ss <- SM(t(X))
  # Sg <- SM(X)
  set.seed(seedd)#?????????????ӣ???֤????????һ??
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
    c <- c*((X%*%t(p)+ lambdaC*c0 +2*lambda2*Sg%*%c))/(c%*%p%*%t(p)+ lambdaC*c+ 2*lambda2*c%*%t(c)%*%c)
    p <- p*(t(c)%*%X+2*lambda1*p%*%Ss)/(t(c)%*%c%*%p+2*lambda1*p%*%t(p)%*%p)
    jump = jump+1
    
    rata_change_c <- Change_rata(c_init,c)
    rata_change_p <- Change_rata(p_init,p)
    if(rata_change_p == Inf || is.na(rata_change_p))rata_change_p = 100
    if(rata_change_c == Inf || is.na(rata_change_c))rata_change_c = 100
    
    if(ReturnL){
          L[jump] <- norm(X-c%*%p, type = "F")^2 + lambdaC*norm(c-c0, type = "F")^2 
                +lambda1*norm(Ss-t(p)%*%p, type = "F")^2 + +lambda2*norm(Sg-c%*%t(c), type = "F")^2
    }

    # setTxtProgressBar(pb, jump/all_number)
    # if(rata_change_c < error){
    
    if(all(rata_change_p < error & rata_change_c < error)){
      # cat("\n????????: ",jump)
      break
    }
    
  }
  # end_time <- Sys.time()
  # close(pb)
  # p <- SumEqual_1(p)
  c[which(c < 0.000001) ] <- 0
  if(ReturnL){
    return(list(c=c,p=p,jump=jump,L=L))
  }else{
    return(list(c=c,p=p,jump=jump))
  }
}

CDSC_2 <- function(X, k, lambda1, lambda2,lambdaC=0,error,seedd,Ss,Sg,all_number = 500,ReturnL = F){
  X <- as.matrix(X)
  if(storage.mode(X) != "double"){
    storage.mode(X) <- "double"
  }
  
  gn <- dim(X)
  g <- gn[1]
  n <- gn[2]
  # Ss <- SM(t(X))
  # Sg <- SM(X)
  set.seed(seedd)#?????????????ӣ???֤????????һ??
  c <- matrix(runif(g*k) ,g,k)
  set.seed(seedd)
  p <- matrix(runif(k*n) ,k,n)
  p <- NORM(p)
  # error = 10^-5
  jump = 0
  L = 0
  
  # pb <- txtProgressBar(style = 3)
  # star_time <- Sys.time()
  while(jump < all_number){
    c_init <- c
    p_init <- p
    c <- c*((X%*%t(p) +2*lambda2*Sg%*%c))/(c%*%p%*%t(p)+ 2*lambda2*c%*%t(c)%*%c)
    p <- p*(t(c)%*%X+2*lambda1*p%*%Ss)/(t(c)%*%c%*%p+2*lambda1*p%*%t(p)%*%p)
    jump = jump+1
    
    rata_change_c <- Change_rata(c_init,c)
    rata_change_p <- Change_rata(p_init,p)
    if(ReturnL){
      L[jump] <- norm(X-c%*%p, type = "F")^2 + lambda1*norm(Ss-t(p)%*%p, type = "F")^2 + +lambda2*norm(Sg-c%*%t(c), type = "F")^2
    }
    # setTxtProgressBar(pb, jump/all_number)
    # if(rata_change_c < error){
    if(all(c(rata_change_p,rata_change_c) < error)){
      # cat("\n????????: ",jump)
      break
    }
    
  }
  # end_time <- Sys.time()
  # close(pb)
  # p <- SumEqual_1(p)
  c[which(c < 0.000001) ] <- 0
  if(ReturnL){
    return(list(c=c,p=p,jump=jump,L=L))
  }else{
    return(list(c=c,p=p,jump=jump))
  }
  
}
CDSC_3_mix <- function(X,c0, k, lambda1, lambda2,lambdaC,error,seedd, all_number = 500,
                       Ss=matrix(1,dim(X)[2],dim(X)[2]),Sg=matrix(1,dim(X)[1],dim(X)[1])){
  X <- as.matrix(X)
  if(storage.mode(X) != "double"){
    storage.mode(X) <- "double"
  }
  
  gn <- dim(X)
  g <- gn[1]
  n <- gn[2]
  # Ss <- SM(t(X))
  # Sg <- SM(X)
  set.seed(seedd)
  c <- matrix(runif(g*k) ,g,k)
  set.seed(seedd)
  p <- matrix(runif(k*n) ,k,n)
  # error = 10^-5
  jump = 0
  L = 0
  # all_number = 500
  # pb <- txtProgressBar(style = 3)
  # star_time <- Sys.time()
  while(jump < all_number){
    c_init <- c
    p_init <- p
    c <- c*((X%*%t(p)+ lambdaC*c0 +2*lambda2*Sg%*%c))/(c%*%p%*%t(p)+ lambdaC*c+ 2*lambda2*c%*%t(c)%*%c)
    p <- p*(t(c)%*%X+2*lambda1*p%*%Ss)/(t(c)%*%c%*%p+2*lambda1*p%*%t(p)%*%p)
    jump = jump+1
    
  }
  # c[which(c < 0.000001) ] <- 0
  return(list(c=c,p=p,jump=jump))
}
########


#######
SM <- function(A){
  ##????cmf-impute?????ģ?d = 1-pearson,s = 1/(1+d)
  score = 1/(2-cor(t(A)))
  
  # score=matrix(1,dim(A)[1],dim(A)[1]) - cor(t(A))
  # for (i in 1:dim(score)[1]) {
  #   for (j in 1:dim(score)[2]) {
  #     score[i,j] <- 1/(1+score[i,j])
  #   }
  # }
  return(score)
}
########
NORM <- function(A){
  v <- colSums(A)
  for (i in 1:length(v)) {
    if(v[i]!=0){
      A[,i] <- A[,i]/v[i]
    }
  }
  # d <- diag(v)
  # A <- A%*%solve(d)
  return(A)
}

SumEqual_1 <- function(p){
  return(apply(p,2,function(x) x/sum(x)))
}
########
Change_rata <- function(x_init,x){
  a = norm(x_init - x,"F")^2
  b = norm(x_init,"F")^2
  return(a/b)
}


