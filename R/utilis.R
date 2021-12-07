### Create autoregressive of order 1 matrix
AR.matrix <- function(n,rho){
  SIGMA<-matrix(0,ncol=n,nrow=n)
  for (i in 1:n){
    for (j in 1:n){
      SIGMA[i,j]<-rho^abs(i-j)
    }
  }
  return(SIGMA)
}

### Calculate value of Bayesian information Criterion
BIC<-function(U,Y.data,X.data){

  N <- nrow(Y.data)
  # sigmahat<-sum(diag((1/nrow(Y.data))*(t(Y.data-X.data%*%U)%*%(Y.data-X.data%*%U))))
  sigmahat<-sum(t(Y.data-X.data%*%U)%*%(Y.data-X.data%*%U))/N
  BIC<-N*log(sigmahat) + length(which(U!=0))*log(N)
  return(BIC)
}

### Calculate projection error
projection.error<-function(x,y){
  sqrt(sum((x%*%solve(t(x)%*%x)%*%t(x) - y%*%solve(t(y)%*%y)%*%t(y))^2))
}

