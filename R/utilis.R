### Calculate value of Bayesian information Criterion
BIC<-function(U,Y.data,X.data){

  N <- nrow(Y.data)
  # sigmahat<-sum(diag((1/nrow(Y.data))*(t(Y.data-X.data%*%U)%*%(Y.data-X.data%*%U))))
  sigmahat<-sum(t(Y.data-X.data%*%U)%*%(Y.data-X.data%*%U))/N
  BIC<-N*log(sigmahat) + length(which(U!=0))*log(N)
  return(BIC)
}

##' @export
### Calculate projection error
projection.error<-function(x,y){
  sqrt(sum((x%*%solve(t(x)%*%x)%*%t(x) - y%*%solve(t(y)%*%y)%*%t(y))^2))
}

##' @export
### Calculate tail probability
Tail.Prob <- function(X,Y){
  N <- nrow(X)
  dfX <- ncol(X)
  dfY <- ncol(Y)
  det.BIG <- det(cov(cbind(X,Y)))
  det.X <- det(cov(X))
  if(ncol(Y) == 1) {det.Y <- abs(var(Y))} else{det.Y <- det(cov(Y))}
  Wilks.LAMBDA <- try( as.vector(det.BIG / ( det.X * det.Y )) )
  scaled.Wilks.LAMBDA <- try( -(N-1-0.5*(dfX+dfY+1))*log(Wilks.LAMBDA) )

  Test.Stat <- try(1 - pchisq(scaled.Wilks.LAMBDA, df=dfX * dfY))
  return(Test.Stat)
}
