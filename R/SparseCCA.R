
#' @title
#' Implement SparseCCA
#'
#' @description
#' Implement an iterative penalized least squares approach to sparse canonical correlation analysis (SparseCCA)
#' with various penalty functions. Modified Wilms, Ines, and Christophe Croux. "Robust sparse canonical correlation analysis." BMC systems biology 10.1 (2016): 1-13.
#' The original code is accessible https://sites.google.com/view/iwilms/software?authuser=0
#'
#'
#' ### INPUT
#' @param X                   : \eqn{n} by \eqn{p} exposure data matrix, where \eqn{n} is sample size and \eqn{p} is number of exposures.
#' @param Y                   : \eqn{n} by \eqn{q} outcome data matrix, where \eqn{n} is sample size and \eqn{q} is number of outcomes.
#' @param Xmethod             : penalty function for the exposure, i.e. penalty function when regressing Y onto X. Possible values are:
#'   - "lasso": Lasso
#'   - "alasso": Adaptive Lasso
#'   - "gglasso": Group Lasso
#'   - "SGL": Sparse Group Lasso
#' @param Ymethod             : penalty function for the outcome, i.e. penalty function when regressing X onto Y. Possible values are:
#'   - "lasso": Lasso
#'   - "alasso": Adaptive Lasso
#'   - "gglasso": Group Lasso
#'   - "SGL": Sparse Group Lasso
#'   - "OLS": Ordinary Least Square
#' @param init.method         : "Possible two values are:
#'   - "lasso": When high-dimensional exposure data (i.e. n<p)
#'   - "OLS": Ordinary Least Square
#'   - "SVD": Singular Value Decomposition
#' @param X.groupidx          : A vector of length \eqn{p} that indicates grouping structure of exposure `X`
#' @param Y.groupidx          : A vector of length \eqn{q} that indicates grouping structure of outcome `Y`
#' @param standardize         : A standardization (center and scale) of exposure `X` and outcome `Y`
#' @param max.iter            : A maximum number of iterations
#' @param conv                : A tolerance value for convergence \eqn{epsilon}
#'
#' @return
#' The function returns a list of 6 objects according to the following order:
#'   - ALPHA          : An estimated canonical vector of length \eqn{p} corresponding to the exposure data \eqn{X}
#'   - BETA           : An estimated canonical vector of length \eqn{q} corresponding to the outcome data \eqn{Y}
#'   - U_ALL          : An estimated canonical variates \eqn{X\alpha} corresponding to the exposure data \eqn{X} which is a vector of length \eqn{n}
#'   - V_ALL          : An estimated canonical variates \eqn{Y\beta} corresponding to the outcome data \eqn{Y} which is a vector of length \eqn{n}
#'   - lamdbaA        : A numeric value of the sparsity parameter lambdaA
#'   - lamdbaB        : A numeric value of the sparsity parameter lambdaB
#'   - it             : Number of iterations until convergence
#'
#' @export
#'
#' @examples
#' data.list <- generate.data(n=500)
#' DATA.X <- data.list$DATA.X
#' DATA.Y <- data.list$DATA.Y
#' SparseCCA.result <- SparseCCA(X=DATA.X,Y=DATA.Y,Xmethod="SGL",Ymethod="lasso",X.groupidx=c(rep(1,4),rep(2,3),rep(3,2),rep(1,3)))
#' str(SparseCCA.result)
#' projection.error(data.list$TRUE.ALPHA,SparseCCA.result$ALPHA)
#' projection.error(data.list$TRUE.BETA,SparseCCA.result$BETA)
#'
#'

SparseCCA<-function(X,Y,Xmethod=c("lasso","alasso","gglasso","SGL"),Ymethod=c("lasso","alasso","SGL","SLR"),init.method="SVD",X.groupidx=NULL,Y.groupidx=NULL,standardize=T,max.iter=50,conv=5*10^-2){
  ### Function to perform Sparse Canonical Correlation Analysis using alternating regressions

  ### AUXILIARY FUNCTIONS
  # NORMALIZATION_UNIT
  # BIC
  # Sparse.alternating

  ### STORE RESULTS
  # ALPHA_ALL<-matrix(NA,ncol=rank,nrow=ncol(X))
  # BETA_ALL<-matrix(NA,ncol=rank,nrow=ncol(Y))
  # U_ALL<-matrix(NA,ncol=rank,nrow=nrow(X))
  # V_ALL<-matrix(NA,ncol=rank,nrow=nrow(Y))
  # cancors<-matrix(NA,ncol=rank,nrow=1)
  # lambdaA_ALL<-matrix(NA,ncol=rank,nrow=max.iter)
  # lambdaB_ALL<-matrix(NA,ncol=rank,nrow=max.iter)
  # iterations<-matrix(NA,ncol=rank,nrow=1)
  # obj.it<-matrix(NA,ncol=rank,nrow=max.iter+1)
  lambdaA_ALL<-lambdaB_ALL<-obj.it<-rep(NA,length=max.iter)

  ### START CODE

  # Sequentially extract the r canonical vector pairs
  # for (i.r in 1:rank){

  # if (i.r==1){# for r=1: start from original data sets
  X_data<-X
  Y_data<-Y
  # }

  # STARTING VALUE
  if (init.method == "OLS"){ # OLS : better for nrow(X)>ncol(X) because initial values won't be sparse
    decomp<-eigen(cov(X_data))
    B.PRINCOMP<-X_data%*%decomp$vectors[,1]
    B.STARTING<-matrix(MASS::ginv(t(Y_data)%*%Y_data)%*%t(Y_data)%*%B.PRINCOMP,ncol=1)
    # B.STARTING<-apply(B.STARTING,2,NORMALIZATION_UNIT)
    A.STARTING<-matrix(decomp$vectors[,1],ncol=1)
  } else if (init.method == "SVD"){ # SVD
    decomp<-svd(cov(X,Y),nu=1,nv=1)
    A.STARTING <- matrix(decomp$u[,1],ncol=1)
    B.STARTING <- matrix(decomp$v[,1],ncol=1)
  } else { # LASSO : high-dimensional (when nrow(X)<ncol(X))
    hugefit<-huge(cov(X_data),lambda=1,verbose=F,cov.output=T,method="glasso")
    decomp<-try(eigs(as.matrix(hugefit$cov[[1]]),k=1))
    if(class(decomp)=="try-error"){
      decomp<-eigen(as.matrix(hugefit$cov[[1]]))
    }
    B.PRINCOMP<-X_data%*%decomp$vectors[,1]
    LASSOFIT_init<-glmnet(y=c(B.PRINCOMP),x=Y_data,family="gaussian",lambda=lambdaBseq)
    B.STARTING<-matrix(LASSOFIT_init$beta[,which(LASSOFIT_init$df!=0)],nrow=ncol(Y_data))
    if(all(B.STARTING==0)){
      LASSOFIT_init<-glmnet(y=c(B.PRINCOMP),x=Y_data,family="gaussian",lambda=10^-6)
      B.STARTING<-matrix(LASSOFIT_init$beta[,which(LASSOFIT_init$df!=0)],nrow=ncol(Y_data))
    }
    B_BIC<-apply(B.STARTING,2,BIC,Y.data=B.PRINCOMP,X.data=Y_data)
    B.STARTING<-matrix(B.STARTING[,which.min(B_BIC)],ncol=1)
    # B.STARTING<-apply(B.STARTING,2,NORMALIZATION_UNIT)
    A.STARTING<-matrix(decomp$vectors[,1],ncol=1)
  }

  # CONVERGENCE CRITERION
  obj.initial<-mean((X_data%*%matrix(A.STARTING,ncol=1)-Y_data%*%matrix(B.STARTING,ncol=1))^2)
  obj.it[1]<-obj.initial

  # INITIALIZE CONVERGENCE PARAMETERS
  it<-1
  diff.obj<-conv*10

  # FROM i until convergence: canonical vectors
  while( (it<max.iter) & (diff.obj>conv) ){

    # Estimating A conditional on B
    FIT.A<-Sparse.alternating(Xreg=X_data,Yreg=Y_data%*%B.STARTING,method=Xmethod,groupidx=X.groupidx)
    AHAT_FINAL<-FIT.A$COEF_FINAL
    lambdaA_ALL[it]<-FIT.A$LAMBDA_FINAL
    if(sum(AHAT_FINAL)==0) {
      BHAT_FINAL <- rep(0,ncol(Y)); break
    }

    # Estimating B conditional on A
    if(Ymethod=="SLR"){
      tmp.dat <- data.frame(Xtmp=X%*%AHAT_FINAL, Ytmp=Y)
      FIT.B<-lm(Xtmp~Y-1,tmp.dat)
      FIT.B$COEF_FINAL <- coef(FIT.B)
      FIT.B$LAMBDA_FINAL <- NA
    } else{
      FIT.B<-Sparse.alternating(Xreg=Y_data,Yreg=X_data%*%AHAT_FINAL,method=Ymethod,groupidx=Y.groupidx)
    }
    BHAT_FINAL<-FIT.B$COEF_FINAL
    if(sum(BHAT_FINAL)==0) {lambdaB_ALL[it]<-NA; break} else{lambdaB_ALL[it]<-FIT.B$LAMBDA_FINAL}

    # Check convergence
    obj.new<-mean((X_data%*%matrix(AHAT_FINAL,ncol=1)-Y_data%*%matrix(BHAT_FINAL,ncol=1))^2)
    obj.it[it+1]<-obj.new
    diff.obj<-abs(obj.new-obj.initial)/obj.initial

    # Updated starting values
    B.STARTING<-BHAT_FINAL
    A.STARTING<-AHAT_FINAL
    obj.initial<-obj.new
    it<-it+1
  } # end while-loop

  # # Number of ITERATIONS
  # iterations[i.r]<-it

  # CANONICAL VARIATES after convergence
  Uhat<-X_data%*%AHAT_FINAL
  Vhat<-Y_data%*%BHAT_FINAL


  # Express canonical vectors in terms of ORIGINAL DATA MATRICES
  # if (i.r==1){ # FIRST DIMENSION

  # # Final estimates of canonical vectors, variates and canonical correlation
  # ALPHA_ALL[,i.r]<-AHAT_FINAL
  # BETA_ALL[,i.r]<-BHAT_FINAL
  # U_ALL[,i.r]<-Uhat
  # V_ALL[,i.r]<-Vhat
  # cancors[1,i.r]<-abs(cor(Uhat,Vhat))

  # # Deflated data matrices
  # X_data<- round(X_data - Uhat%*%solve(t(Uhat)%*%Uhat)%*%t(Uhat)%*%X_data,10)
  # Y_data<- round(Y_data - Vhat%*%solve(t(Vhat)%*%Vhat)%*%t(Vhat)%*%Y_data,10)
  # Final estimates of canonical vectors, variates and canonical correlation
  if(sum(AHAT_FINAL)==0 | sum(BHAT_FINAL)==0) {
    cancors.spearman <- cancors.pearson <- NA
  } else{
    cancors.spearman<-cor(Uhat,Vhat,method="spearman")
    cancors.pearson<-cor(Uhat,Vhat,method="pearson")
  }

  # Final estimates of canonical vectors, variates and canonical correlation
  if(sum(AHAT_FINAL)==0 | sum(BHAT_FINAL)==0) {
    cancors.spearman <- cancors.pearson <- NA
  } else{
    cancors.spearman<-cor(Uhat,Vhat,method="spearman")
    cancors.pearson<-cor(Uhat,Vhat,method="pearson")
  }
  names(AHAT_FINAL) <- colnames(X)
  names(BHAT_FINAL) <- colnames(Y)

  # lambdaA_FINAL<-FIT.A$LAMBDA_FINAL
  # lambdaB_FINAL<-FIT.B$LAMBDA_FINAL
  lambdaA_FINAL<-lambdaA_ALL[it-1]
  lambdaB_FINAL<-lambdaB_ALL[it-1]
  # } else { # HIGHER ORDER DIMENSIONS
  #
  #   # A expressed in terms of original data set X
  #   FIT.Aorig<-Sparse.alternating(Yreg=Uhat,Xreg=X,lambdaseq=lambdaAseq)
  #   ALPHAhat<-FIT.Aorig$COEF_FINAL
  #   lambdaA_FINAL<-FIT.Aorig$LAMBDA_FINAL
  #
  #   # B expressed in terms of original data set Y
  #   FIT.Borig<-Sparse.alternating(Yreg=Vhat,Xreg=Y,lambdaseq=lambdaBseq)
  #   BETAhat<-FIT.Borig$COEF_FINAL
  #   lambdaB_FINAL<-FIT.Borig$LAMBDA_FINAL
  #
  #   # Final estimates of canonical vectors, variates and canonical correlation
  #   ALPHA_ALL[,i.r]<-ALPHAhat
  #   BETA_ALL[,i.r]<-BETAhat
  #   Uhat<-X%*%ALPHAhat
  #   Vhat<-Y%*%BETAhat
  #   U_ALL[,i.r]<-Uhat
  #   V_ALL[,i.r]<-Vhat
  #   cancors[1,i.r]<-abs(cor(Uhat,Vhat))
  #
  #   # Deflated data matrices: regress original data sets on all previously found canonical variates
  #   X_data<- X  - U_ALL[,1:i.r]%*%ginv(t(U_ALL[,1:i.r])%*%U_ALL[,1:i.r])%*%t(U_ALL[,1:i.r])%*%X
  #   Y_data<- Y -  V_ALL[,1:i.r]%*%ginv(t(V_ALL[,1:i.r])%*%V_ALL[,1:i.r])%*%t(V_ALL[,1:i.r])%*%Y
  # }

  # } # END FOR-LOOP

  ## OUTPUT
  out<-list(ALPHA=c(AHAT_FINAL),BETA=c(BHAT_FINAL),cancors.spearman=cancors.spearman,cancors.pearson=cancors.pearson,Uhat=Uhat,Vhat=Vhat,lambdaA=lambdaA_FINAL,lambdaB=lambdaB_FINAL,it=it-1)
}



