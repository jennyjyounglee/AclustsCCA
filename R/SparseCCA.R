### R-script containing Functions to perform Sparse CCA###

### Reference:
### SparseCCA developed by
### Wilms, Ines, and Christophe Croux. "Robust sparse canonical correlation analysis." BMC systems biology 10.1 (2016): 1-13.
###
### Downloaded From:
### https://sites.google.com/view/iwilms/software?authuser=0

SparseCCA<-function(X,Y,lambdaAseq=seq(from=1,to=0.01,by=-0.01),lambdaBseq=seq(from=1,to=0.01,by=-0.01),rank,max.iter=50,conv=5*10^-2){
  ### Function to perform Sparse Canonical Correlation Analysis using alternating regressions

  ### AUXILIARY FUNCTIONS
  # NORMALIZATION_UNIT
  # BIC
  # Sparse.alternating

  ### INPUT
  # X                   : (nxp) data matrix
  # Y                   : (nxq) data matrix
  # lambdaAseq          : grid of sparsity parameters for lambdaA
  # lambdaBseq          : grid of sparsity parameters for lambdaB
  # rank                : number of canonical vector pairs to extract
  # max.iter            : maximum number of iterations
  # conv                : tolerance value for convergence

  ### OUTPUT
  # ALPHA               : (pxr) estimated canonical vectors correpsonding to the first data set
  # BETA                : (qxr) estimated canonical vectors correpsonding to the second data set
  # cancors             : r estimated canonical correlations
  # U_ALL               : (nxr) estimated canonical variates corresponding to the first data set
  # V_ALL               : (nxr) estimated canonical variates corresponding to the second data set
  # lamdbaA             : value of the sparsity parameter lambdaA
  # lamdbaB             : value of the sparsity parameter lambdaB
  # it                  : number of iterations



  ### STORE RESULTS
  ALPHA_ALL<-matrix(NA,ncol=rank,nrow=ncol(X))
  BETA_ALL<-matrix(NA,ncol=rank,nrow=ncol(Y))
  U_ALL<-matrix(NA,ncol=rank,nrow=nrow(X))
  V_ALL<-matrix(NA,ncol=rank,nrow=nrow(Y))
  cancors<-matrix(NA,ncol=rank,nrow=1)
  lambdaA_ALL<-matrix(NA,ncol=rank,nrow=max.iter)
  lambdaB_ALL<-matrix(NA,ncol=rank,nrow=max.iter)
  iterations<-matrix(NA,ncol=rank,nrow=1)
  obj.it<-matrix(NA,ncol=rank,nrow=max.iter+1)

  ### START CODE

  # Sequentially extract the r canonical vector pairs
  for (i.r in 1:rank){

    if (i.r==1){# for r=1: start from original data sets
      X_data<-X
      Y_data<-Y
    }

    # STARTING VALUE
    if (nrow(X_data)>ncol(X_data)){ # OLS
      decomp<-eigen(cov(X_data))
      B.PRINCOMP<-X_data%*%decomp$vectors[,1]
      B.STARTING<-matrix(ginv(t(Y_data)%*%Y_data)%*%t(Y_data)%*%B.PRINCOMP,ncol=1)
      B.STARTING<-apply(B.STARTING,2,NORMALIZATION_UNIT)
      A.STARTING<-matrix(decomp$vectors[,1],ncol=1)
    } else { # LASSO
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
      B.STARTING<-apply(B.STARTING,2,NORMALIZATION_UNIT)
      A.STARTING<-matrix(decomp$vectors[,1],ncol=1)
    }

    # CONVERGENCE CRITERION
    obj.initial<-mean((X_data%*%matrix(A.STARTING,ncol=1)-Y_data%*%matrix(B.STARTING,ncol=1))^2)
    obj.it[1,i.r]<-obj.initial

    # INITIALIZE CONVERGENCE PARAMETERS
    it<-1
    diff.obj<-conv*10

    # FROM i until convergence: canonical vectors
    while( (it<max.iter) & (diff.obj>conv) ){

      # Estimating A conditional on B
      FIT.A<-Sparse.alternating(Xreg=X_data,Yreg=Y_data%*%B.STARTING,lambdaseq=lambdaAseq)
      AHAT_FINAL<-FIT.A$COEF_FINAL
      lambdaA_ALL[it,i.r]<-FIT.A$LAMBDA_FINAL

      # Estimating B conditional on A
      FIT.B<-Sparse.alternating(Xreg=Y_data,Yreg=X_data%*%AHAT_FINAL,lambdaseq=lambdaBseq)
      BHAT_FINAL<-FIT.B$COEF_FINAL
      lambdaB_ALL[it,i.r]<-FIT.B$LAMBDA_FINAL

      # Check convergence
      obj.new<-mean((X_data%*%matrix(AHAT_FINAL,ncol=1)-Y_data%*%matrix(BHAT_FINAL,ncol=1))^2)
      obj.it[it+1,i.r]<-obj.new
      diff.obj<-abs(obj.new-obj.initial)/obj.initial

      # Updated starting values
      B.STARTING<-BHAT_FINAL
      A.STARTING<-AHAT_FINAL
      obj.initial<-obj.new
      it<-it+1
    } # end while-loop

    # Number of ITERATIONS
    iterations[1,i.r]<-it

    # CANONICAL VARIATES after convergence
    Uhat<-X_data%*%AHAT_FINAL
    Vhat<-Y_data%*%BHAT_FINAL


    # Express canonical vectors in terms of ORIGINAL DATA MATRICES
    if (i.r==1){ # FIRST DIMENSION

      # Final estimates of canonical vectors, variates and canonical correlation
      ALPHA_ALL[,i.r]<-AHAT_FINAL
      BETA_ALL[,i.r]<-BHAT_FINAL
      U_ALL[,i.r]<-Uhat
      V_ALL[,i.r]<-Vhat
      cancors[1,i.r]<-abs(cor(Uhat,Vhat))

      # Deflated data matrices
      X_data<- round(X_data  - Uhat%*%solve(t(Uhat)%*%Uhat)%*%t(Uhat)%*%X_data,10)
      Y_data<- round(Y_data - Vhat%*%solve(t(Vhat)%*%Vhat)%*%t(Vhat)%*%Y_data,10)

      lambdaA_FINAL<-FIT.A$LAMBDA_FINAL
      lambdaB_FINAL<-FIT.B$LAMBDA_FINAL

    } else { # HIGHER ORDER DIMENSIONS

      # A expressed in terms of original data set X
      FIT.Aorig<-Sparse.alternating(Yreg=Uhat,Xreg=X,lambdaseq=lambdaAseq)
      ALPHAhat<-FIT.Aorig$COEF_FINAL
      lambdaA_FINAL<-FIT.Aorig$LAMBDA_FINAL

      # B expressed in terms of original data set Y
      FIT.Borig<-Sparse.alternating(Yreg=Vhat,Xreg=Y,lambdaseq=lambdaBseq)
      BETAhat<-FIT.Borig$COEF_FINAL
      lambdaB_FINAL<-FIT.Borig$LAMBDA_FINAL

      # Final estimates of canonical vectors, variates and canonical correlation
      ALPHA_ALL[,i.r]<-ALPHAhat
      BETA_ALL[,i.r]<-BETAhat
      Uhat<-X%*%ALPHAhat
      Vhat<-Y%*%BETAhat
      U_ALL[,i.r]<-Uhat
      V_ALL[,i.r]<-Vhat
      cancors[1,i.r]<-abs(cor(Uhat,Vhat))

      # Deflated data matrices: regress original data sets on all previously found canonical variates
      X_data<- X  - U_ALL[,1:i.r]%*%ginv(t(U_ALL[,1:i.r])%*%U_ALL[,1:i.r])%*%t(U_ALL[,1:i.r])%*%X
      Y_data<- Y -  V_ALL[,1:i.r]%*%ginv(t(V_ALL[,1:i.r])%*%V_ALL[,1:i.r])%*%t(V_ALL[,1:i.r])%*%Y
    }

  } # END FOR-LOOP

  ## OUTPUT
  out<-list(ALPHA=ALPHA_ALL,BETA=BETA_ALL,cancors=cancors,U_ALL=U_ALL,V_ALL=V_ALL,lambdaA=lambdaA_FINAL,lambdaB=lambdaB_FINAL,it=it)

}

### Sparse CCA AUXILIARY FUNCTIONS ###
Sparse.alternating<-function(Xreg,Yreg,lambdaseq){
  # AUX FUNCTION
  ### Function to perform sparse alternating regression

  ### INPUT
  #Xreg               : design matrix
  #Yreg               : response
  #lambdaseq          : sequence of sparsity parameters

  ### OUTPUT
  #COEF_FINAL         : estimated coefficients
  #LAMBDA_FINAL       : optimal sparsity parameter

  ## Standardize
  Xreg_st<-matrix(stdize(Xreg),ncol=ncol(Xreg))
  for (i.variable in 1:ncol(Xreg)){
    if (is.na(apply(Xreg_st,2,sum)[i.variable])==T) {
      Xreg_st[,i.variable]<-0}
  }

  ## LASSO Fit
  LASSOFIT<-glmnet(y=Yreg,x=Xreg_st,family="gaussian",lambda=lambdaseq,intercept=T)
  if (is.integer(which(LASSOFIT$df!=0)) && length(which(LASSOFIT$df!=0)) == 0L) {
    # Smaller lambda sequence necessary
    LASSOFIT<-glmnet(y=Yreg,x=Xreg_st,family="gaussian",intercept=T)
    COEFhat<-matrix(LASSOFIT$beta[,which(LASSOFIT$df!=0)[1:length(lambdaseq)]],nrow=ncol(Xreg_st)) # estimated coefficients
    LAMBDA<-LASSOFIT$lambda[which(LASSOFIT$df!=0)[1:length(lambdaseq)]] # lambda values
  } else {
    COEFhat<-matrix(LASSOFIT$beta[,which(LASSOFIT$df!=0)],nrow=ncol(Xreg_st)) # estimated coefficients
    LAMBDA<-LASSOFIT$lambda[which(LASSOFIT$df!=0)] # lambda values
  }

  # Selection of sparsity parameter
  BICvalues<-apply(COEFhat,2,BIC,Y.data=Yreg,X.data=Xreg_st) # BIC
  COEF_FINAL<-matrix(COEFhat[,which.min(BICvalues)],ncol=1) # Final coefficient estimates
  COEF_FINAL[which(apply(Xreg,2,sd)!=0),]<-COEF_FINAL[which(apply(Xreg,2,sd)!=0),]/apply(Xreg,2,sd)[which(apply(Xreg,2,sd)!=0)]
  COEF_FINAL<-apply(COEF_FINAL,2,NORMALIZATION_UNIT)
  LAMBDA_FINAL<-LAMBDA[which.min(BICvalues)]


  ## OUTPUT
  out<-list(COEF_FINAL=COEF_FINAL,LAMBDA_FINAL=LAMBDA_FINAL)
}

NORMALIZATION_UNIT<-function(U){
  # AUX FUNCTION
  # output: normalized vector U
  length.U<-as.numeric(sqrt(t(U)%*%U))
  if(length.U==0){length.U<-1}
  Uunit<-U/length.U
}

BIC<-function(U,Y.data,X.data){
  # AUX FUNCTION
  ### Calculate value of Bayesian information Criterion

  sigmahat<-sum(diag((1/nrow(Y.data))*(t(Y.data-X.data%*%U)%*%(Y.data-X.data%*%U))))
  BIC<-nrow(Y.data)*log(sigmahat) + length(which(U!=0))*log(nrow(Y.data))
  return(BIC)
}

### Sparse CCA with robust initalization ###
SparseCCA_Robinit<-function(X,Y,lambdaAseq=seq(from=1,to=0.01,by=-0.01),lambdaBseq=seq(from=1,to=0.01,by=-0.01),rank,max.iter=50,conv=5*10^-2,mode="fraction",ncores=1,tol=10^-8,alpha.sparselts=0.75){
  ### Function to perform Sparse Canonical Correlation Analysis using alternating regressions, with robust initialization

  ### INPUT
  # X                   : (nxp) data matrix
  # Y                   : (nxq) data matrix
  # lambdaAseq          : grid of sparsity parameters for lambdaA
  # lambdaBseq          : grid of sparsity parameters for lambdaB
  # rank                : number of canonical vector pairs to extract
  # max.iter            : maximum number of iterations
  # conv                : tolerance value for convergence
  # mode                : mode argument sparseLTS function
  # ncores              : number of cores to use
  # tol                 : tolerance argument sparse LTS
  # alpha.sparselts     : alpha argument sparse LTS

  ### OUTPUT
  # ALPHA               : (pxr) estimated canonical vectors correpsonding to the first data set
  # BETA                : (qxr) estimated canonical vectors correpsonding to the second data set
  # cancors             : r estimated canonical correlations
  # U_ALL               : (nxr) estimated canonical variates corresponding to the first data set
  # V_ALL               : (nxr) estimated canonical variates corresponding to the second data set
  # lamdbaA             : value of the sparsity parameter lambdaA
  # lamdbaB             : value of the sparsity parameter lambdaB
  # it                  : number of iterations


  ### STORE RESULTS
  ALPHA_ALL<-matrix(NA,ncol=rank,nrow=ncol(X))
  BETA_ALL<-matrix(NA,ncol=rank,nrow=ncol(Y))
  U_ALL<-matrix(NA,ncol=rank,nrow=nrow(X))
  V_ALL<-matrix(NA,ncol=rank,nrow=nrow(Y))
  cancors<-matrix(NA,ncol=rank,nrow=1)
  lambdaA_ALL<-matrix(NA,ncol=rank,nrow=max.iter)
  lambdaB_ALL<-matrix(NA,ncol=rank,nrow=max.iter)
  iterations<-matrix(NA,ncol=rank,nrow=1)
  obj.it<-matrix(NA,ncol=rank,nrow=max.iter+1)

  ### START CODE

  # Sequentially extract the r canonical vector pairs
  for (i.r in 1:rank){

    if (i.r==1){# for r=1: start from original data sets
      X_data<-X
      Y_data<-Y
    }

    decomp<-eigs(SCov(X_data),k=1)
    B.PRINCOMP<-X_data%*%decomp$vectors[,1]
    LASSOFIT_init<-sparseLTS(x=Y_data,y=c(B.PRINCOMP),lambda=0.001,mode = mode,intercept=F,normalize=T,initial="sparse", tol = tol,ncores=ncores,alpha=alpha.sparselts)
    B.STARTING<-matrix(LASSOFIT_init$coefficients,ncol=1)
    B.STARTING<-apply(B.STARTING,2,NORMALIZATION_UNIT)
    A.STARTING<-matrix(decomp$vectors[,1],ncol=1)

    # CONVERGENCE CRITERION
    obj.initial<-mean((X_data%*%matrix(A.STARTING,ncol=1)-Y_data%*%matrix(B.STARTING,ncol=1))^2)
    obj.it[1,i.r]<-obj.initial

    # INITIALIZE CONVERGENCE PARAMETERS
    it<-1
    diff.obj<-conv*10

    # FROM i until convergence: canonical vectors
    while( (it<max.iter) & (diff.obj>conv) ){

      # Estimating A conditional on B
      FIT.A<-Sparse.alternating(Xreg=X_data,Yreg=Y_data%*%B.STARTING,lambdaseq=lambdaAseq)
      AHAT_FINAL<-FIT.A$COEF_FINAL
      lambdaA_ALL[it,i.r]<-FIT.A$LAMBDA_FINAL

      # Estimating B conditional on A
      FIT.B<-Sparse.alternating(Xreg=Y_data,Yreg=X_data%*%AHAT_FINAL,lambdaseq=lambdaBseq)
      BHAT_FINAL<-FIT.B$COEF_FINAL
      lambdaB_ALL[it,i.r]<-FIT.B$LAMBDA_FINAL

      # Check convergence
      obj.new<-mean((X_data%*%matrix(AHAT_FINAL,ncol=1)-Y_data%*%matrix(BHAT_FINAL,ncol=1))^2)
      obj.it[it+1,i.r]<-obj.new
      diff.obj<-abs(obj.new-obj.initial)/obj.initial

      # Updated starting values
      B.STARTING<-BHAT_FINAL
      A.STARTING<-AHAT_FINAL
      obj.initial<-obj.new
      it<-it+1
    } # end while-loop

    # Number of ITERATIONS
    iterations[1,i.r]<-it

    # CANONICAL VARIATES after convergence
    Uhat<-X_data%*%AHAT_FINAL
    Vhat<-Y_data%*%BHAT_FINAL

    # Express canonical vectors in terms of ORIGINAL DATA MATRICES

    if (i.r==1){ # FIRST DIMENSION

      # Final estimates of canonical vectors, variates and canonical correlation
      ALPHA_ALL[,i.r]<-AHAT_FINAL
      BETA_ALL[,i.r]<-BHAT_FINAL
      U_ALL[,i.r]<-Uhat
      V_ALL[,i.r]<-Vhat
      cancors[1,i.r]<-abs(cor(Uhat,Vhat))

      # Deflated data matrices
      X_data<- round(X_data  - Uhat%*%solve(t(Uhat)%*%Uhat)%*%t(Uhat)%*%X_data,10)
      Y_data<- round(Y_data - Vhat%*%solve(t(Vhat)%*%Vhat)%*%t(Vhat)%*%Y_data,10)

      lambdaA_FINAL<-FIT.A$LAMBDA_FINAL
      lambdaB_FINAL<-FIT.B$LAMBDA_FINAL

    } else { # HIGHER ORDER DIMENSIONS

      # A expressed in terms of original data set X
      FIT.Aorig<-Sparse.alternating(Yreg=Uhat,Xreg=X,lambdaseq=lambdaAseq)
      ALPHAhat<-FIT.Aorig$COEF_FINAL
      lambdaA_FINAL<-FIT.Aorig$LAMBDA_FINAL

      # B expressed in terms of original data set Y
      FIT.Borig<-Sparse.alternating(Yreg=Vhat,Xreg=Y,lambdaseq=lambdaBseq)
      BETAhat<-FIT.Borig$COEF_FINAL
      lambdaB_FINAL<-FIT.Borig$LAMBDA_FINAL

      # Final estimates of canonical vectors, variates and canonical correlation
      ALPHA_ALL[,i.r]<-ALPHAhat
      BETA_ALL[,i.r]<-BETAhat
      Uhat<-X%*%ALPHAhat
      Vhat<-Y%*%BETAhat
      U_ALL[,i.r]<-Uhat
      V_ALL[,i.r]<-Vhat
      cancors[1,i.r]<-abs(cor(Uhat,Vhat))

      # Deflated data matrices: regress original data sets on all previously found canonical variates
      X_data<- X  - U_ALL[,1:i.r]%*%ginv(t(U_ALL[,1:i.r])%*%U_ALL[,1:i.r])%*%t(U_ALL[,1:i.r])%*%X
      Y_data<- Y -  V_ALL[,1:i.r]%*%ginv(t(V_ALL[,1:i.r])%*%V_ALL[,1:i.r])%*%t(V_ALL[,1:i.r])%*%Y
    }

  } # END FOR-LOOP

  ## OUTPUT
  out<-list(ALPHA=ALPHA_ALL,BETA=BETA_ALL,cancors=cancors,U_ALL=U_ALL,V_ALL=V_ALL,lambdaA=lambdaA_FINAL,lambdaB=lambdaB_FINAL,it=it)

}
