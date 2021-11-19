### R-script containing Functions to perform  CCA based on alternating regressions###

### Reference:
### CCA based on alternating regressions developed by
### Wilms, Ines, and Christophe Croux. "Robust sparse canonical correlation analysis." BMC systems biology 10.1 (2016): 1-13.
###
### Downloaded From:
### https://sites.google.com/view/iwilms/software?authuser=0


LS_alternating<-function(X,Y,rank,max.iter=5,conv=10^-2){
  ### Function to perform Standard Canonical Correlation Analysis using alternating regressions

  ### AUXILIARY FUNCTIONS
  # NORMALIZATION_UNIT

  ### INPUT
  # X                   : (nxp) data matrix
  # Y                   : (nxq) data matrix
  # rank                : number of canonical vector pairs to extract
  # max.iter            : maximum number of iterations
  # conv                : tolerance value for convergence

  ### OUTPUT
  # ALPHA               : (pxr) estimated canonical vectors correpsonding to the first data set
  # BETA                : (qxr) estimated canonical vectors correpsonding to the second data set
  # cancors             : r estimated canonical correlations
  # U_ALL               : (nxr) estimated canonical variates corresponding to the first data set
  # V_ALL               : (nxr) estimated canonical variates corresponding to the second data set
  # it                  : number of iterations

  ### STORE RESULTS
  ALPHA_ALL<-matrix(NA,ncol=rank,nrow=ncol(X))
  BETA_ALL<-matrix(NA,ncol=rank,nrow=ncol(Y))
  U_ALL<-matrix(NA,ncol=rank,nrow=nrow(X))
  V_ALL<-matrix(NA,ncol=rank,nrow=nrow(Y))
  cancors<-matrix(NA,ncol=rank,nrow=1)
  iterations<-matrix(NA,ncol=rank,nrow=1)
  obj.it<-matrix(NA,ncol=rank,nrow=max.iter+1)

  ### START CODE

  # Sequentially extract the r canonical vector pairs
  for (i.r in 1:rank){

    if (i.r==1){ # for r=1: start from original data sets
      X_data<-X
      Y_data<-Y
    }

    # STARTING VALUES
    B.PRINCOMP<-matrix(princomp(X_data)$scores[,1],ncol=1)
    B.STARTING<-matrix(ginv(t(Y_data)%*%Y_data)%*%t(Y_data)%*%B.PRINCOMP,ncol=1)
    B.STARTING<-apply(B.STARTING,2,NORMALIZATION_UNIT)
    A.PRINCOMP<-matrix(princomp(Y_data)$scores[,1],ncol=1)
    A.STARTING<-matrix(ginv(t(X_data)%*%X_data)%*%t(X_data)%*%A.PRINCOMP,ncol=1)
    A.STARTING<-apply(A.STARTING,2,NORMALIZATION_UNIT)

    # CONVERGENCE CRITERION
    obj.initial<-mean((X_data%*%matrix(A.STARTING,ncol=1)-Y_data%*%matrix(B.STARTING,ncol=1))^2)
    obj.it[1,i.r]<-obj.initial

    # INITIALIZE CONVERGENCE PARAMETERS
    it<-1
    diff.obj<-conv*10

    # FROM i until convergence: canonical vectors
    while( (it<max.iter) & (diff.obj>conv) ){

      # Estimating A conditional on B
      AHAT_FINAL<-matrix(ginv(t(X_data)%*%X_data)%*%t(X_data)%*%Y_data%*%B.STARTING,ncol=1)
      AHAT_FINAL<-apply(AHAT_FINAL,2,NORMALIZATION_UNIT)

      # Estimating B conditional on A
      BHAT_FINAL<-matrix(ginv(t(Y_data)%*%Y_data)%*%t(Y_data)%*%X_data%*%AHAT_FINAL,ncol=1)
      BHAT_FINAL<-apply(BHAT_FINAL,2,NORMALIZATION_UNIT)

      # Check convergence
      obj.new<-mean((X_data%*%matrix(AHAT_FINAL,ncol=1)-Y_data%*%matrix(BHAT_FINAL,ncol=1))^2)
      obj.it[it+1,i.r]<-obj.new
      diff.obj<-abs(obj.new-obj.initial)/obj.initial

      # updated starting values
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

    } else { # HIGHER ORDER DIMENSIONS

      # A expressed in terms of original data set X
      ALPHA_FINAL<-matrix(ginv(t(X)%*%X)%*%t(X)%*%Uhat,ncol=1)
      ALPHAhat<-apply(ALPHA_FINAL,2,NORMALIZATION_UNIT)

      # B expressed in terms of original data set Y
      BETA_FINAL<-matrix(ginv(t(Y)%*%Y)%*%t(Y)%*%Vhat,ncol=1)
      BETAhat<-apply(BETA_FINAL,2,NORMALIZATION_UNIT)

      # Final estimates of canonical vectors, variates, and canonical correlation
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
  out<-list(ALPHA=ALPHA_ALL,BETA=BETA_ALL,U_ALL=U_ALL,V_ALL=V_ALL,cancors=cancors,it=it)
}

### CCA AUXILIARY FUNCTIONS ###

NORMALIZATION_UNIT<-function(U){
  # AUX FUNCTION
  # output: normalized vector U
  length.U<-as.numeric(sqrt(t(U)%*%U))
  if(length.U==0){length.U<-1}
  Uunit<-U/length.U
}
