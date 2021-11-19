Sparse.alternating<-function(Xreg,Yreg,method,group.idx=NULL){
  # AUX FUNCTION
  ### Function to perform sparse alternating regression

  ### INPUT
  #Xreg               : design matrix
  #Yreg               : response
  #method             : sparse function
  #group.idx          : grouping index for group lasso or SGL

  ### OUTPUT
  #COEF_FINAL         : estimated coefficients
  #LAMBDA_FINAL       : optimal sparsity parameter

  ## LASSO Fit
  if(method == "lasso"){
    LASSOFIT<-glmnet(y=Yreg,x=Xreg_st,family="gaussian",intercept=F,standardize=F)
  # if (is.integer(which(LASSOFIT$df!=0)) && length(which(LASSOFIT$df!=0)) == 0L) {
  #   # Smaller lambda sequence necessary
  #   LASSOFIT<-glmnet(y=Yreg,x=Xreg_st,family="gaussian",intercept=T)
  #   COEFhat<-matrix(LASSOFIT$beta[,which(LASSOFIT$df!=0)[1:length(lambdaseq)]],nrow=ncol(Xreg_st)) # estimated coefficients
  #   LAMBDA<-LASSOFIT$lambda[which(LASSOFIT$df!=0)[1:length(lambdaseq)]] # lambda values
  # } else {
  #   COEFhat<-matrix(LASSOFIT$beta[,which(LASSOFIT$df!=0)],nrow=ncol(Xreg_st)) # estimated coefficients
  #   LAMBDA<-LASSOFIT$lambda[which(LASSOFIT$df!=0)] # lambda values
  # }
  } else if (method == "alasso"){
    ridge.cv.FIT <- cv.glmnet(x = Xreg,y = Yreg, type.measure = "mse",nfold = 5,alpha = 0,intercept=F,standardize=F)
    best_ridge_coef <- as.numeric(coef(ridge.cv.FIT, s = ridge.cv.FIT$lambda.min))[-1]
    FIT <- glmnet(x =Xreg,y =Yreg,alpha = 1,penalty.factor = 1 / abs(best_ridge_coef),intercept=F,standardize=F)
  } else if (method == "SGL"){
    FIT<-SGL(list(x=Xreg,y=Yreg), index = group.idx, type = "linear", nlam = 20,standardize=F) # ,intercept=F
    FIT$df <- apply(FIT$beta,2,function(x) sum(x !=0))
    if (is.integer(which(FIT$df!=0)) && length(which(FIT$df!=0)) == 0L){
      FIT<-SGL(list(x=Xreg,y=Yreg), index = group.idx, type = "linear", nlam = 100,standardize=F)}
    FIT$lambda <-FIT$lambdas
  } else if (method == "gglasso"){ # automatically standardize? group idx should be ordered. not working properly
    FIT<-gglasso(y=Yreg,x=Xreg,group=group.idx,loss="ls",intercept=F) # ,standardize=F
    if (is.integer(which(FIT$df!=0)) && length(which(FIT$df!=0)) == 0L){
      FIT<-gglasso(y=Yreg,x=Xreg,group=group.idx,loss="ls",nlambda=300)}
  } else{
    print("Condition is wrong")
  }
  COEFhat<-matrix(FIT$beta[,which(FIT$df!=0)],nrow=ncol(Xreg)) # estimated coefficients
  LAMBDA<-FIT$lambda[which(FIT$df!=0)] # lambda values

  # Selection of sparsity parameter
  BICvalues<-apply(COEFhat,2,BIC,Y.data=Yreg,X.data=Xreg) # BIC
  COEF_FINAL<-matrix(COEFhat[,BICvalues==min(BICvalues)],ncol=1) # Coefficient estimates: before scaling
  COEF_FINAL <- COEF_FINAL / sd(Xreg %*% COEF_FINAL) # Final coefficient estimates
  # COEF_FINAL<-matrix(COEFhat[,which.min(BICvalues)],ncol=1) # Final coefficient estimates
  # COEF_FINAL[which(apply(Xreg,2,sd)!=0),]<-COEF_FINAL[which(apply(Xreg,2,sd)!=0),]/apply(Xreg,2,sd)[which(apply(Xreg,2,sd)!=0)]
  # COEF_FINAL<-apply(COEF_FINAL,2,NORMALIZATION_UNIT)
  # LAMBDA_FINAL<-LAMBDA[which.min(BICvalues)]
  LAMBDA_FINAL<-LAMBDA[BICvalues==min(BICvalues)]
  BIC_FINAL <- which.min(BICvalues)


  ## OUTPUT
  out<-list(COEF_FINAL=COEF_FINAL,LAMBDA_FINAL=LAMBDA_FINAL,BIC=BIC_FINAL)
}
