## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- warning=F, message=F, eval=F--------------------------------------------
#  devtools::install_github("https://github.com/jennyjyounglee/AclustsCCA")
#  library("AclustsCCA")

## ---- eval=F------------------------------------------------------------------
#  data(annot) # row: CpG sites
#  data.list <- generate.data(n=500)
#  
#  DATA.X <- data.list$DATA.X # row: subjects (n), column: Metals (p)
#  DATA.Y <- data.list$DATA.Y # row: subjects (n), column: CpG sites (q)
#  
#  
#  dist.type <- "spearman"
#  Aclust.method <-"average"
#  thresh.dist <-0.2
#  max.dist <-1000
#  bp.thresh.dist <-999
#  
#  Xmethod <- "lasso"
#  Ymethod <- "OLS"
#  maxB <- 300
#  nthread <- 2
#  
#  AclustsCCA.result <- AclustsCCA(X=DATA.X,
#                                  Y=DATA.Y,
#                                  clusters.list=NULL,
#                                  annot=annot,
#                                  # parameters for A-clustering
#                                  dist.type = dist.type,
#                                  Aclust.method = Aclust.method,
#                                  thresh.dist = thresh.dist,
#                                  max.dist = max.dist,
#                                  bp.thresh.dist = bp.thresh.dist,
#                                  # parameters for SparseCCA
#                                  Xmethod=Xmethod,
#                                  Ymethod=Ymethod,
#                                  # parameters for permutation test for AclustsCCA
#                                  maxB=maxB,
#                                  nthread=nthread)

## ----eval=F-------------------------------------------------------------------
#  # Implement A-clustering
#  all.clusters.list <- assign.to.clusters(betas = t(DATA.Y),
#                                          annot = annot,
#                                          dist.type = dist.type,
#                                          method = Aclust.method,
#                                          thresh.dist = thresh.dist,
#                                          bp.thresh.dist = bp.thresh.dist,
#                                          max.dist = max.dist)
#  # AclustsCCA only considers clusters with at least two probes
#  clusters.list <- all.clusters.list[sapply(all.clusters.list,length)!=1]

## ----eval=F-------------------------------------------------------------------
#  X.resid <- partial.residual(data=DATA.X,Z=Z,nthread=nthread)
#  Y.resid <- partial.residual(data=DATA.Y,Z=Z,nthread=nthread)

## ----eval=F-------------------------------------------------------------------
#  AclustsCCA.result <- AclustsCCA(X=X.resid,
#                                  Y=Y.resid,
#                                  clusters.list=clusters.list,
#                                  # parameters for SparseCCA
#                                  Xmethod=Xmethod,
#                                  Ymethod=Ymethod,
#                                  # parameters for permutation test for AclustsCCA
#                                  maxB=maxB,
#                                  permute=T,
#                                  nthread=3)
#  

