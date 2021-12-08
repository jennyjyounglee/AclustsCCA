SparseCCA.permute <- function(clusters.list,X,Y,Xmethod="lasso",Ymethod="OLS",standardize=T,X.groupidx=NULL,Y.groupidx=NULL,init.method="SVD",max.iter=100,conv=10^-2,maxnum=NULL,maxB=10000,FDR.thresh=0.05,permute.tmp.filepath=NULL){
  if(is.null(maxnum)) maxnum <- length(clusters.list)*1000

  settings <- list(Xmethod=Xmethod,
                   Ymethod=Ymethod,
                   standardize=standardize,
                   X.groupidx=X.groupidx,
                   Y.groupidx=Y.groupidx,
                   init.method=init.method,
                   max.iter=max.iter,
                   conv=conv,
                   maxnum=maxnum,
                   maxB=maxB)

  # (1) Implement Sparse CCA on observed data
  AclustsCCA.observed <- lapply(1:length(clusters.list), function(cluster.idx) {
    cat("[Observed data: SparseCCA] cluster.idx = ", cluster.idx, "\n")
    Y.subset <- Y[, clusters.list[[cluster.idx]]]
    return(SparseCCA(X=X,Y=Y.subset,standardize=settings$standardize,Xmethod=settings$Xmethod,Ymethod=settings$Ymethod,X.groupidx=settings$X.groupidx,Y.groupidx=settings$Y.groupidx,init.method=settings$init.method,max.iter=settings$max.iter,conv=settings$conv))
  })
  ALPHA.observed <- lapply(AclustsCCA.observed, function(x) x$ALPHA)
  BETA.observed <- lapply(AclustsCCA.observed, function(x) x$BETA)
  cancors.observed <- sapply(AclustsCCA.observed, function(x) x$cancors.spearman)

  # (2) Implement Sparse CCA on permuted data
  # (2-1) Define sampler to run: test statistic as canonical correlation (spearman)
  data.sCCA <- list(X=X,Y=Y,
                    clusters.list=clusters.list,
                    TestStat.observed=cancors.observed,
                    settings=settings)
  SCCASampler.run <- new("SCCASampler", data=data.sCCA)
  alg <- mmctest(h=hBH, threshold=FDR.thresh, m=length(clusters.list))

  # (2-2) Run mcctest using scca.sampler and define temporary file path
  sampler.result <- run(alg, SCCASampler.run, maxsteps=list(undecided=0,maxnum=maxnum,maxB=maxB))

  # (2-3) Remove unnecessary output
  sampler.result@gensample@data$X <- NULL
  sampler.result@gensample@data$Y <- NULL
  sampler.result@gensample@data$clusters.list <- NULL
  sampler.result@gensample@data$TestStat.observed <- NULL

  # (3) Return
  return(list("clusters.list"=clusters.list,
              "ALPHA.observed"=ALPHA.observed,
              "BETA.observed"=BETA.observed,
              "cancors.observed"=cancors.observed,
              "sampler.result"=sampler.result,
              "settings"=settings))
}

