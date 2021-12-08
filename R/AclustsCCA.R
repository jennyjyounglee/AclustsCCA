AclustsCCA <- function(clusters.list=NULL,X,Y,annot=NULL,dist.type="spearman",Aclust.method="average",thresh.dist=0.2,max.dist=1000,bp.thresh.dist=999,Xmethod="lasso",Ymethod="OLS",standardize=T,X.groupidx=NULL,Y.groupidx=NULL,init.method="SVD",max.iter=100,conv=10^-2,maxnum=NULL,maxB=10000,FDR.thresh=0.05,permute.tmp.filepath=NULL,permute=T,nthread=2){
  if(is.null(clusters.list)){
    ##########################################################################
    # (1) Implement Aclustering
    ##########################################################################
    clusters.list <- assign.to.clusters(betas=t(DATA.Y), annot=annot,
                                        dist.type = dist.type,
                                        method = Aclust.method,
                                        thresh.dist = thresh.dist,
                                        bp.thresh.dist = bp.thresh.dist,
                                        max.dist = max.dist)
    # We are only interested in clusters (non-sigletons)
    clusters.list <- clusters.list[sapply(clusters.list,length)!=1]
  }

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
                   maxB=maxB,
                   permute=permute,
                   nthread=nthread)
  ##########################################################################
  # (2) Implement Sparse CCA on observed data
  ##########################################################################
  if(settings$nthread>1){ # Run in parallel

    cl <- parallel::makeCluster(settings$nthread)
    doSNOW::registerDoSNOW(cl)

    # Print out the progress for every iteration
    progress <- function(n) cat(sprintf("[Observed data: SparseCCA] cluster.idx = %d is complete\n", n))
    opts <- list(progress=progress)

    # Implement Sparse CCA on observed data
    AclustsCCA.observed <- foreach(cluster.idx = 1:length(clusters.list),.options.snow = opts,.packages = "AclustsCCA") %dopar% {
      Y.subset <- Y[, clusters.list[[cluster.idx]]]
      return(SparseCCA(X=X,Y=Y.subset,standardize=settings$standardize,Xmethod=settings$Xmethod,Ymethod=settings$Ymethod,X.groupidx=settings$X.groupidx,Y.groupidx=settings$Y.groupidx,init.method=settings$init.method,max.iter=settings$max.iter,conv=settings$conv))
    }
    parallel::stopCluster(cl)
  }else{  # Don't run in parallel
    AclustsCCA.observed <- lapply(1:length(clusters.list), function(cluster.idx) {
      cat("[Observed data: SparseCCA] cluster.idx = ", cluster.idx, "\n")
      Y.subset <- Y[, clusters.list[[cluster.idx]]]
      return(SparseCCA(X=X,Y=Y.subset,standardize=settings$standardize,Xmethod=settings$Xmethod,Ymethod=settings$Ymethod,X.groupidx=settings$X.groupidx,Y.groupidx=settings$Y.groupidx,init.method=settings$init.method,max.iter=settings$max.iter,conv=settings$conv))
    })
}

  ALPHA.observed <- lapply(AclustsCCA.observed, function(x) x$ALPHA)
  BETA.observed <- lapply(AclustsCCA.observed, function(x) x$BETA)
  cancors.observed <- sapply(AclustsCCA.observed, function(x) x$cancors.spearman)

  ##########################################################################
  # (3) Implement Sparse CCA on permuted data
  ##########################################################################

  if(isTRUE(permute)){ # run permutation test
    # (3-1) Define sampler to run: test statistic as canonical correlation (spearman)
    data.sCCA <- list(X=X,Y=Y,
                      clusters.list=clusters.list,
                      TestStat.observed=cancors.observed,
                      settings=settings)
    AclustsCCASampler.run <- new("AclustsCCASampler", data=data.sCCA)
    alg <- mmctest(h=hBH, threshold=FDR.thresh)

    # (3-2) Run mcctest using AclustsCCASampler and define temporary file path
    permute.result <- run(alg, AclustsCCASampler.run, maxsteps=list(undecided=0,maxnum=maxnum,maxB=maxB))

    # (3-3) Remove unnecessary output to save memory
    permute.result@gensample@data$X <- NULL
    permute.result@gensample@data$Y <- NULL
    permute.result@gensample@data$clusters.list <- NULL
    permute.result@gensample@data$TestStat.observed <- NULL
  } else{ # don't run permutation test
    permute.result <- NA
  }

  # (3) Return
  return(list("clusters.list"=clusters.list,
              "ALPHA.observed"=ALPHA.observed,
              "BETA.observed"=BETA.observed,
              "cancors.observed"=cancors.observed,
              "permutation.result"=permute.result,
              "settings"=settings))
}

