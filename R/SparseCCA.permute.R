SparseCCA.permute <- function(clusters.list,Y,X,Xmethod,Ymethod,group.idx,maxnum,maxB,permute.tmp.filepath){

  # (1) implement Sparse CCA on observed data
  AclustsCCA.observed <- lapply(1:length(clusters.list), function(cluster.idx) {
    cat("[Observed data: SparseCCA] cluster.idx = ", cluster.idx, "\n")
    Y.subset <- as.matrix(Y[, colnames(Y) %in% clusters.list[[cluster.idx]]])
    return(SparseCCA(X, Y.subset, standardize=T,
                     method=Xmethod,Ymethod=Ymethod,group.idx=group.idx,init.method="SVD",max.iter=100,conv=10^-2))
  })
  ALPHA.observed <- lapply(AclustsCCA.observed, function(x) x$ALPHA)
  BETA.observed <- lapply(AclustsCCA.observed, function(x) x$BETA)
  cancors.spearman.observed <- sapply(AclustsCCA.observed, function(x) x$cancors.spearman)

  # (2) Implement Sparse CCA on permuted data
  # (2-1) Define sampler to run: test statistic as canonical correlation (spearman)
  Y <- data.table(Y)
  data.sCCA <- list(X=X,Y=Y,clusters.list=clusters.list,TestStat.observed=cancors.spearman.observed)
  SCCASampler.run <- new("SCCASampler.cancors", data=data.sCCA)
  m <- mmctestres(h=hBH, threshold=0.05, m=length(clusters.list))

  # (2-2) Run mcctest using scca.sampler and define temporary filepath
  sampler.result <- run(m, SCCASampler.run, maxsteps=list(undecided=0,maxnum=maxnum,maxB=maxB))

  # (2-3) Remove unnecessary output
  sampler.result@gensample@data$X <- NULL
  sampler.result@gensample@data$Y <- NULL
  sampler.result@gensample@data$clusters.list <- NULL
  sampler.result@gensample@data$TestStat.observed <- NULL

  # (3) Return
  return(list("clusters.list"=clusters.list,
              "ALPHA.observed"=ALPHA.observed,
              "BETA.observed"=BETA.observed,
              "cancors.observed"=cancors.spearman.observed,
              "sampler.result"=sampler.result,
              "Xmethod"=Xmethod,
              "Ymethod"=Ymethod,
              "maxnum"=maxnum,
              "maxB"=maxB,
              "group.idx"=group.idx))
}

