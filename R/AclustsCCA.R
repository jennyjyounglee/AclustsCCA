
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
#' @param X` \eqn{n} by \eqn{p} exposure data matrix, where \eqn{n} is sample size and \eqn{p} is number of exposures.
#' @param Y` \eqn{n} by \eqn{q} outcome data matrix, where \eqn{n} is sample size and \eqn{q} is number of outcomes.
#' @param Z` \eqn{n} by \eqn{e} confounder data matrix, where \eqn{n} is sample size and \eqn{r} is number of confounders. If `NULL`, partial residuals are used for SparseCCA analysis.
#' @param clusters.list A list of clusters with CpG sites obtained using A-clustering, each item is a cluster that contains a set of probes. A-clustering is implemented if `NULL` or can be provided by users.
#' @param annot A preloaded annotation file that includes columns "IlmnID", "Coordinate_37", "Islands_Name", "Relation_to_Island", "UCSC_RefGene_Name". Only needed if \cite{clusters.list} is `NULL`.
#' @param dist.type A type of similarity distance function. Options are "spearman" (default), "pearson" (correlation measures) or "euclid".
#' @param Aclust.method A type of clustering function. Options are "single", "complete" or "average" (default).
#' @param thresh.dist A similarity distance threshold. Two neighboring clusters are merged to a single cluster if the similarity distance between them is above dist.thresh. Corresponds to $\bar{D}$ in the paper and the default is $0.2$
#' @param max.dist Optional maximum length between neighboring variables permitting to cluster them together. Corresponds to $\bar{d}_{bp}$ in the paper and the default is $1000$.
#' @param bp.thresh.dist A distance in chromosomal location. Any set of methylation sites within an interval smaller or equal to bp.dist will be potentially merged, depending on the similarity between sites at the ends of the interval. Corresponds to $\underline{d}_{bp}$ in the paper and the default is $999$.
#' @param Xmethod A penalty function for the exposure, i.e. penalty function when regressing Y onto X. Options are "lasso", "alasso","gglasso", and "SGL" (default).
#' @param Ymethod A penalty function for the outcome, i.e. penalty function when regressing X onto Y. Options are "lasso", "alasso","gglasso", "SGL", and "OLS" (default).
#' @param init.method         : Initialization method. Options are "lasso", "OLS", and "SVD" (default).
#' @param X.groupidx          : A vector of length \eqn{p} that indicates grouping structure of exposure \eqn{X}.
#' @param Y.groupidx          : A vector of length \eqn{q} that indicates grouping structure of outcome \eqn{Y}.
#' @param standardize         : A logical flag for exposure \eqn{X} and outcome \eqn{Y} standardization, prior to fitting the model.
#' @param max.iter            : A maximum number of iterations of SparseCCA. The default is $100$.
#' @param conv                : A tolerance value for convergence \eqn{epsilon} of SparseCCA. The default is $10e-2$.
#' @param maxnum A maximal total number of permutations across all the clusters.
#' @param maxB A maximal number of permutations for a single cluster.
#' @param permute.tmp.filepath Filepath to save intermittent permutation results.
#' @param permute A logical flag for whether to run permutation test or not.
#' @param nthread A number of threads to parallelize permutation test and implementation of SparseCCA across all the clusters.
#' @param FDR.thresh FDR threshold. The default is $0.05$.
#'
#' @return
#' The function returns a list of 6 objects according to the following order:
#'   - clusters.list          : A list of clusters with CpG sites obtained using A-clustering, each item is a cluster that contains a set of probes. If A-clustering is not implemented inside AclustsCCA, return NA.
#'   - ALPHA.observed         : A list of estimated canonical vector of length \eqn{p} corresponding to the exposure data \eqn{X} for each cluster
#'   - BETA.observed          : A list of estimated canonical vector of length \eqn{q} corresponding to the outcome data \eqn{Y} for each cluster
#'   - cancors.observed       : A vector of estimated canonical correlation for each cluster
#'   - permutation.result     : A \code{mmctest} object that contains permutation results
#'   - settings               : A settings used for the analysis
#'
#' @export
#'
#' @examples
#'
#'
AclustsCCA <- function(clusters.list=NULL,X,Y,Z=NULL,X.resid,Y.resid,annot=NULL,dist.type="spearman",Aclust.method="average",thresh.dist=0.2,max.dist=1000,bp.thresh.dist=999,Xmethod="lasso",Ymethod="OLS",standardize=T,X.groupidx=NULL,Y.groupidx=NULL,init.method="SVD",max.iter=100,conv=10^-2,maxnum=NULL,maxB=10000,FDR.thresh=0.05,permute.tmp.filepath=NULL,permute=T,nthread=2){
  if(is.null(clusters.list)){
    ##########################################################################
    # (1) Implement Aclustering
    ##########################################################################
    clusters.list <- Aclust::assign.to.clusters(betas=t(Y), annot=annot,
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
  # (2) Take partial residuals to take into account of confounders
  ##########################################################################
  if(!is.null(Z)){
    X.resid <- partial.residual(X,Z,nthread)
    Y.resid <- partial.residual(Y,Z,nthread)
  }
  if(is.null(Z) & (is.null(X.resid) | is.null(Y.resid))){
    X.resid <- X
    Y.resid <- Y
  }

  ##########################################################################
  # (3) Implement Sparse CCA on observed data
  ##########################################################################
  if(settings$nthread>1){ # Run in parallel

    cl <- parallel::makeCluster(settings$nthread)
    doSNOW::registerDoSNOW(cl)

    # Print out the progress for every iteration
    progress <- function(n) cat(sprintf("[Observed data: SparseCCA] cluster.idx = %d is complete\n", n))
    opts <- list(progress=progress)

    # Implement Sparse CCA on observed data
    AclustsCCA.observed <- foreach(cluster.idx = 1:length(clusters.list),.options.snow = opts,.packages = "AclustsCCA") %dopar% {
      Y.resid.subset <- Y.resid[, clusters.list[[cluster.idx]]]
      return(SparseCCA(X=X.resid,Y=Y.resid.subset,standardize=settings$standardize,Xmethod=settings$Xmethod,Ymethod=settings$Ymethod,X.groupidx=settings$X.groupidx,Y.groupidx=settings$Y.groupidx,init.method=settings$init.method,max.iter=settings$max.iter,conv=settings$conv))
    }
    parallel::stopCluster(cl)
  }else{  # Don't run in parallel
    AclustsCCA.observed <- lapply(1:length(clusters.list), function(cluster.idx) {
      cat("[Observed data: SparseCCA] cluster.idx = ", cluster.idx, "\n")
      Y.resid.subset <- Y[, clusters.list[[cluster.idx]]]
      return(SparseCCA(X=X.resid,Y=Y.resid.subset,standardize=settings$standardize,Xmethod=settings$Xmethod,Ymethod=settings$Ymethod,X.groupidx=settings$X.groupidx,Y.groupidx=settings$Y.groupidx,init.method=settings$init.method,max.iter=settings$max.iter,conv=settings$conv))
    })
}

  ALPHA.observed <- lapply(AclustsCCA.observed, function(x) x$ALPHA)
  BETA.observed <- lapply(AclustsCCA.observed, function(x) x$BETA)
  cancors.observed <- sapply(AclustsCCA.observed, function(x) x$cancors.spearman)

  ##########################################################################
  # (4) Implement Sparse CCA on permuted data
  ##########################################################################

  if(isTRUE(permute)){ # run permutation test
    # (3-1) Define sampler to run: test statistic as canonical correlation (spearman)
    data.sCCA <- list(X=X.resid,
                      Y=Y.resid.subset,
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

  if(is.null(annot)){ # clusters.list provided initially
    clusters.list <- NULL
  }

  # (3) Return
  return(list("clusters.list"=clusters.list,
              "ALPHA.observed"=ALPHA.observed,
              "BETA.observed"=BETA.observed,
              "cancors.observed"=cancors.observed,
              "permutation.result"=permute.result,
              "settings"=settings))
}

