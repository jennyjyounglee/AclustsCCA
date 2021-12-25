
#' @title
#' Implement additional permutations of AclustsCCA
#'
#' @description
#' Implement an iterative penalized least squares approach to
#' sparse canonical correlation analysis (SparseCCA)
#' with various penalty functions.
#'
#' @param obj A result of `AclustsCCA` function.
#' @param X \eqn{n} by \eqn{p} exposure data matrix, where \eqn{n} is sample size and \eqn{p} is number of exposures.
#' @param Y \eqn{n} by \eqn{q} outcome data matrix, where \eqn{n} is sample size and \eqn{q} is number of outcomes.
#' @param maxnum A maximal total number of permutations across all the clusters.
#' @param maxB A maximal number of permutations for a single cluster.
#'
#' @return
#' The function returns a list of 6 objects according to the following order:
#'   - clusters.list          : A list of clusters with CpG sites obtained using A-clustering, each item is a cluster that contains a set of probes. If A-clustering is not implemented inside AclustsCCA, return `NA`.
#'   - ALPHA.observed         : A list of estimated canonical vector of length \eqn{p} corresponding to the exposure data \eqn{X} for each cluster.
#'   - BETA.observed          : A list of estimated canonical vector of length \eqn{q} corresponding to the outcome data \eqn{Y} for each cluster.
#'   - cancors.observed       : A vector of estimated canonical correlation for each cluster.
#'   - permutation.result     : A \code{mmctest} object that contains permutation results.
#'   - settings               : A settings used for the analysis.
#'
#' @export
#'
#' @examples
#'
#'
#'
AclustsCCA.cont <- function(obj,X,Y,maxnum=NULL,maxB=10000){
  if(!obj$settings$permute){
    message("No permutation result found in AclustsCCA result Re-run AclustsCCA with permute=T")
    return();
  }
  if(is.null(maxnum)) maxnum <- length(obj$clusters.list)*1000
  if(obj$settings$maxB > maxB) message("maxB is smaller than previous maxB which is ",obj$settings$maxB,". Increase maxB.")
  if(obj$settings$maxnum > maxnum) message("maxnum is smaller than previous maxnum which is ",obj$settings$maxnum,". Increase maxnum.")

  sampler.result <- obj$permutation.result
  # (1) Define data to be used
  sampler.result@gensample@data$X <- X
  sampler.result@gensample@data$Y <- Y
  sampler.result@gensample@data$clusters.list <- obj$clusters.list
  sampler.result@gensample@data$TestStat.observed <- obj$cancors.observed

  # (2) Continue running
  sampler.result <- cont(sampler.result, steps=list(undecided=0,maxnum=maxnum,maxB=maxB))

  # (3) Return
  return(list("clusters.list"=obj$clusters.list,
              "ALPHA.observed"=obj$ALPHA.observed,
              "BETA.observed"=obj$BETA.observed,
              "cancors.observed"=obj$cancors.observed,
              "permutation.result"=sampler.result,
              "settings"=obj$settings))
}
