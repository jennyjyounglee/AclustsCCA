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
