AclustsCCA.cont <- function(AclustsCCA.result,X,Y,maxnum,maxB,permute.tmp.filepath){
  if(AclustsCCA.result$permute){
    message("No permutation result found in AclustsCCA function. Re-run AclustsCCA with permute=T")
    break;
  }

  sampler.result <- AclustsCCA.result$permute.result
  # (1) Define data to be used
  Y <- data.table(Y)
  sampler.result@gensample@data$X <- X
  sampler.result@gensample@data$Y <- Y
  sampler.result@gensample@data$clusters.list <- AclustsCCA.result$clusters.list
  sampler.result@gensample@data$TestStat.observed <- AclustsCCA.result$cancors.observed

  # (2) Continue running
  sampler.result <- cont(sampler.result, steps=list(undecided=0,maxnum=maxnum,maxB=maxB))

  # (3) Return
  return(list("clusters.list"=AclustsCCA.result$clusters.list,
              "ALPHA.observed"=AclustsCCA.result$ALPHA.observed,
              "BETA.observed"=AclustsCCA.result$BETA.observed,
              "cancors.observed"=AclustsCCA.result$cancors.observed,
              "permutation.result"=sampler.result,
              "settings"=AclustsCCA.result$settings))
}
