SparseCCA.permute.cont <- function(SparseCCA.permute.result,maxnum,maxB,Y,X,permute.tmp.filepath){
  sampler.result <- SparseCCA.permute.result$sampler.result

  # (1) Define data to be used
  Y <- data.table(Y)
  sampler.result@gensample@data$X <- X
  sampler.result@gensample@data$Y <- Y
  sampler.result@gensample@data$clusters.list <- SparseCCA.permute.result$clusters.list
  sampler.result@gensample@data$TestStat.observed <- SparseCCA.permute.result$cancors.observed
  # (2) Continue running
  sampler.result <- cont(sampler.result, steps=list(undecided=0,maxnum=maxnum,maxB=maxB))

  # (3) Return
  return(list("clusters.list"=SparseCCA.permute.result$clusters.list,
              "ALPHA.observed"=SparseCCA.permute.result$ALPHA.observed,
              "BETA.observed"=SparseCCA.permute.result$BETA.observed,
              "cancors.observed"=SparseCCA.permute.result$cancors.observed,
              "sampler.result"=sampler.result,
              "Xmethod"=SparseCCA.permute.result$Xmethod,
              "Ymethod"=SparseCCA.permute.result$Ymethod,
              "maxnum"=maxnum,
              "maxB"=maxB,
              "group.idx"=SparseCCA.permute.result$group.idx))
}
