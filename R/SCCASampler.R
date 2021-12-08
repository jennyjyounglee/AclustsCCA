
#' @title
#' Build SCCA.SAMPLER
#'
#' @description
#' Build SCCA.SAMPLER
#'
#' @include simctest.R mmctest.R
#'
#' @param gen: function that performs one sampling step. Returns 0 (sampled test statistic does not exceed the observation) or 1 (sampled test static exceeds the observation)
#' @param method: which method to use for stopping
#'
#' @return
#' \code{mctest}, \code{mctest.RL} and \code{mctest.simctest} all return an object of class type \code{mctestres}, which has a print function (\code{print.mctestres}).
#' An object of class \code{mctestres} is a list with the following components: step (total batched number of samples drawn), decision.interval (interval for the p-value), decision (expressing significance), est.p (an estimate of the p-value) and realn (the actual number of samples taken without batching).
#'
#' @examples
#' Example use
############################################################################

# sCCA sampler
setClass("SCCASampler", contains="mmctSamplerGeneric",
         representation=representation(data="list"))

setMethod("getSamples", signature(obj="SCCASampler"),
          function(obj, ind, n) {
            X <- obj@data$X # n X p
            Y <- obj@data$Y # n X q
            clusters.list <- obj@data$clusters.list
            TestStat.observed <- obj@data$TestStat.observed
            settings <- obj@data$settings

            result <- matrix(NA,nrow=length(ind),ncol=n[1]) # number of hypothesis X number of permutations
            for(b in 1:n[1]) { # for each permutation
              # cat("permutation b = ",b,", hypothesis = ",length(ind),"\n")
              index <- sample(1:nrow(X), size = nrow(X), replace = FALSE)
              X_permute<-X[index,]

              for(i in 1:length(ind)){
                Y.subset <- Y[, clusters.list[[i]]] # n X q matrix
                result[i,b] <- SparseCCA(X=X_permute,Y=Y.subset,standardize=settings$standardize,Xmethod=settings$Xmethod,Ymethod=settings$Ymethod,X.groupidx=settings$X.groupidx,Y.groupidx=settings$Y.groupidx,init.method=settings$init.method,max.iter=settings$max.iter,conv=settings$conv)$cancors.spearman
              }
            }

            # RETURN: the number of exceedness for each hypothesis
            if(length(ind) <= 1){
              exceed <- sum(result > TestStat.observed[ind])
            }else{
              exceed <- apply(apply(result,2,function(x) 1*(x > TestStat.observed[ind])),1,sum)
            }

            return(exceed)
          }
)


setMethod("getNumber", signature(obj="SCCASampler"),
          function(obj) {
            return(length(obj@data$clusters.list));
          }
)
