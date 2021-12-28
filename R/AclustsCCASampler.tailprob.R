
#' @title
#' Build AclustsCCA sampler based on canonical correlation
#'
#' @description
#' class AclustsCCASampler.tailprob, inherited from mmctSamplerGeneric
#' For more information: \url{https://cran.r-project.org/web/packages/simctest/vignettes/simctest-mmctest-intro.pdf}
#'
#' @include simctest.R mmctest.R
#'
#' @return
#'
############################################################################

# sCCA sampler
setClass("AclustsCCASampler.tailprob", contains="mmctSamplerGeneric",
         representation=representation(data="list"))

setMethod("getSamples", signature(obj="AclustsCCASampler.tailprob"),
          function(obj, ind, n) {
            X <- obj@data$X # n X p
            Y <- obj@data$Y # n X q
            num.clusters <- obj@data$num.clusters
            clusters.list <- obj@data$clusters.list # list of clusters obtained from Aclust
            TestStat.observed <- obj@data$TestStat.observed # Test statistics for permutation test
            settings <- obj@data$settings

            result <- matrix(NA,nrow=length(ind),ncol=max(n)) # number of hypothesis X number of permutations
            for(i in 1:length(ind)) { # for each hypothesis
              cat("[Permuted data: SparseCCA] running ",n[i],"permutations for hypothesis b = ",ind[i],", remaining hypothesis = ",length(ind),"\n")
              Y.subset <- Y[, clusters.list[[i]]] # n X q matrix

              if(settings$nthread>1 & (n[i] > 100)){ # Run in parallel
                cl <- parallel::makeCluster(settings$nthread)
                doParallel::registerDoParallel(cl)
                result[i,] <- foreach(b = 1:n[i], .combine = 'c') %dopar% {
                  index <- sample(1:nrow(X), size = nrow(X), replace = FALSE)
                  X_permute<-X[index,]
                  SparseCCA(X=X_permute,Y=Y.subset,standardize=settings$standardize,Xmethod=settings$Xmethod,Ymethod=settings$Ymethod,X.groupidx=settings$X.groupidx,init.method=settings$init.method,max.iter=settings$max.iter,conv=settings$conv)$tail.prob
                }
                parallel::stopCluster(cl)
              } else{  # Don't run in parallel
                for(b in 1:n[i]){
                  index <- sample(1:nrow(X), size = nrow(X), replace = FALSE)
                  X_permute<-X[index,]
                  result[i,b] <- SparseCCA(X=X_permute,Y=Y.subset,standardize=settings$standardize,Xmethod=settings$Xmethod,Ymethod=settings$Ymethod,X.groupidx=settings$X.groupidx,init.method=settings$init.method,max.iter=settings$max.iter,conv=settings$conv)$tail.prob
                }
              }
            }

            # RETURN: the number of exceeds for each hypothesis
            if(length(ind) <= 1){
              exceed <- sum(result < TestStat.observed[ind])
            }else{
              exceed <- apply(apply(result,2,function(x) 1*(x < TestStat.observed[ind])),1,function(y) sum(y,na.rm=T))
            }

            return(exceed)
          }
)


setMethod("getNumber", signature(obj="AclustsCCASampler.tailprob"),
          function(obj) {
            return(obj@data$num.clusters);
          }
)
