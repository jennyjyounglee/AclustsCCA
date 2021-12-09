
#' @title
#' Build AclustsCCA sampler
#'
#' @description
#' class AclustsCCASampler, inherited from mmctSamplerGeneric
#' For more information: \url{https://cran.r-project.org/web/packages/simctest/vignettes/simctest-mmctest-intro.pdf}
#'
#' @include simctest.R mmctest.R
#'
#' @return
#'
############################################################################

# sCCA sampler
setClass("AclustsCCASampler", contains="mmctSamplerGeneric",
         representation=representation(data="list"))

setMethod("getSamples", signature(obj="AclustsCCASampler"),
          function(obj, ind, n) {
            X <- obj@data$X # n X p
            Y <- obj@data$Y # n X q
            clusters.list <- obj@data$clusters.list # list of clusters obtained from Aclust
            TestStat.observed <- obj@data$TestStat.observed # Test statistics for permutation test
            settings <- obj@data$settings

            result <- matrix(NA,nrow=length(ind),ncol=n[1]) # number of hypothesis X number of permutations
            for(b in 1:n[1]) { # for each permutation
              cat("[Permuted data: SparseCCA] permutation b = ",b,", remaining hypothesis = ",length(ind),"\n")
              index <- sample(1:nrow(X), size = nrow(X), replace = FALSE)
              X_permute<-X[index,]

              if(settings$nthread>1 & (length(ind) > 100)){ # Run in parallel
                cl <- parallel::makeCluster(settings$nthread)
                doParallel::registerDoParallel(cl)
                result[,b] <- foreach(i = 1:length(ind), .combine = 'c') %dopar% {
                  Y.subset <- Y[, clusters.list[[i]]] # n X q matrix
                  SparseCCA(X=X_permute,Y=Y.subset,standardize=settings$standardize,Xmethod=settings$Xmethod,Ymethod=settings$Ymethod,X.groupidx=settings$X.groupidx,Y.groupidx=settings$Y.groupidx,init.method=settings$init.method,max.iter=settings$max.iter,conv=settings$conv)$cancors.spearman
                }
                parallel::stopCluster(cl)
              } else{  # Don't run in parallel
                for(i in 1:length(ind)){
                  Y.subset <- Y[, clusters.list[[i]]] # n X q matrix
                  result[i,b] <- SparseCCA(X=X_permute,Y=Y.subset,standardize=settings$standardize,Xmethod=settings$Xmethod,Ymethod=settings$Ymethod,X.groupidx=settings$X.groupidx,Y.groupidx=settings$Y.groupidx,init.method=settings$init.method,max.iter=settings$max.iter,conv=settings$conv)$cancors.spearman
                }
              }
            }

            # RETURN: the number of exceeds for each hypothesis
            if(length(ind) <= 1){
              exceed <- sum(result > TestStat.observed[ind])
            }else{
              exceed <- apply(apply(result,2,function(x) 1*(x > TestStat.observed[ind])),1,sum)
            }

            return(exceed)
          }
)


setMethod("getNumber", signature(obj="AclustsCCASampler"),
          function(obj) {
            return(length(obj@data$clusters.list));
          }
)
