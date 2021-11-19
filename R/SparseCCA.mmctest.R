#########################################################################################
# simctest.R

#Class sampalg
setClass("sampalg", representation(internal="environment"));

#Class sampalgPrecomp
setClass("sampalgPrecomp",contains="sampalg");

getalgprecomp <- function(level=0.05,epsilon=1e-3,halfspend=1000){
  sampalg <- new("sampalgPrecomp",internal=new.env());
  evalq({
    spending <-  function(n) {epsilon*n/(n+halfspend)};
  },envir=sampalg@internal)
  sampalg@internal$level <- level
  sampalg@internal$epsilon <- epsilon
  sampalg@internal$halfspend <- halfspend
  sampalg@internal$len <- 0;
  sampalg
}

setGeneric("getboundaryandprob",def=function(alg,n) {standardGeneric("getboundaryandprob")})
setMethod("getboundaryandprob", signature(alg="sampalgPrecomp"),
          function(alg,n) {
            extendbounds(alg,n);
            w <- (alg@internal$n<=n)
            list(U=alg@internal$U[1:n],
                 L=alg@internal$L[1:n],
                 n=alg@internal$n[w],
                 Sn=alg@internal$Sn[w],
                 prob=alg@internal$prob[w])
          })


setGeneric("extendbounds",def=function(alg,n) {standardGeneric("extendbounds")})
setMethod("extendbounds", signature(alg="sampalgPrecomp"),
          function(alg, n){
            if (is.null(alg@internal$U)){
              evalq({
                U <- c(2)
                L <- c(-1)
                n <- c()
                Sn <- c()
                prob <- c()
                preverr=c(0.,0.);
                len <- 1;
              },envir=alg@internal)
              assign("porig",c(1-get("level",envir=alg@internal),get("level",envir=alg@internal)), envir=alg@internal)
            }
            n <- as.integer(max(n,alg@internal$len))
            if (n>alg@internal$len){
              x <- .Call(simctest.extendboundsC, n=n-alg@internal$len,
                         level=alg@internal$level,
                         U=alg@internal$U[alg@internal$len],
                         L=alg@internal$L[alg@internal$len],
                         porig=alg@internal$porig,
                         preverr=alg@internal$preverr,
                         maxerr=alg@internal$spending((alg@internal$len+1):n),
                         returnstopbounds=as.integer(1)
              )
              alg@internal$U <- c(alg@internal$U,x$U)
              alg@internal$L <- c(alg@internal$L,x$L)
              alg@internal$n <- c(alg@internal$n,x$n+alg@internal$len+1)
              alg@internal$Sn <- c(alg@internal$Sn,x$Sn)
              alg@internal$prob <- c(alg@internal$prob,x$prob)
              alg@internal$preverr <- x$preverr
              alg@internal$porig <- x$porig
              alg@internal$len <- length(alg@internal$U)
            }
          }
)

setGeneric("getU",def=function(alg,ind) {standardGeneric("getU")})
setMethod("getU", signature(alg="sampalgPrecomp"),
          function(alg,ind) {
            extendbounds(alg,max(ind));alg@internal$U[ind]
          })
setGeneric("getL",def=function(alg,ind) {standardGeneric("getL")})
setMethod("getL", signature(alg="sampalgPrecomp"),
          function(alg,ind) {
            extendbounds(alg,max(ind));alg@internal$L[ind]
          })


setGeneric("run",def=function(alg,gensample,maxsteps=1e4) {standardGeneric("run")})
setMethod("run", signature(alg="sampalgPrecomp"),
          function(alg, gensample,maxsteps){
            data=new("sampalgres",
                     steps=0,pos=0,p.value=as.double(NA),alg=alg,gen=gensample
            )
            cont(data,maxsteps)
          }
)


#Class sampalgres
setClass("sampalgres", representation(p.value="numeric", steps="numeric", pos="numeric",
                                      alg="sampalg",gen="function"));




setGeneric("confint")
setMethod("confint", signature(object="sampalgres",parm="missing"),
          function(object,parm,level=0.95,...){
            if (is.na(object@p.value)) {
              #           warning("Algorithm has not stopped yet -  confidence interval
              #                    based on getbounds(...).\n");
              tmp <- getbounds(object);
              estup <- tmp[2]
              estlo <- tmp[1]
            }else
            {
              estup <- object@pos/object@steps
              estlo <- estup
            }


            nsteps=as.integer(max(400,object@steps+5*max(1/object@alg@internal$level,
                                                         1/(1-object@alg@internal$level))))
            while(1){
              x <- getboundaryandprob(object@alg, nsteps)
              if (length(x$U)-x$U[length(x$U)]>=3&&max(x$L)>=2)
                break;
              nsteps=as.integer(nsteps*1.5)
            }

            B.Sn<- x$Sn
            B.n <- x$n
            B.prob <- x$prob


            transf<-function(p,alpha,Sn,n,prob,target){
              #   exp/log is used to avoid problem with over/underruns
              x <- exp(log(prob)+ifelse(Sn>0,log(p/alpha)*Sn,0)
                       +ifelse(n-Sn>0,log((1-p)/(1-alpha))*(n-Sn),0))
              sum(sort(x))-target
            }
            #lower conf limit
            if (estlo>object@alg@internal$level){ #use upper branch
              target<-(1-level)/2;
              w<-( B.Sn/B.n)>=estlo
            } else {#use lower branch
              target<-1-(1-level)/2;
              w<-( B.Sn/B.n)<estlo
            }
            if (sum(w)==0)
              lower=0
            else
              lower=uniroot(transf,c(0,1),
                            alpha=object@alg@internal$level,
                            Sn=B.Sn[w],
                            n=B.n[w],
                            prob=B.prob[w],
                            target=target,tol=0.00001)$root
            #upper conf limit
            if (estup<object@alg@internal$level){ #use lower branch
              target<-(1-level)/2;
              w<-( B.Sn/B.n)<=estup
            }else{ #use upper branch
              target<-1-(1-level)/2;
              w<-( B.Sn/B.n)>estup
            }
            if (sum(w)==0)
              upper=1
            else
              upper=uniroot(transf,c(0,1),
                            alpha=object@alg@internal$level,
                            Sn=B.Sn[w],
                            n=B.n[w],
                            prob=B.prob[w],
                            target=target,tol=0.00001)$root
            res <- c(lower,upper);
            dim(res) <- c(1,2)
            colnames(res) <- paste(c((1-level)/2,1-(1-level)/2)*100,"%")
            rownames(res) <- "p.value"
            res
          }
)


setGeneric("cont",def=function(data,steps) {standardGeneric("cont")})
setMethod("cont", signature(data="sampalgres"),
          function(data,steps) {
            if (!is.na(data@p.value)) return(data);
            if (steps <=0){return(data)}
            pos <- data@pos
            if (length(formals(data@gen))==0){
              for (i in (data@steps+1):(data@steps+steps)){
                if (i>data@alg@internal$len) extendbounds(data@alg,max(data@alg@internal$len*1.5,100));
                pos <- pos+data@gen()
                if (pos>=data@alg@internal$U[i]){
                  data@p.value=pos/i;
                  data@steps=i
                  data@pos=pos
                  data@gen=function(){}
                  return(data)
                  #return(new("sampalgres",p.value=pos/i,steps=i,pos=pos,alg=data@alg));
                }
                if (pos<=data@alg@internal$L[i]){
                  data@p.value=pos/i;
                  data@steps=i;
                  data@pos=pos;
                  data@gen=function(){}
                  return(data)
                  #return(new("sampalgres",p.value=pos/i,steps=i,pos=pos,alg=data@alg));
                }
              }
            }else{            ## taking larger steps
              if (is.symbol(formals(data@gen)[[1]])){
                batchsize <- steps
              }else{
                batchsize <- formals(data@gen)[[1]]
              }
              i <- data@steps
              while (i < data@steps+steps){
                batchsize <- min(batchsize, data@steps+steps-i)
                if ((i+batchsize)>data@alg@internal$len) extendbounds(data@alg,max(i+batchsize,data@alg@internal$len*1.5,100));

                posvec <- pos+cumsum(data@gen(batchsize))
                U <- data@alg@internal$U[(i+1):(i+batchsize)]
                L <- data@alg@internal$L[(i+1):(i+batchsize)]
                stopind <- (U<=posvec|L>=posvec)
                if (any(stopind)){
                  stopind <- min(which(stopind))
                  i<- i+stopind
                  pos <- posvec[stopind]

                  data@p.value=pos/i
                  data@steps=i
                  data@pos=pos
                  data@gen=function(){}
                  return(data)
                  #return(new("sampalgres",p.value=pos/i,steps=i,pos=pos,alg=data@alg));
                }
                pos <- posvec[batchsize]
                i <- i+ batchsize
              }
            }
            data@p.value=as.double(NA)
            data@steps=i
            data@pos=pos
            return(data)
            #return(new("sampalgres",p.value=as.double(NA),steps=i,pos=pos,alg=data@alg,gen=data@gen));
          }
)

setGeneric("getbounds",def=function(data) {standardGeneric("getbounds")})
getindmax_getbounds <- function(data,start,Uest,Lest){
  alpha <- data@alg@internal$level
  epsilon <- data@alg@internal$epsilon
  k <-data@alg@internal$halfspend
  n <-data@steps+1
  target <- min(Uest-alpha,alpha - Lest)
  Delta <-function(x)sqrt(x*log((k+x)*(k+x-1)/(epsilon*k)))
  f<- function(x)(Delta(x)+1)/x-target
  xupper <- 2*n
  while (f(xupper)>0) xupper<-xupper*2
  indmax <- ceiling(uniroot(f,interval=c(n,xupper))$root)
  if (start<=1e6&&epsilon==1e-3&&alpha==0.05&&k==1000)
    indmax <- min(indmax,start+2000) #precomputed - see article
  indmax
}
getbounds_conservativebound <- function(data,n){
  alpha <- data@alg@internal$level
  epsilon <- data@alg@internal$epsilon
  k <-data@alg@internal$halfspend
  gamma <- (sqrt(-max(n)/2*log(epsilon*k/(k+max(n))/(k+max(n)-1)))+1)/max(n)
  c(alpha-gamma,alpha+gamma)
}
setMethod("getbounds", signature(data="sampalgres"),
          function(data){
            if (!is.na(data@p.value))
              return(c(data@p.value,data@p.value))
            else
            {
              n <-data@steps
              Sn <- data@pos
              #find first potential stopping times
              nteststart <- n+1
              while(TRUE){
                ntest <- max(ceiling(nteststart*1.2),nteststart+100)
                U <- getU(data@alg,nteststart:ntest)
                stoppos <- U<=Sn+nteststart:ntest-n
                if (any(stoppos)) break;
                nteststart <- ntest+1
              }
              kU <- min(which(stoppos))+nteststart-1;


              nteststart <- n+1
              while(TRUE){
                ntest <- max(ceiling(nteststart*1.2),nteststart+100)
                L <- getL(data@alg,nteststart:ntest)
                stoppos <- L>=Sn
                if (any(stoppos)) break;
                nteststart <- ntest+1
              }
              kL <- min(which(stoppos))+nteststart-1;


              indmaxorig <- getindmax_getbounds(data,max(kU,kL),getU(data@alg,kU-1)/kU,(getL(data@alg,kL-1)+1)/kL)
              indmax <- min(indmaxorig,max(2*n,n+1e5)) #bound on the number of steps looked forward
              indL <- kL:indmax
              indU <- kU:indmax
              if (indmax<indmaxorig)
                bcons <- getbounds_conservativebound(data,indmax)
              else
                bcons <- c(1,0)
              L <- getL(data@alg,c(indL-1,indmax))
              SL <- which(c(FALSE, L[-1]>L[-length(L)])) #indices where stopping is possible
              U <- getU(data@alg,c(indU-1,indmax))
              SU <- which(c(FALSE,U[-1]<=U[-length(U)])) # indicies wher stopping is possible
              return(c(max(0,min((L[SL-1]+1)/indL[SL-1],bcons[1])),max(U[SU-1]/indU[SU-1],bcons[2])))
            }
          }
)



setMethod("show", signature(object="sampalgres"),
          function(object){
            if (!is.na(object@p.value)){
              cat("p.value:",object@p.value,"\n")
            }
            else{
              cat("No decision reached.\n")
              b <- getbounds(object)
              cat("Final estimate will be in [",b[1], ",", b[2],"]\n")
              cat("Current estimate of the p.value:",object@pos/object@steps,"\n")
            }
            cat("Number of samples:",object@steps,"\n")
          }
)


#Class sampalgonthefly
setClass("sampalgonthefly", contains="sampalg");


getalgonthefly <- function(level=0.05,epsilon=1e-3,halfspend=1000){
  sampalg <- new("sampalgonthefly",internal=new.env());
  evalq({
    spending <-  function(n) {epsilon*n/(n+halfspend)};
  },envir=sampalg@internal)
  sampalg@internal$level <- level
  sampalg@internal$epsilon <- epsilon
  sampalg@internal$halfspend <- halfspend
  sampalg@internal$stepsize <- as.integer(500);
  sampalg
}


setMethod("getboundaryandprob", signature(alg="sampalgonthefly"),
          function(alg,n) {
            x <-.Call(simctest.extendboundsC, n=n-1,
                      level=as.double(alg@internal$level),
                      U=as.integer(2),
                      L=as.integer(-1),
                      porig=as.double(c(1-alg@internal$level,
                                        alg@internal$level)),
                      preverr=as.double(c(0.,0.)),
                      maxerr=alg@internal$spending(2:n),
                      returnstopprob=as.integer(1) )


            list(U=x$U,L=x$L,n=x$n+2,Sn=x$Sn,prob=x$prob)
          })

setMethod("run", signature(alg="sampalgonthefly"),
          function(alg, gensample,maxsteps){
            data=new("sampalgontheflyres",
                     U=2,L=-1,porig=c(1-alg@internal$level,alg@internal$level),ind=0,preverr=c(0,0),
                     steps=0,pos=0,p.value=as.double(NA),alg=alg,gen=gensample
            )
            cont(data,maxsteps)
          }
)


#Class sampalgontheflyres
setClass("sampalgontheflyres",
         representation(porig="numeric", U="numeric", L="numeric", ind="numeric", preverr="numeric"),
         contains="sampalgres");


setMethod("cont", signature(data="sampalgontheflyres"),
          function(data,steps){
            if (!is.na(data@p.value))return(data);
            if (steps<=0){return(data)}
            U <- data@U
            L <- data@L
            pos <- data@pos
            ind <- data@ind
            for (i in (data@steps+1):(data@steps+steps)){
              ind <- ind+1;
              if (ind >length(U)){
                x <- .Call(simctest.extendboundsC, n=data@alg@internal$stepsize,
                           level=as.double(data@alg@internal$level),
                           U=as.integer(U[length(U)]),
                           L=as.integer(L[length(L)]),
                           porig=as.double(data@porig),
                           preverr=as.double(data@preverr),
                           maxerr=data@alg@internal$spending((i):(i+data@alg@internal$stepsize)),
                           returnstopprob=as.integer(0)
                )
                U <- x$U
                L <- x$L
                data@preverr <- x$preverr
                data@porig <- x$porig
                ind <- 1;
              }
              pos <- pos+data@gen()
              if (pos>=U[ind]){
                data@p.value=pos/i;
                break;
              }
              if (pos<=L[ind]){
                data@p.value=pos/i;
                break;
              }
            }
            data@U <- U;
            data@L <- L;
            data@pos <- pos;
            data@ind <- ind;
            data@steps <- i
            return(data);
          })


setMethod("getbounds", signature(data="sampalgontheflyres"),
          function(data){
            if (!is.na(data@p.value))
              return(c(data@p.value,data@p.value))
            else
            { #not a very efficient implementation...
              Uorig <- c()
              Lorig <- c()
              norig <- c()
              nstart <- data@steps+length(data@U)-data@ind+1; #for further computations
              if (data@ind<length(data@U)){
                Uorig <- data@U[(data@ind+1):length(data@U)]
                Lorig <- data@L[(data@ind+1):length(data@L)]
                norig <- (data@steps+1):(nstart-1)
              }

              extendtoind <- function(ind){
                nsteps<-ind-nstart+1;
                if (nsteps>0){
                  x <- .Call(simctest.extendboundsC, n=nsteps,
                             level=as.double(data@alg@internal$level),
                             U=as.integer(data@U[length(data@U)]),
                             L=as.integer(data@L[length(data@L)]),
                             porig=as.double(data@porig),
                             preverr=as.double(data@preverr),
                             maxerr=data@alg@internal$spending((nstart):(nstart+nsteps)),
                             returnstopprob=as.integer(0)
                  )
                  U <- c(Uorig,x$U)
                  L <- c(Lorig,x$L)
                  n <- c(norig,nstart:(nstart+nsteps-1))
                  list(U=U,L=L,n=n)
                }
                else list(U=Uorig,L=Lorig,n=norig)
              }
              getU <- function(allind){
                x <- extendtoind(max(allind))
                x$U[allind-min(x$n)+1]
              }
              getL <- function(allind){
                x <- extendtoind(max(allind))
                x$L[allind-min(x$n)+1]
              }
              n <-data@steps
              Sn <- data@pos
              #find first potential stopping times
              nteststart <- n+1
              while(TRUE){
                ntest <- max(ceiling(nteststart*1.2),nteststart+100)
                U <- getU(nteststart:ntest)
                stoppos <- U<=Sn+nteststart:ntest-n
                if (any(stoppos)) break;
                nteststart <- ntest+1
              }
              kU <- min(which(stoppos))+nteststart-1;


              nteststart <- n+1
              while(TRUE){
                ntest <- max(ceiling(nteststart*1.2),nteststart+100)
                L <- getL(nteststart:ntest)
                stoppos <- L>=Sn
                if (any(stoppos)) break;
                nteststart <- ntest+1
              }
              kL <- min(which(stoppos))+nteststart-1;


              indmaxorig <- getindmax_getbounds(data,max(kU,kL),getU(kU-1)/kU,(getL(kL-1)+1)/kL)
              indmax <- min(indmaxorig,max(2*n,n+1e5)) #bound on the number of steps looked forward
              indL <- kL:indmax
              indU <- kU:indmax
              if (indmax<indmaxorig)
                bcons <- getbounds_conservativebound(data,indmax)
              else
                bcons <- c(1,0)
              L <- getL(c(indL-1,indmax))
              SL <- which(c(FALSE, L[-1]>L[-length(L)])) #indices where stopping is possible
              U <- getU(c(indU-1,indmax))
              SU <- which(c(FALSE,U[-1]<=U[-length(U)])) # indicies wher stopping is possible
              return(c(max(0,min((L[SL-1]+1)/indL[SL-1],bcons[1])),max(U[SU-1]/indU[SU-1],bcons[2])))
            }
          }
)


#A shortcut:
simctest <- function(gensample,level=0.05,epsilon=1e-3,maxsteps=1e4) run(getalgonthefly(level=level,epsilon=epsilon), gensample,maxsteps=maxsteps)


#########################################################################################

# mmctest -- Multiple Monte-Carlo Tests





# export: mmctest, hBH, hPC, hIndep
# some FDR controlling procedures

# FDR control by Benjamini-Hochberg
# This takes into account of ties
hBH <- function(p, threshold) {
  p.rank <- match(p, sort(unique(p)))
  sort.p.rank <- sort(p.rank)
  return(
    p.rank <= max( c(sort.p.rank[sort(p)<=(sort.p.rank*threshold/length(p))],-1))
  );
}

# Oridinal code
# hBH <- function(p, threshold) {
#
#   return(rank(p) <= max( c(which(sort(p)<=(1:length(p))*threshold/length(p)),-1) ));
# }

# FDR control by Pounds&Cheng
hPC <- function(p, threshold) {

  m <- length(p);
  pi0_hat <- min(1, 2/m*sum(p));
  return( rank(p) <= max( c(which(sort(p)<=(1:length(p))*threshold/ (pi0_hat*m) ),-1) ));
}

# independent testing via Bonferroni
hBonferroni <- function(p, threshold) {

  m <- length(p);
  return(p<=threshold/m);
}





# class mmctSamplerGeneric
# derive a class from mmctSamplerGeneric to implement the interface used to draw new samples
# or use class mmctSampler to directly pass a function and the number of hypotheses
setClass("mmctSamplerGeneric")
setGeneric("getSamples", def=function(obj, ind, n){standardGeneric("getSamples")})
setGeneric("getNumber", def=function(obj){standardGeneric("getNumber")})





# class mmctSampler
# directly pass a function "f" and the number of hypotheses "n", class has a slot for additional data
setClass("mmctSampler", contains="mmctSamplerGeneric", representation=representation(f="function", num="numeric", data="numeric"))

# implemented method
# for each hypotheses in vector ind[i], request n[i] new samples
setMethod("getSamples", signature(obj="mmctSampler"), function(obj, ind, n) {
  return(obj@f(ind, n, obj@data));
}
)

# implemented method
setMethod("getNumber", signature(obj="mmctSampler"), function(obj) {
  return(obj@num);
}
)

# pseudo constructor
mmctSampler <- function(f, num, data=NULL) {
  obj <- new("mmctSampler");
  obj@f <- f;
  obj@num <- num;
  obj@data <- data;
  return(obj);
}





# class mmctestres
# exportMethods: cont, show, pEstimate
setClass("mmctestres", contains="sampalgPrecomp", representation=representation(epsilon="numeric", threshold="numeric", r="numeric",batchN="numeric",pu="numeric",pl="numeric",m="numeric",
                                                                                h="function", gensample="mmctSamplerGeneric", g="numeric", num="numeric", A="numeric", B="numeric", C="numeric"))


setGeneric("mainalg", def=function(obj, stopcrit){standardGeneric("mainalg")})
setMethod("mainalg", signature(obj="mmctestres"), function(obj, stopcrit) {

  g <- obj@g;
  num <- obj@num;
  batchN <- obj@batchN;
  pu <- obj@pu;
  pl <- obj@pl;

  A <- obj@A;
  B <- obj@B;
  C <- obj@C;

  # m <- getNumber(obj@gensample);
  m <- obj@m
  currentBatch <- max(batchN);
  it <- 0;
  timer <- proc.time()[[3]];
  copying <- 0;

  tryCatch({
    while(length(B)>stopcrit$undecided) {
      cat("(1) # max permutations (max B) = ", max(num),"\n")
      # sample all \hat{p}_i in B
      currentBatch <- floor(1.25*currentBatch);
      if(currentBatch > 100000) { currentBatch <- 100000; }
      batchN[B] <- currentBatch;
      g[B] <- g[B] + getSamples(obj@gensample, B, currentBatch);
      num[B] <- num[B] + batchN[B];

      # compute upper "exact" confidence level (Clopper-Pearson)
      a <- (num/(num+obj@r) - (num-batchN)/(num-batchN+obj@r)) * (obj@epsilon/length(B));
      pu_ <- pu;
      pl_ <- pl;
      qindex <- (g>0) & (g<num);
      pu[qindex] <- 1 - qbeta(a[qindex]/2, num[qindex]-g[qindex], g[qindex]+1);
      pu[g==0] <- 1-(a[g==0]/2)**(1/num[g==0]);
      pu[g==num] <- 1;

      # ...and lower "exact" confidence level (Clopper-Pearson)
      pl[qindex] <- 1 - qbeta(1-a[qindex]/2, num[qindex]+1-g[qindex], g[qindex]);
      pl[g==0] <- 0;
      pl[g==num] <- (a[g==num]/2)**(1/num[g==num]);
      pu[setdiff(1:m,B)] <- pu_[setdiff(1:m,B)];
      pl[setdiff(1:m,B)] <- pl_[setdiff(1:m,B)];

      A <- which(obj@h(pu, obj@threshold));
      C <- which(obj@h(pl, obj@threshold));
      B <- setdiff(C, A);

      copying <- 1;
      obj@g <- g;
      obj@num <- num;
      obj@batchN <- batchN;
      obj@pu <- pu;
      obj@pl <- pl;

      obj@m <- m;
      obj@A <- A;
      obj@B <- B;
      obj@C <- C;
      copying <- 0;

      it <- it+1;
      if((stopcrit$maxnum>0) && (sum(num)+length(B)*floor(1.25*currentBatch)>=stopcrit$maxnum)) { break; }
      if((stopcrit$maxit>0) && (it>=stopcrit$maxit)) { break; }
      if((stopcrit$elapsedsec>0) && (proc.time()[[3]]-timer>=stopcrit$elapsedsec)) { break; }
      if((stopcrit$maxB>0) && (max(num)>=stopcrit$maxB)) { break; }

      if(!is.null(permute.tmp.filepath)){
        SAMPLER.TMP.RESULT <- obj
        SAMPLER.TMP.RESULT@gensample@data$X <- NULL
        SAMPLER.TMP.RESULT@gensample@data$Y <- NULL
        SAMPLER.TMP.RESULT@gensample@data$clusters.list <- NULL
        SAMPLER.TMP.RESULT@gensample@data$TailProb.observed <- NULL

        save(SAMPLER.TMP.RESULT,file=permute.tmp.filepath)
      }

    } # end while
  }, interrupt = function(interrupt){
    # finish copying if appropriate
    if(copying) {
      obj@g <- g;
      obj@num <- num;
      obj@batchN <- batchN;
      obj@pu <- pu;
      obj@pl <- pl;

      obj@m <- m;
      obj@A <- A;
      obj@B <- B;
      obj@C <- C;
      return(obj);
    }
  }) # end tryCatch

  return(obj);
}
)

# cont offers four stopping criteria given as list "steps"
# steps=list(maxit=0, maxnum=0, undecided=0, elapsedsec=0)
# corresonding to a maximal number of iterations, maximal number of samples drawn, number of yet undecided hypotheses, elapsed time
setMethod("cont", signature(data="mmctestres"), function(data, steps=list(maxit=0, maxnum=0, undecided=0, elapsedsec=0, maxB=0)) {

  if(length(setdiff(ls(steps),c("maxit","maxnum","undecided","elapsedsec","maxB")))>0) { stop("parameter steps contains invalid data"); }

  if(length(steps)==0) { stopcrit <- list(maxit=0, maxnum=0, undecided=0, elapsedsec=0, maxB=0); }
  else {
    stopcrit <- list();
    if("maxit" %in% ls(steps)) { stopcrit<-c(stopcrit,maxit=steps$maxit); } else { stopcrit<-c(stopcrit,maxit=0); }
    if("maxnum" %in% ls(steps)) { stopcrit<-c(stopcrit,maxnum=steps$maxnum); } else { stopcrit<-c(stopcrit,maxnum=0); }
    if("undecided" %in% ls(steps)) { stopcrit<-c(stopcrit,undecided=steps$undecided); } else { stopcrit<-c(stopcrit,undecided=0); }
    if("elapsedsec" %in% ls(steps)) { stopcrit<-c(stopcrit,elapsedsec=steps$elapsedsec); } else { stopcrit<-c(stopcrit,elapsedsec=0); }
    if("maxB" %in% ls(steps)) { stopcrit<-c(stopcrit,maxB=steps$maxB); } else { stopcrit<-c(stopcrit,maxB=0); }
  }

  return(mainalg(data, stopcrit));
}
)

setMethod("show", signature(object="mmctestres"), function(object) {

  # m <- getNumber(object@gensample);
  m <- object@m
  cat(paste("Number of rejected hypotheses: ",length(object@A),"\n",sep=""));
  cat(paste("Number of non-rejected hypotheses: ",length(setdiff(1:m, object@C)),"\n",sep=""));
  cat(paste("Number of unclassified hypotheses: ",length(object@B),"\n",sep=""));
  cat(paste("Total number of samples: ",sum(object@num),"\n",sep=""));
}
)

setGeneric("pEstimate", def=function(obj){standardGeneric("pEstimate")})
setMethod("pEstimate", signature(obj="mmctestres"), function(obj) {
  # return((obj@g+1)/(obj@num+1));
  return(obj@g/obj@num); # EDITED BY JENNY
}
)

setGeneric("confidenceLimits", def=function(obj){standardGeneric("confidenceLimits")})
setMethod("confidenceLimits", signature(obj="mmctestres"), function(obj) {

  # m <- getNumber(obj@gensample);
  m <- obj@m
  k <- 1000;
  g <- obj@g;
  num <- obj@num;
  batchN <- obj@batchN;
  pu <- rep(0, m);
  pl <- rep(0, m);

  # compute upper "exact" confidence level (Clopper-Pearson)
  a <- (num/(num+k) - (num-batchN)/(num-batchN+k)) * (1-(1-obj@epsilon)**(1/m));
  qindex <- (g>0) & (g<num);
  pu[qindex] <- 1 - qbeta(a[qindex]/2, num[qindex]-g[qindex], g[qindex]+1);
  pu[g==0] <- 1-(a[g==0]/2)**(1/num[g==0]);
  pu[g==num] <- 1;

  # ...and lower "exact" confidence level (Clopper-Pearson)
  pl[qindex] <- 1 - qbeta(1-a[qindex]/2, num[qindex]+1-g[qindex], g[qindex]);
  pl[g==0] <- 0;
  pl[g==num] <- (a[g==num]/2)**(1/num[g==num]);

  l <- list();
  l$lowerLimits <- pl;
  l$upperLimits <- pu;
  return(l);
}
)

setGeneric("testResult", def=function(obj){standardGeneric("testResult")})
setMethod("testResult", signature(obj="mmctestres"), function(obj) {

  # m <- getNumber(obj@gensample);
  m <- obj@m
  l <- list();
  l$rejected <- obj@A;
  l$nonrejected <- setdiff(1:m, obj@C);
  l$undecided <- obj@B;
  return(l);
}
)

setGeneric("summary.mmctestres", def=function(object,...){standardGeneric("summary.mmctestres")})
setMethod("summary.mmctestres", signature(object="mmctestres"), function(object) {
  # m <- getNumber(object@gensample)
  m <- object@m
  cat(paste("Number of hypotheses: ",m,sep=""));
  cat(strwrap(paste("Indices of rejected hypotheses:", paste(object@A,collapse=" ")),prefix="\n"));
  cat(strwrap(paste("Indices of unclassified hypotheses:", paste(object@B,collapse=" ")),prefix="\n"));
  cat("\nAll hypotheses not listed are classified as not rejected.\n");
}
)





# class mmctest
# exportMethods: run
# setClass("mmctest", contains="mmctestres", representation=representation(internal="environment"))

setMethod("run", signature(alg="mmctestres", gensample="mmctSamplerGeneric"), function(alg, gensample, maxsteps=list(maxit=0, maxnum=0, undecided=0, elapsedsec=0, maxB=0)) {

  obj <- new("mmctestres");

  obj@epsilon <- alg@epsilon;
  obj@threshold <- alg@threshold;
  obj@r <- alg@r;
  obj@h <- alg@h;
  obj@gensample <- gensample;

  # m <- getNumber(obj@gensample);
  m <-alg@m
  obj@g <- rep(0, m);
  obj@num <- rep(0, m);
  obj@batchN <- rep(10, m);
  obj@pu <- rep(0, m);
  obj@pl <- rep(0, m);

  obj@m <- m;
  obj@A <- 0;
  obj@B <- 1:m;
  obj@C <- 0;

  return(cont(obj, steps=maxsteps));
}
)


# pseudo constructor
mmctestres <- function(epsilon=0.01, threshold=0.1, r=10000, m, h) {

  obj <- new("mmctestres");
  obj@m=m;
  obj@epsilon=epsilon;
  obj@threshold=threshold;
  obj@r=r;
  obj@h=h;
  return(obj);
}


############################################################################
# sCCA sampler
setClass("SCCASampler.cancors", contains="mmctSamplerGeneric",
         representation=representation(data="list"))

setMethod("getSamples", signature(obj="SCCASampler.cancors"),
          function(obj, ind, n) {
            X <- obj@data$X # n X p
            Y <- obj@data$Y # n X q (data.table)
            clusters.list <- obj@data$clusters.list
            TestStat.observed <- obj@data$TestStat.observed

            result <- matrix(NA,nrow=length(ind),ncol=n[1]) # number of hypothesis X number of permutations
            for(b in 1:n[1]) { # for each permutation
              # cat("permutation b = ",b,", hypothesis = ",length(ind),"\n")
              index <- sample(1:nrow(X), size = nrow(X), replace = FALSE)
              X_permute<-X[index,]

              for(i in 1:length(ind)){
                Y.subset <- Y[, clusters.list[[i]], with=F] # n X q matrix
                result[i,b] <- SparseCCA(X_permute,as.matrix(Y.subset),standardize=T,method=Xmethod,Ymethod=Ymethod,group.idx=group.idx,init.method="SVD",max.iter=100,conv=10^-2)$cancors.spearman
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
