n <- 500
p <- 20
q <- 100

# Create autoregressive of order 1 matrix
AR.matrix <- function(n,rho){
  SIGMA<-matrix(0,ncol=n,nrow=n)
  for (i in 1:n){
    for (j in 1:n){
      SIGMA[i,j]<-rho^abs(i-j)
    }
  }
  return(SIGMA)
}

set.seed(12345)
# Load annotation file
annot <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annot <- data.table(data.frame(annot))
annot <- annot[!chr%in%c("chrY","chrX"),]
annot <- annot[,list("IlmnID"=Name,"Coordinate_37"=pos,"CHR"=as.numeric(gsub("chr", "", chr)),Islands_Name,Relation_to_Island,UCSC_RefGene_Name)]
annot <- annot[CHR==7,]
annot <- annot[order(CHR,Coordinate_37),]
annot <- annot[UCSC_RefGene_Name!="",]
annot <- annot[1:q,]

# Create clusters
annot[,"ClustIdx":= bumphunter::clusterMaker(annot$CHR, annot$Coordinate_37, maxGap = 1000)]
annot[,"ClustSize":=.N, by="ClustIdx"]
TRUE.ClustIdx <- sapply(c(2,3,4,7,12),function(x) (annot[ClustSize==x,][sample(.N,1),]$ClustIdx))
TRUE.CpGs <- annot[ClustIdx%in%TRUE.ClustIdx,IlmnID]
TRUE.Exposures <- sprintf("Exposure%d",c(1:3,9,10))

# Define correlation structure of betas
SIGMA.YY.data <- annot[,.("corr.mat"=list(AR.matrix(n=unique(ClustSize),rho=0.9))),by="ClustIdx"]
SIGMA.YY <- as.matrix(Matrix::bdiag(SIGMA.YY.data$corr.mat))
colnames(SIGMA.YY) <- rownames(SIGMA.YY) <-annot$IlmnID

# Define correlation structure of exposure
SIGMA.XX.data <- lapply(rep(0.7,p/5),function(x) AR.matrix(n=5,rho=x))
SIGMA.XX <- as.matrix(Matrix::bdiag(SIGMA.XX.data))
colnames(SIGMA.XX) <- rownames(SIGMA.XX) <- sprintf("Exposure%d",1:p)

# Define true elements
TRUE.ALPHA <- (rownames(SIGMA.XX) %in% TRUE.Exposures)*1
TRUE.BETA <- (colnames(SIGMA.YY) %in% TRUE.CpGs)*1
names(TRUE.ALPHA) <- rownames(SIGMA.XX)
names(TRUE.BETA) <- rownames(SIGMA.YY)

# Normalize loading
TRUE.ALPHA <- TRUE.ALPHA / as.numeric(sqrt(t(TRUE.ALPHA) %*% SIGMA.XX %*% TRUE.ALPHA))
TRUE.BETA<- TRUE.BETA / as.numeric(sqrt(t(TRUE.BETA) %*% SIGMA.YY %*% TRUE.BETA))

# DEFINE COVARIANCE MATRIX OF X AND Y (BASED ON CHEN'S PAPER)
SIGMA.XY <- 0.9 * SIGMA.XX %*% TRUE.ALPHA %*% t(TRUE.BETA) %*% SIGMA.YY

# CREATE COVARIANCE MATRIX FOR X AND Y
SIGMA_BIG<-rbind(cbind(SIGMA.XX,SIGMA.XY),cbind(t(SIGMA.XY),SIGMA.YY))

# Generate exposure and outcome data
DATA_clean<-mvtnorm::rmvnorm(floor(n),mean=rep(0,p+q),sigma=SIGMA_BIG)
DATA.X<-DATA_clean[,1:p]
DATA.Y<-DATA_clean[,(p+1):(p+q)]
DATA.Z <- data.frame("confounder1"=rnorm(n,0,0.0001),"confounder2"=rnorm(n))

colnames(DATA.X) <- colnames(SIGMA.XX)
colnames(DATA.Y) <- colnames(SIGMA.YY)

rownames(DATA.X) <- rownames(DATA.Y) <- rownames(DATA.Z)<- paste("Subject",1:n,sep="")

# Generate confounder data
confounding.effect <- as.matrix(DATA.Z) %*% matrix(c(0.01,0.01),nrow=2,ncol=1)
DATA.X <- DATA.X + matrix(rep(confounding.effect,p),ncol=p)
DATA.Y <- DATA.Y + matrix(rep(confounding.effect,q),ncol=q)

# Implement A-clustering
all.clusters.list <- Aclust::assign.to.clusters(betas = t(DATA.Y),
                                                annot = annot,
                                                dist.type = "spearman",
                                                method = "average",
                                                thresh.dist = 0.2,
                                                bp.thresh.dist = 999,
                                                max.dist = 1000)
clusters.list <- all.clusters.list[sapply(all.clusters.list,length)!=1]

sample.data <- list("DATA.X"=DATA.X,
                    "DATA.Y"=DATA.Y,
                    "DATA.Z"=DATA.Z,
                    "TRUE.CpGs"=TRUE.CpGs,
                    "TRUE.Exposures"=TRUE.Exposures,
                    "TRUE.ALPHA"=TRUE.ALPHA,
                    "TRUE.BETA"=TRUE.BETA,
                    "SIGMA.XX"=SIGMA.XX,
                    "SIGMA.YY"=SIGMA.YY,
                    "SIGMA.XY"=SIGMA.XY,
                    "clusters.list"=clusters.list)


usethis::use_data(sample.data, compress = "xz",overwrite = TRUE)
