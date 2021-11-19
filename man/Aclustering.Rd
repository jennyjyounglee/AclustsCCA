\name{Aclustering}
\alias{Aclustering}
\title{Implement Adjacent Clustering}
\usage{
assign.to.clusters(betas, annot = NULL, annotation.file.name = NULL, dist.type = "spearman", method = "single", dist.thresh = 0.5, bp.merge = 0, bp.thresh.clust = Inf, return.chroms = NULL)
}
\arguments{
  \item{Y}{An q X n matrix of DNA methylation values (either beta-values or M-values) where \eqn{q} is number of DNA methylation sites and \eqn{n} is number of subjects .}
  \item{annot}{A pre-loaded data frame containing annotation file. An annotation file is used in a function \code{\link[Aclust]{order.betas.by.chrom.location}} within \code{\link[Aclust]{assign.to.clusters}} to organize a matrix of methylation values to be ordered by chromosomal location. An annotation file should contain below three columns:
    \enumerate{
      \item IlmnID: Unique identifier from the Illumina CG database. (The probe ID).
      \item CHR: Chromosome containing the CpG (Build 37).
      \item MAPINFO: Chromosomal coordinates of the CpG (Build 37).
    }
  For example, I downloaded from Infinium Global Screening Array v1.0 Support Files (https://support.illumina.com/downloads/infinium-global-screening-array-v1-0-support-files.html).}
  \item{dist.type}{either "spearman" or "pearson". The deafult is "spearman".}
  \item{method}{either "simple","average",or "complete". The default is "average".}
  \item{dist.thresh}{A similarity distance threshold. Two neighboring clusters are merged to a single cluster if the similarity distance, based on correlation, between them is above dist.thresh. \eqn{\bar{D}} based on the original A-clustering paper. The default is \eqn{\bar{D}=0.2}}
  \item{bp.thresh.clust}{A maximum base pair difference between two neighboring sites permitting to cluster them together. \eqn{\bar{d}_{bp}} based on the orignal A-clustering paper. The default is \eqn{\bar{d}_{bp}=1000}.}
  \item{bp.merge}{A distance for dbp-merge initialization. Any set of methylation sites within an interval smaller or equal to \code{bp.merge} will be potentially merged at intiliazation, depending on the similarity between sites at the ends of the interval. \eqn{\underbar{d}_{bp}} based on the orignal A-clustering paper. The default is \eqn{\underbar{d}_{bp}=999}. }
}

\value{
Implements Adjacent Clustering (A-Clustering) based on DNA methylation values \code{Y} using annotation file \code{annot}.
}
\description{
Add together two numbers
}
\examples{
data(betas.7) ## upload methylation data
exposure <- rbinom(ncol(betas.7), 1,prob = 0.5) ## generate random exposure
covariates <- matrix(rnorm(2*ncol(betas.7)), ncol = 2)
### If wants to use Tim Triche's Bioconductor annotation package
# source("http://bioconductor.org/biocLite.R")
# biocLite("IlluminaHumanMethylation450k.db")
# require(IlluminaHumanMethylation450k.db)
# annot <- create.annot.triche(rownames(betas.7))

## else a version created on July 2013 is already in the package:
data(annot.7)

clusters.list<- assign.to.clusters(betas=betas.7, annot=annot.7,
                                     annotation.file.name = NULL,
                                     dist.type = "spearman",
                                     method = "average",
                                     dist.thresh = 0.2,
                                     bp.thresh.clust = 1000,
                                     bp.merge = 999,
                                     return.chroms = NULL)
}
