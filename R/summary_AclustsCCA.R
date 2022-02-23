
#' @title
#' Outputs summary of AclustsCCA result
#'
#' @description
#' This function creates summary of AclustsCCA result table that summarizes A-clustering result.
#'
#' @param obj Output of **`AclustsCCA`**
#' @param annot A preloaded annotation file that includes columns the below name. For more information, \href{https://bioconductor.org/packages/release/data/annotation/html/IlluminaHumanMethylation450kanno.ilmn12.hg19.html}{IlluminaHumanMethylation450kanno.ilmn12.hg19}.
#'   - IlmnID: a character vector of IlmnID
#'   - CHR: a numeric vector of chromosome containing the CpG
#'   - Coordinate_37: a numeric vector of chromosomal coordinates of the CpG
#'   - Islands_Name: a character vector of chromosomal coordinates of the CpG Island
#'   - Relation_to_Island: a character vector of the location of the CpG relative to the CpG island
#'   - UCSC_RefGene_Name: a character vector of Target gene name(s)
#' @param n.top A number of top clusters with smallest p.value to show
#'
#' @return
#' The function returns a table of summary characteristics for each cluster:
#'   - ClustIdx: Index number of cluster corresponding to **`clusters.list`**
#'   - CHR: a numeric value of chromosome containing the CpG
#'   - UCSC_RefGene_Name: a list of target gene name(s)
#'   - Islands_Name: a list of chromosomal coordinates of the CpG Island
#'   - Significant: whether or not that cluster is significant
#'   - Cancors: a numeric value of canonical correlation
#'   - N.Probes: a number of probes in that cluster
#'   - DMR.length: a length of DMR in base pairs
#'   - Num.Exceed: a number of exceedness in permutation test
#'   - Num.Perm: a number of permutation
#'   - P.value: P.value corresponding to this cluster
#'   - CpGs: CpG sites with non-zero loadings
#'   - Exposures: Exposures with non-zero loadings
#'   - CpGs.Loading: A loading vector of CpGs sites with non-zero loadings
#'   - Exposures.Loading:  A loading vector of exposures with non-zero loadings
#'
#' @export
#'
#'
#'
summary_AclustsCCA <- function(obj,annot,n.top=NULL){
  sampler.result<-obj$permutation.result
  clusters.list<-obj$clusters.list
  ALPHA.observed<-obj$ALPHA.observed
  BETA.observed<-obj$BETA.observed
  cancors.observed<-obj$cancors.observed
  settings <- obj$settings

  rejected.cluster <- unique(c(sampler.result@A,
                               which(settings$h(pEstimate(sampler.result), settings$FDR.thresh))))
  num.clusters <- length(sampler.result@num)
  num.rejected <- length(rejected.cluster)
  num.perm <- sampler.result@num
  num.exceed <- sampler.result@g
  p.value <- pEstimate(sampler.result)

  # (2-0) Print the result
  # cat(paste("Total number of clusters: ",num.clusters,"\n",sep=""));
  # cat(paste("Number of rejected clusters: ",num.rejected,"\n",sep=""));
  # cat(paste("Maximum number of permutation among all clusters: ",max(sampler.result@num),"\n",sep=""));
  # cat(paste("Total number of permutation across all clusters: ",sum(sampler.result@num),"\n",sep=""));

  # (2-1) DEFINE TOP CLUSTERS WITH SAMLLEST PVALUE
  p.rank <- match(pEstimate(sampler.result), sort(unique(pEstimate(sampler.result)))) # 1 2 3 3 4 4 4 5 6 ...
  sort.p.rank <- sort(p.rank)
  order.p.rank <- order(p.rank)
  if(is.null(n.top)){ # if null, then extract significant results
    n.top <- num.rejected
    if(n.top==0){
      cat("None are significant")
      return(NULL)
    }
  }
  cluster.toprank <- order.p.rank[sort.p.rank <= n.top]


  # (2-2) Create Table
  nonzero.name <- function(x) names(x)[x!=0]
  nonzero.vector <- function(x) x[x!=0]

  TABLE2 <- lapply(cluster.toprank, function(x){
    data.table(
      "ClustIdx"=x,
      "CHR"=annot[IlmnID %in% nonzero.name(BETA.observed[[x]]),unique(CHR)],
      "UCSC_RefGene_Name"=list(annot[IlmnID %in% nonzero.name(BETA.observed[[x]]),unique(unlist(strsplit(UCSC_RefGene_Name,";")))]),
      "Islands_Name"=list(annot[IlmnID %in% nonzero.name(BETA.observed[[x]]),unique(unlist(strsplit(Islands_Name,";")))]),
      "Significant"=ifelse(x %in% rejected.cluster,"Yes","No"),
      "Cancors"=cancors.observed[[x]],
      "N.Probes"=annot[IlmnID %in% nonzero.name(BETA.observed[[x]]),.N],
      "DMR.Length"=annot[IlmnID %in% nonzero.name(BETA.observed[[x]]),max(Coordinate_37)-min(Coordinate_37)],
      "Num.Exceed"=num.exceed[x],
      "Num.Perm"=num.perm[x],
      "P.Value"=p.value[x],
      "CpGs"=list(nonzero.name(BETA.observed[[x]])),
      "Exposures"=list(nonzero.name(ALPHA.observed[[x]])),
      "CpGs.Loading"=list(nonzero.vector(BETA.observed[[x]])),
      "Exposures.Loading"=list(nonzero.vector(ALPHA.observed[[x]])))
  })
  TABLE2 <- do.call("rbind",TABLE2)
  TABLE2[,"rank":=match(P.Value, sort(unique(P.Value)))]

  return(data.frame(TABLE2))
}
