
#' @title
#' Summarize the A-clustering result
#'
#' @description
#' Implement an iterative penalized least squares approach to sparse canonical correlation analysis (SparseCCA)
#' with various penalty functions
#'
#'
#'
#' ### INPUT
#' @param clusters.list       : A list of clusters with CpG sites, each item is a cluster that contains a set of probes
#' @param annot               : A preloaded annotation file that includes columns the below name. For more information, \url{https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/methylationepic/infinium-methylationepic-manifest-column-headings.pdf}.
#'   - IlmnID: The IlmnID
#'   - Coordinate_37: Chromosomal coordinates of the CpG
#'   - Islands_Name: Chromosomal coordinates of the CpG Island
#'   - Relation_to_Island: The location of the CpG relative to the CpG island
#'   - UCSC_RefGene_Name: Target gene name(s)
#' @param digits              : Show rounded values to the specified number of decimal places
#'
#' @return
#' The function returns a table of summary characteristics:
#'   - Total number of sites used for the analysis
#'   - Number of singletons (a cluster with one CpG site)
#'   - Number of clusters (a cluster with at least two sites)
#'   - Number of CpG sites that are clustered
#'   - Number of CpG sites in cluster (min, median, max)
#'   - Base-pair distance between extremes (min, median, max)
#'   - (1) Clusters associated with a single refernece gene
#'   - (2) Clusters associated with a single CpG Island
#'   - (3-1) Clusters associated with a CpG resort : at least one resort regions
#'   - (3-2) Clusters associated with a CpG resort : at least two resort regions
#'   - (3-3) Clusters associated with a CpG resort : one resort regions
#'   - (3-4) Clusters associated with a CpG resort : two resort regions
#'   - (3-5) Clusters associated with a CpG resort : three resort regions
#'
#' @export
#'
#' @examples
#' library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
#' annot <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
#' annot <- data.table(data.frame(annot))
#' annot <- annot[!chr%in%c("chrY","chrX"),]
#' annot <- annot[,list("IlmnID"=Name,"Coordinate_37"=pos,"CHR"=as.numeric(gsub("chr", "", chr)),Islands_Name,Relation_to_Island,UCSC_RefGene_Name)]
#' annot <- annot[order(CHR,Coordinate_37),]
#'
#' clusters.list <- assign.to.clusters(betas.7, annot)
#' summary.Aclustering(clusters.list=clusters.list,annot=annot,digits=2)
#'
summary.Aclustering <- function(clusters.list,annot,digits=2){
  DAT.clusters <- lapply(seq_along(clusters.list), function(i){
    data.frame("IlmnID"=clusters.list[[i]], "ClustIdx"=i)
  })
  DAT.clusters <- do.call("rbind",DAT.clusters)
  DAT.clusters$IlmnID <- as.character(DAT.clusters$IlmnID)
  DAT.clusters <- data.table(dplyr::left_join(DAT.clusters,annot,by="IlmnID"))
  DAT.clusters[,"clust.size":=.N,by="ClustIdx"]
  DAT.clusters[,"n.CGI":=uniqueN(Islands_Name),by="ClustIdx"]

  # (2-1) Total Number of Sites
  n.sites <- DAT.clusters[,unique(clust.size),by="ClustIdx"][,sum(V1)]
  # (2-2) Number of Singletons and Clusters (at leat two sites)
  n.clust <- DAT.clusters[clust.size!=1,unique(clust.size),by="ClustIdx"][,uniqueN(ClustIdx)]
  n.single <- DAT.clusters[clust.size==1,unique(clust.size),by="ClustIdx"][,uniqueN(ClustIdx)]
  # (2-3) Number of CpG sites in clusters (min, median, max)
  n.CpGs.in.clust <- DAT.clusters[clust.size!=1,unique(clust.size),by="ClustIdx"]
  n.CpGs.in.clust.summary <- n.CpGs.in.clust[, quantile(V1,c(0,0.5,1))]
  # (2-4) Base-pair distance between extremes (min, median, max)
  summary.dist <- DAT.clusters[clust.size!=1,list("distance"=max(Coordinate_37)-min(Coordinate_37)),by=c("ClustIdx")]
  bp.dist <- quantile(summary.dist$distance,c(0,0.5,1))
  # (2-5) Clusters associated with one reference gene
  summary.gene <- DAT.clusters[clust.size!=1,list("n.gene"=uniqueN(UCSC_RefGene_Name)),by=c("ClustIdx")]
  n.gene.1 <- summary.gene[n.gene==1,.N]
  # (2-6) Clusters associated with one CpG Island
  summary.island <- DAT.clusters[clust.size!=1,list("n.CpGisland"=unique(n.CGI)),by=c("ClustIdx")]
  n.island.1 <- summary.island[n.CpGisland==1,.N]
  # (2-7) Clusters associated with a CpG resort: island, north/south shelf, or north/south shore within a single CGI
  CpG.resort <- c("Island", "N_Shelf", "N_Shore", "S_Shelf", "S_Shore")
  summary.resort <- DAT.clusters[n.CGI==1 & clust.size!=1,list("n.resort"=sum(CpG.resort %in% Relation_to_Island)),by=c("ClustIdx")]
  n.resort.1 <- summary.resort[n.resort==1,.N]  # Clusters associated with a predefined CpG Island
  n.resort.2 <- summary.resort[n.resort==2,.N]
  n.resort.3 <- summary.resort[n.resort==3,.N]
  n.resort.atleast1 <- summary.resort[n.resort>0,.N]
  n.resort.atleast2 <- summary.resort[n.resort>1,.N]


  ###########################################################################
  ##### (3) CREATE A TABLE OF SUMMARY OF ACLUSTERING RESULT
  ###########################################################################
  TABLE1 <-  data.frame(rbind(
    c("Total number of sites", n.sites),
    c("Number of singletons", n.single),
    c("Number of clusters (at least two sites)", n.clust),
    c("Number of CpG sites that are clustered ", paste0(sum(n.CpGs.in.clust$V1)," (",round(sum(n.CpGs.in.clust$V1)/n.sites * 100,digits),"%)",sep="")),
    c("Number of CpG sites in cluster (min, median, max)", paste0("(",paste0(n.CpGs.in.clust.summary,collapse=","),")",sep="")),
    c("Base-pair distance between extremes (min, median, max) ", paste0("(",paste0(bp.dist,collapse=","),")",sep="")),
    c("(1) Clusters associated with a single refernece gene ", paste0(n.gene.1,"/",n.clust," (",round(n.gene.1/n.clust * 100,digits),"%)",sep="")),
    c("(2) Clusters associated with a single CpG Island ", paste0(n.island.1,"/",n.clust," (",round(n.island.1/n.clust * 100,digits),"%)",sep="")),
    c("(3-1) Clusters associated with a CpG resort : at least one resort regions", paste0(n.resort.atleast1," (",round(n.resort.atleast1/n.clust * 100,digits),"%)",sep="")),
    c("(3-2) Clusters associated with a CpG resort : at least two resort regions", paste0(n.resort.atleast2," (",round(n.resort.atleast2/n.clust * 100,digits),"%)",sep="")),
    c("(3-3) Clusters associated with a CpG resort : one resort regions", paste0(n.resort.1," (",round(n.resort.1/n.clust * 100,digits),"%)",sep="")),
    c("(3-4) Clusters associated with a CpG resort : two resort regions", paste0(n.resort.2," (",round(n.resort.2/n.clust * 100,digits),"%)",sep="")),
    c("(3-5) Clusters associated with a CpG resort : three resort regions", paste0(n.resort.3," (",round(n.resort.3/n.clust * 100,digits),"%)",sep=""))
  ))
  colnames(TABLE1) <- c("Characteristics","Quality")
  return(TABLE1)
}
