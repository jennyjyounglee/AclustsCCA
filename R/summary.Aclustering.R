summary.Aclustering <- function(clusters.list,annot,digits=2){
  DAT.clusters <- lapply(seq_along(clusters.list), function(i){
    data.frame("IlmnID"=clusters.list[[i]], "ClustIdx"=i)
  })
  DAT.clusters <- do.call("rbind",DAT.clusters)
  DAT.clusters$IlmnID <- as.character(DAT.clusters$IlmnID)
  DAT.clusters <- data.table(left_join(DAT.clusters,annot,by="IlmnID"))
  DAT.clusters[,"n.CGI":=uniqueN(UCSC_CpG_Islands_Name),by="ClustIdx"]

  clust.size.vec <- sapply(clusters.list,length)
  # (2-1) Total Number of Sites
  n.sites <- sum(clust.size.vec)
  # (2-2) Number of Singletons and Clusters (at leat two sites))
  n.clust <- length(clusters.list[clust.size.vec != 1])
  n.single <- sum( clust.size.vec == 1)
  # (2-3) Number of CpG sites in clusters (min, median, max)
  n.CpGs.in.clust <- quantile(clust.size.vec[clust.size.vec!=1],c(0,0.5,1))
  # (2-4) Base-pair distance between extremes (min, median, max)
  summary.dist <- DAT.clusters[,list("distance"=max(Coordinate_37)-min(Coordinate_37)),by=c("ClustIdx")]
  bp.dist <- quantile(summary.dist$distance,c(0,0.5,1))
  # (2-5) Clusters associated with a CpG resort: island, north/south shelf, or north/south shore within a single CGI
  CpG.resort <- c("Island", "N_Shelf", "N_Shore", "S_Shelf", "S_Shore")
  summary.resort <- DAT.clusters[n.CGI==1,list("n.resort"=sum(CpG.resort %in% Relation_to_UCSC_CpG_Island)),by=c("ClustIdx")]
  n.resort.1 <- summary.resort[,sum(n.resort==1)]  # Clusters associated with a predefined CpG Island
  n.resort.2 <- summary.resort[,sum(n.resort==2)]
  n.resort.3 <- summary.resort[,sum(n.resort==3)]
  n.resort.atleast1 <- summary.resort[,sum(n.resort>0)]
  n.resort.atleast2 <- summary.resort[,sum(n.resort>1)]


  ###########################################################################
  ##### (3) CREATE A TABLE OF SUMMARY OF ACLUSTERING RESULT
  ###########################################################################
  TABLE1 <-  data.frame(rbind(
    c("Total number of sites", n.sites),
    c("Number of singletons", n.single),
    c("Number of clusters (at least two sites)", n.clust),
    c("Number of CpG sites that are clustered ", paste0(sum(clust.size.vec[clust.size.vec!=1])," (",round(sum(clust.size.vec[clust.size.vec!=1])/n.sites * 100,digits),"%)",sep="")),
    c("Number of CpG sites in cluster (min, median, max)", paste0("(",paste0(n.CpGs.in.clust,collapse=","),")",sep="")),
    c("Base-pair distance between extremes (min, median, max) ", paste0("(",paste0(bp.dist,collapse=","),")",sep="")),
    c("(1) Clusters associated with a CpG resort : at least one resort regions", paste0(n.resort.atleast1," (",round(n.resort.atleast1/n.clust * 100,digits),"%)",sep="")),
    c("(2) Clusters associated with a CpG resort : at least two resort regions", paste0(n.resort.atleast2," (",round(n.resort.atleast2/n.clust * 100,digits),"%)",sep="")),
    c("(3) Clusters associated with a CpG resort : one resort regions", paste0(n.resort.1," (",round(n.resort.1/n.clust * 100,digits),"%)",sep="")),
    c("(4) Clusters associated with a CpG resort : two resort regions", paste0(n.resort.2," (",round(n.resort.2/n.clust * 100,digits),"%)",sep="")),
    c("(5) Clusters associated with a CpG resort : three resort regions", paste0(n.resort.3," (",round(n.resort.3/n.clust * 100,digits),"%)",sep=""))
  ))
  colnames(TABLE1) <- c("Characteristics","Quality")
  return(TABLE1)
}
