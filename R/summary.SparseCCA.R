summary.SparseCCA <- function(obj,annot,n.top){
  sampler.result<-obj$sampler.result
  ALPHA.observed<-obj$ALPHA.observed
  BETA.observed<-obj$BETA.observed
  cancors.observed<-obj$cancors.observed

  # (2-1) DEFINE SIGNIFICANT CLUSTERS
  cluster.rejected <- unique(c(sampler.result@A, which(hBH(pEstimate(sampler.result), 0.05))))
  p.rank <- match(pEstimate(sampler.result), sort(unique(pEstimate(sampler.result)))) # 1 2 3 3 4 4 4 5 6 ...
  sort.p.rank <- sort(p.rank)
  order.p.rank <- order(p.rank)
  cluster.toprank<- order.p.rank[sort.p.rank <= n.top]

  # (2-2) WHICH PM25 AND CPGS ARE ASSOCIATED WITH THESE SIGNIFICANT CLUSTERS
  PM25 <- lapply(cluster.toprank,function(x) {
    names(ALPHA.observed[[x]])[ALPHA.observed[[x]]!=0]
  })
  CpGs <- lapply(cluster.toprank,function(x) {
    names(BETA.observed[[x]])[BETA.observed[[x]]!=0]
  })
  cancors <- sapply(cluster.toprank,function(x) {cancors.observed[[x]]})

  TABLE2 <- lapply(1:length(cluster.toprank), function(x){
    data.table(
      "Gene"=paste0(annot[IlmnID %in% CpGs[[x]],unique(unlist(strsplit(UCSC_RefGene_Name,";")))],collapse=","),
      "GeneGroup"=paste0(annot[IlmnID %in% CpGs[[x]],unique(unlist(strsplit(UCSC_RefGene_Group,";")))],collapse=","),
      "ClustIdx"=cluster.toprank[x],
      "CHR"=annot[IlmnID %in% CpGs[[x]],unique(CHR)],
      "PM25"=paste(PM25[[x]],collapse=","),
      "Loading"=paste(as.vector(ALPHA.observed[[cluster.toprank[x]]]),collapse=","),
      "cancors"=cancors[x],
      "n.sites"=annot[IlmnID %in% CpGs[[x]],.N],
      "DMR.length"=annot[IlmnID %in% CpGs[[x]],max(Coordinate_37)-min(Coordinate_37)],
      "num.exceed"=sampler.result@g[cluster.toprank[x]],
      "num.perm"=sampler.result@num[cluster.toprank[x]],
      "p.value"=sampler.result@g[cluster.toprank[x]]/sampler.result@num[cluster.toprank[x]],
      "num.hypothesis"=length(sampler.result@num),
      "num.undecided"=length(sampler.result@B),
      "CpGs"=paste0(CpGs[[x]],collapse=","))
  })
  TABLE2 <- do.call("rbind",TABLE2)
  TABLE2[ClustIdx %in% cluster.rejected,rejected:=1]
  TABLE2[,"rank":=match(p.value, sort(unique(p.value)))]

  return(data.frame(TABLE2))
}
