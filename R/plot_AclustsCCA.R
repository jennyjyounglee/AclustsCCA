
#' @title
#' Creates a heatmap of significant CpG sites/Genes and exposures
#'
#' @description
#' This function creates a heatmap of significant CpG sites/Genes and exposures
#'
#' @param obj Output of **`AclustsCCA`**
#' @param annot A preloaded annotation file that includes columns the below name. For more information, \href{https://bioconductor.org/packages/release/data/annotation/html/IlluminaHumanMethylation450kanno.ilmn12.hg19.html}{IlluminaHumanMethylation450kanno.ilmn12.hg19}.
#'   - IlmnID: a character vector of IlmnID
#'   - CHR: a numeric vector of chromosome containing the CpG
#'   - Coordinate_37: a numeric vector of chromosomal coordinates of the CpG
#'   - Islands_Name: a character vector of chromosomal coordinates of the CpG Island
#'   - Relation_to_Island: a character vector of the location of the CpG relative to the CpG island
#'   - UCSC_RefGene_Name: a character vector of Target gene name(s)
#' @param abs A logical flag for showing absolute value of loadings or not
#'
#' @return
#' A heatmap of significant CpG sites/Genes and exposures
#' where columns are CpG sites/Genes and rows are exposures.
#' Selected exposures are colored otherwise remains white.
#'
#' @import ggplot2 reshape2
#' @export
#'
#'
#'

plot_AclustsCCA <- function(obj,annot,abs=TRUE){
  sampler.result<-obj$permutation.result
  clusters.list<-obj$clusters.list
  ALPHA.observed<-obj$ALPHA.observed
  BETA.observed<-obj$BETA.observed
  cancors.observed<-obj$cancors.observed

  TABLE.AclustsCCA <- summary_AclustsCCA(obj=AclustsCCA.result,annot=annot)
  if(is.null(TABLE.AclustsCCA)){
    cat("None are significant")
    return(NULL)
  }
  TABLE.AclustsCCA <- TABLE.AclustsCCA[TABLE.AclustsCCA$Significant=="Yes",]

  # (2-1) DEFINE SIGNIFICANT CLUSTERS
  rejected.cluster <- rev(TABLE.AclustsCCA$ClustIdx)

  # (2-2) WHICH PM25 AND CPGS ARE ASSOCIATED WITH THESE SIGNIFICANT CLUSTERS
  CpGs <- lapply(rejected.cluster,function(x) {names(BETA.observed[[x]])[BETA.observed[[x]]!=0]})
  Genes <- sapply(1:length(rejected.cluster),function(i) paste0(annot[IlmnID %in% CpGs[[i]],unique(unlist(strsplit(UCSC_RefGene_Name,";")))],collapse=","))
  Genes <- sapply(Genes,function(x) sub("\\-.*","",x))
  PM25 <- sapply(rejected.cluster,function(x) ALPHA.observed[[x]])
  colnames(PM25) <- Genes
  tmp <- data.frame(t(PM25))
  PM25 <- tmp[order(abs(tmp$Ni),abs(tmp$V),abs(tmp$Zn),decreasing=T),]
  PM25 <- t(PM25)
  if(isTRUE(abs)){
    PM25 <- abs(PM25)
    xmin <- round(min(PM25),digits=1)
  }else{
    xmin <- round(min(PM25)-0.1,digits=1)
  }
  xmax <- round(max(PM25)+0.1,digits=1)

  # (2-3) Create a heatmap
  DAT.plot <- reshape2::melt(PM25)
  # DAT.plot$value <- abs(DAT.plot$value)
  PLOT.HEATMAP <- ggplot2::ggplot(data = DAT.plot, aes(x=Var1, y=Var2, fill=value))+
    ggplot2::ggtitle("Heatmap of Significant DMRs and PM2.5 Components")+
    ggplot2::geom_tile(color = "gray", size=1.5)+
    ggplot2::theme_bw() +
    ggplot2::scale_fill_gradient2(low = "blue", high = "red",
                                  limit = c(xmin,xmax), space = "Lab",
                                  name="Loadings") +
    ggplot2::theme(plot.title = element_text(hjust = 0.5, size = 20),
                   axis.text.x = element_text(size=15,face="bold"),
                   axis.text.y = element_text(size=15,face="bold"),
                   axis.title.x = element_blank(),
                   axis.title.y = element_blank(),
                   # legend.title = element_blank(),
                   # legend.text = element_blank(),
                   # legend.position = "none",
                   panel.grid.major = element_blank(),
                   panel.border = element_blank(),
                   panel.background = element_blank(),
                   axis.ticks = element_blank())+
    ggplot2::coord_fixed()

  ####################################################################
  ####################################################################
  ####################################################################

  # TABLE.AclustsCCA <- summary_AclustsCCA(obj=AclustsCCA.result,annot=annot,n.top=n.top)
  # TABLE.AclustsCCA <- TABLE.AclustsCCA[TABLE.AclustsCCA$Significant=="Yes",]
  # TABLE.AclustsCCA <- data.table(TABLE.AclustsCCA)
  # unique.CpGs <- TABLE.AclustsCCA[,as.character(unique(do.call("c",CpGs)))]
  # annot.subset <- annot[IlmnID %in% unique.CpGs,]
  # annot.subset <- annot.subset[order(CHR,Coordinate_37),]
  # annot.subset[,IlmnID := as.character(IlmnID)]
  # annot.subset[,Exposures := sapply(unique.CpGs, function(x) {
  #   unlist(TABLE.AclustsCCA[sapply(CpGs,function(y) {sum(y%in% x)!=0}),Exposures])}
  #   )]
  # annot.subset[,CHR := as.factor(CHR)]
  # annot.subset[,UCSC_RefGene_Name:=sapply(UCSC_RefGene_Name,function(x) paste0(unique(unlist(strsplit(x,";"))),collapse=";"))]
  # annot.subset[,UCSC_RefGene_Name:=as.character(UCSC_RefGene_Name)]
  #
  # PM25.Names <- names(AclustsCCA.result$ALPHA.observed[[1]])
  #
  # DAT.FOR.PLOT <- t(sapply(unique.CpGs,function(y){ # for AclustsCCA
  #   PM25.Names %in% unlist(strsplit(unlist(annot.subset[IlmnID%in%y,Exposures]),","))
  # }))
  # DAT.FOR.PLOT <- 1*DAT.FOR.PLOT
  # DAT.FOR.PLOT[DAT.FOR.PLOT==1] <- "Selected"
  # DAT.FOR.PLOT[DAT.FOR.PLOT==0] <- "Not Selected"
  # colnames(DAT.FOR.PLOT) <- PM25.Names
  # rownames(DAT.FOR.PLOT) <- unique.CpGs
  #
  # title <- "Heatmap of AclustsCCA"
  # row_fontsize <- 3
  # # (1) Re-format Data
  # library("ComplexHeatmap")
  # library("circlize")
  # col_fun = colorRamp2(c(0, 5, 10), c("blue", "white", "red"))
  # # (2) Create pHeatmap
  # PLOT.HEATMAP <- Heatmap(DAT.FOR.PLOT,
  #                         column_title = title, # title
  #                         column_title_gp = gpar(fontsize = 15, fontface = "bold"), # title size
  #                         col = structure(c("Red","White"),
  #                                         names = c("Selected","Not Selected")), # block color
  #                         rect_gp = gpar(col = "bisque4", lwd = 1), # border color
  #                         name = " ", # legend title
  #                         cluster_columns=FALSE,
  #                         cluster_rows=TRUE,
  #                         column_names_rot = 0, # x-axis text no rotation
  #                         column_names_side = "top", # x-axis text on top
  #                         row_names_centered = TRUE,
  #                         column_names_centered = TRUE,
  #                         row_names_gp=gpar(fontsize=row_fontsize),
  #                         # row_title_gp = gpar(fill = c("green", "orange", "purple","yellow")),
  #                         row_split = annot.subset$UCSC_RefGene_Name,
  #                         row_title = NULL,
  #                         show_heatmap_legend = FALSE,
  #                         left_annotation = rowAnnotation(Gene=annot.subset$UCSC_RefGene_Name))
  # # col = list(Gene = c("HAPLN1"="red", "C6orf27"="green", "CCDC8"="blue", "COX4I2"="yellow"))))
  return(PLOT.HEATMAP)
}

