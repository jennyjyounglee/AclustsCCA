\name{annot}
\alias{annot}
\docType{data}
\title{
  Annotation for the probes on chromosome 7 in the Illumina 450K chip
}
\description{
  This data contains annotation for the probes on chromosome 7 in the Illumina 450K chip, created using the package IlluminaHumanMethylation450kanno.ilmn12.hg19 from Bioconductor. A code to create this data is provided below. This data is a data.table with annotation for probes, sorted by chromosome and chromosomal coordinates of the CpG.
}
\usage{data(annot)}
\format{
  A data frame with 28396 observations on the following 8 variables.
  \describe{
    \item{\code{IlmnID}}{a character vector of IlmnID}
    \item{\code{CHR}}{a numeric vector of chromosome containing the CpG}
    \item{\code{Coordinate_37}}{a numeric vector of chromosomal coordinates of the CpG}
    \item{\code{UCSC_RefGene_Name}}{a character vector of Target gene name(s)}
    \item{\code{Islands_Name}}{a character vector of chromosomal coordinates of the CpG Island}
    \item{\code{Relation_to_Island}}{a character vector of the location of the CpG relative to the CpG island
}
\examples{
  # This annotation file is created using the code below:
  # library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
  # annot <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  # annot <- data.table(data.frame(annot))
  # annot <- annot[!chr%in%c("chrY","chrX"),]
  # annot <- annot[,list("IlmnID"=Name,"Coordinate_37"=pos,"CHR"=as.numeric(gsub("chr", "", chr)),Islands_Name,Relation_to_Island,UCSC_RefGene_Name)]
  # annot <- annot[order(CHR,Coordinate_37),]
  data(annot)
  head(annot)
}


\examples{

}
