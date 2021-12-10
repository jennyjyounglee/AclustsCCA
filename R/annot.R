#' @title Sample Annotation Data for Analysis
#'
#' @description
#' A \code{data.table} that contains annotation for the probes on chromosome 7 in the Illumina 450K chip,
#' sorted by chromosome and chromosomal coordinates of the CpG.
#' created using the R package IlluminaHumanMethylation450kanno.ilmn12.hg19 from Bioconductor.
#'
#'
#' @format A \code{data.table} with 6 column, which are:
#' \describe{
#' \item{IlmnID}{a character vector of IlmnID}
#' \item{CHR}{a numeric vector of chromosome containing the CpG}
#' \item{Coordinate_37}{a numeric vector of chromosomal coordinates of the CpG}
#' \item{UCSC_RefGene_Name}{a character vector of Target gene name(s)}
#' \item{Islands_Name}{a character vector of chromosomal coordinates of the CpG Island}
#' \item{Relation_to_Island}{a character vector of the location of the CpG relative to the CpG island}
#' }
#'
"annot"


