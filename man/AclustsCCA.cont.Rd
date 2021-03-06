% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AclustsCCA.cont.R
\name{AclustsCCA.cont}
\alias{AclustsCCA.cont}
\title{Implement additional permutations of AclustsCCA}
\usage{
AclustsCCA.cont(obj, X, Y, maxnum = NULL, maxB = 10000)
}
\arguments{
\item{obj}{A result of \code{AclustsCCA} function.}

\item{X}{\eqn{n} by \eqn{p} exposure data matrix, where \eqn{n} is sample size and \eqn{p} is number of exposures.}

\item{Y}{\eqn{n} by \eqn{q} outcome data matrix, where \eqn{n} is sample size and \eqn{q} is number of outcomes.}

\item{maxnum}{A maximal total number of permutations across all the clusters.}

\item{maxB}{A maximal number of permutations for a single cluster.}
}
\value{
The function returns a list of 6 objects according to the following order:
\itemize{
\item clusters.list          : A list of clusters with CpG sites obtained using A-clustering, each item is a cluster that contains a set of probes. If A-clustering is not implemented inside AclustsCCA, return \code{NA}.
\item ALPHA.observed         : A list of estimated canonical vector of length \eqn{p} corresponding to the exposure data \eqn{X} for each cluster.
\item BETA.observed          : A list of estimated canonical vector of length \eqn{q} corresponding to the outcome data \eqn{Y} for each cluster.
\item cancors.observed       : A vector of estimated canonical correlation for each cluster.
\item permutation.result     : A \code{mmctest} object that contains permutation results.
\item settings               : A settings used for the analysis.
}
}
\description{
Implement an iterative penalized least squares approach to
sparse canonical correlation analysis (SparseCCA)
with various penalty functions.
}
\examples{



}
