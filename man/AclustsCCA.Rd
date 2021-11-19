\name{AclustsCCA}
\alias{AclustsCCA}
\title{Implement Sparse Canonical Correlation Analysis on each cluster identified by A-clustering}
\usage{
AclustsCCA(X,Y,clusters.list,method="adaptive lasso",Ypenalty="SGL",init.method="SVD",max.iter=100,conv=10^-2)}

\arguments{
  \item{X}{an n X p matrix of p exposure varaible for n subjects.}
  \item{Y}{an q X n matrix of DNA methylation values (either beta-values or M-values) where \eqn{q} is number of DNA methylation sites and \eqn{n} is number of subjects.}
  \item{method}{either "lasso", "adaptive lasso", "gglasso", or "SGL". If Ypenalty is NULL, then this penalty is applied to both exposure and outcome.}
  \item{Ypenalty}{either "lasso", "adaptive lasso", "gglasso", or "SGL". This is a penalty for outcome (when regressing outcome Y onto exposure X to enable variable selection on outcome). If Ypenalty If Ypenalty is NULL, then penalty defined by method is applied to both exposure and outcome.}
  \item{init.method}{either "OLS" or "SVD". Initialization method}
  \item{max.iter}{maximum iteration of Sparse CCA.}
  \item{conv}{Convergence criteria of Sparse CCA}
}

\value{
Implements Sparse Canonical Correlation Analysis (sCCA) on each clusters identifed by A-clustering.
}
\description{
Add together two numbers
}
\examples{
sCCA.observed.FIT <- sCCA.observed(X=X,Y=Y,clusters.list=clusters.list)
}
