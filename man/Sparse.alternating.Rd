% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Sparse.alternating.R
\name{Sparse.alternating}
\alias{Sparse.alternating}
\title{Penalized least squares function used for SparseCCA}
\usage{
Sparse.alternating(Xreg, Yreg, method, groupidx = NULL)
}
\arguments{
\item{Xreg}{: A data matrix of \eqn{n} rows}

\item{Yreg}{: A vector of length \eqn{n}}

\item{Xmethod}{: penalty function for the exposure, i.e. penalty function when regressing Yreg onto Xreg. Possible values are:
\itemize{
\item "lasso": Lasso
\item "alasso": Adaptive Lasso
\item "gglasso": Group Lasso
\item "SGL": Sparse Group Lasso
\item "OLS": Ordinary Least Square
}}
}
\description{
Implement an penalized least squares needed to run sparse canonical correlation analysis (SparseCCA)
with various penalty functions. Modified Wilms, Ines, and Christophe Croux. "Robust sparse canonical correlation analysis." BMC systems biology 10.1 (2016): 1-13.
The original code is accessible https://sites.google.com/view/iwilms/software?authuser=0
\subsection{INPUT}{
}
}
\examples{
data.list <- generate.data(n=500)
DATA.X <- data.list$DATA.X
DATA.Y <- data.list$DATA.Y
Sparse.alternating.result <- Sparse.alternating(Xreg=DATA.X,Yreg=DATA.Y[,1],method="SGL")
str(Sparse.alternating.result)



}
