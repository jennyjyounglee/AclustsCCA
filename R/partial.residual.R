#' Title
#'
#' @title
#' partial.residual
#'
#' @description
#' Regress `data` onto `Z`
#'
#' @param data    : A data matrix with \eqn{n} rows.
#' @param Z       : A \eqn{n} by \eqn{r} confounder data matrix, where \eqn{n} is sample size and \eqn{r} is number of potential confounders
#' @param nthread : A number of threads to parallelize regression.
#'
#' @export
#' @return
#' Returns a data matrix with same dimension as `data`
#'
#'
#'
#'
partial.residual <- function(data,Z,nthread){
  p <- ncol(data)

  cl <- parallel::makeCluster(nthread)
  doSNOW::registerDoSNOW(cl)
  out <- foreach(j=1:p, .combine=rbind) %dopar% {
    exposure <- data[,j]
    dat <- data.frame(exposure,Z)
    return(lm(exposure ~. ,data=dat)$residuals)
  }
  out <- t(out)
  colnames(out) <- colnames(data)
  rownames(out) <- rownames(data)

  parallel::stopCluster(cl)
  return(out)
}
