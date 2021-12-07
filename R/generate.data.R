
#' @title
#' Generate Synthetic Data for AclustsCCA R Package
#'
#' @description
#' Generates synthetic data with five exposures (p=5) and DNA methylation data with 14 CpG sites (q=14) in chromosome 7.
#' The data is generated from multivariate normal distribution where first three exposures (exposure 1-3) are associated with a single DMR region with three CpG sites (“cg03120555”,“cg04952627”,“cg10251229”).
#'
#'
#' ### INPUT
#' @param n : A number of subjects
#'
#' @return
#' The function returns a list of objects:
#'   - DATA.X: A \eqn{n} by \eqn{p=12} matrix where \eqn{n} is a number of subjects and \eqn{p=12} is a number of exposures.
#'   - DATA.Y: A \eqn{n} by \eqn{q=14} matrix where \eqn{n} is a number of subjects and \eqn{q=14} is a number of CpG sites.
#'   - TRUE.ALPHA: A vector of length \eqn{p=12} which represents true loading vector \eqn{alpha}
#'   - TRUE.BETA: A vector of length \eqn{q=14} which represents true loading vector \eqn{beta}
#'   - SIGMA.XX: A \eqn{p} by \eqn{p} exposure correlation matrix used to generated data
#'   - SIGMA.YY: A \eqn{q} by \eqn{q} outcome correlation matrix used to generated data
#'   - SIGMA.XY: A \eqn{p} by \eqn{q} correlation matrix of exposure and outcome used to generated data

#' @export
#'
#' @examples
#' data.list <- generate.data(n=500)
#' str(data.list)
#'
generate.data <- function(n=500){

  # Define correlation structure of outcome
  SIGMA.YY <- list(AR.matrix(n=6,rho=0.9),
                   AR.matrix(n=2,rho=0.9),
                   AR.matrix(n=3,rho=0.9),
                   AR.matrix(n=1,rho=1),
                   AR.matrix(n=1,rho=1),
                   AR.matrix(n=1,rho=1))
  SIGMA.YY <- as.matrix(Matrix::bdiag(SIGMA.YY))
  SIGMA.YY[SIGMA.YY==0] <- 0.1
  colnames(SIGMA.YY) <- rownames(SIGMA.YY) <- c("cg00725145","cg16080333","cg23568068","cg05107246","cg15535638","cg20044143",
                                                "cg16625683","cg04650322",
                                                "cg03120555","cg04952627","cg10251229",
                                                "cg03888078",
                                                "cg27406091",
                                                "cg27406091")
  q <- ncol(SIGMA.YY)

  # Define correlation structure of exposure
  SIGMA.XX <- list(AR.matrix(n=4,rho=0.9),
                   AR.matrix(n=3,rho=0.9),
                   AR.matrix(n=2,rho=0.9),
                   AR.matrix(n=1,rho=1),
                   AR.matrix(n=1,rho=1),
                   AR.matrix(n=1,rho=1))
  SIGMA.XX <- as.matrix(Matrix::bdiag(SIGMA.XX))
  SIGMA.XX[SIGMA.XX==0] <- 0.1
  p <- ncol(SIGMA.XX)
  colnames(SIGMA.XX) <- rownames(SIGMA.XX) <- paste("exposure",1:p,sep="")

  # Define true elements
  TRUE.ALPHA <- (rownames(SIGMA.XX) %in% c("exposure1","exposure2","exposure3"))
  TRUE.BETA <- (colnames(SIGMA.YY) %in% c("cg03120555","cg04952627","cg10251229"))

  # Normalize loading
  TRUE.ALPHA <- TRUE.ALPHA / as.numeric(sqrt(t(TRUE.ALPHA) %*% SIGMA.XX %*% TRUE.ALPHA))
  TRUE.BETA<- TRUE.BETA / as.numeric(sqrt(t(TRUE.BETA) %*% SIGMA.YY %*% TRUE.BETA))

  names(TRUE.ALPHA) <- rownames(SIGMA.XX)
  names(TRUE.BETA) <- rownames(SIGMA.YY)

  # DEFINE COVARIANCE MATRIX OF X AND Y (BASED ON CHEN'S PAPER)
  SIGMA.XY <- 0.9 * SIGMA.XX %*% TRUE.ALPHA %*% t(TRUE.BETA) %*% SIGMA.YY

  # CREATE COVARIANCE MATRIX FOR X AND Y
  SIGMA_BIG<-rbind(cbind(SIGMA.XX,SIGMA.XY),cbind(t(SIGMA.XY),SIGMA.YY))

  # Generate exposure data
  set.seed(20210322)
  DATA_clean<-mvtnorm::rmvnorm(floor(n),mean=rep(0,p+q),sigma=SIGMA_BIG)
  DATA.X<-DATA_clean[,1:p]
  DATA.Y<-DATA_clean[,(p+1):(p+q)]

  colnames(DATA.X) <- colnames(SIGMA.XX)
  colnames(DATA.Y) <- colnames(SIGMA.YY)

  rownames(DATA.X) <- rownames(DATA.Y) <- paste("Subject",1:n,sep="")
  return(list("DATA.X"=DATA.X,"DATA.Y"=DATA.Y,"TRUE.ALPHA"=TRUE.ALPHA,"TRUE.BETA"=TRUE.BETA,"SIGMA.XX"=SIGMA.XX,"SIGMA.YY"=SIGMA.YY,"SIGMA.XY"=SIGMA.XY))
}
