
#' Compute CI for bootstrap-calibrated methods
#'
#' @param beta estimated betas
#' @param TB recalibrated bootstrap statistics
#' @param G indices of subgroups
#' @param alpha confidence level
#' @return
#'  \item{LowerBound}{Lower confidence bound}
#'  \item{UpperBound}{Upper confidence bound}
#'  \item{betaMax}{debiased maximum beta estimate}
BSciCoverfun <- function(beta,
                         TB = NULL,
                         G = NULL,
                         alpha = 0.95){

  if(is.null(TB)){

    stop("Specify test statistic.")

  }else if(is.null(G)){
    stop("Specify number of subgroups")
  }

  LowerBound <- max(beta) - quantile(TB, alpha)

  UpperBound <- max(beta) + quantile(TB, alpha)

  result <- list(LowerBound = LowerBound,

                 UpperBound = UpperBound,

                 betaMax = max(beta)-mean(TB))

  return(result)
}
