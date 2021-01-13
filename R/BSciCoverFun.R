
#' Compute CI for bootstrap-calibrated methods
#' @param beta: estimated betas
#' @param TB: correction term
#' @param G: indices of subgroups
#' @param alpha0: confidence level
#' @return
#' coverage: boolean value
#'  LowerBound: Lower bound of the estimates
#'  UpperBound: Upper bound of the estimates
#'  betaMax: debiased maximum beta estimate
#'  @export
BSciCoverfun <- function(beta, TB, G, alpha0){

  LowerBound <- max(beta) - quantile(TB, alpha0)

  UpperBound <- max(beta) + quantile(TB, alpha0)

  result <- list(LowerBound = LowerBound,

                 betaMax = max(beta)-mean(TB))

  return(result)
}
