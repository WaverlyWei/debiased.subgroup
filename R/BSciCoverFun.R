
#' Compute CI for bootstrap-calibrated methods
#' @param beta: estimated betas
#' @param TB: correction term
#' @param beta0: true beta values
#' @param G: indices of subgroups
#' @param alpha0: confidence level
#' @return
#' coverage: boolean value
#'  LowerBound: Lower bound of the estimates
#'  Length: length of the lower bound
#'  Bias: bias of the calibrated beta estimates
#'  @export
BSciCoverfun <- function(beta, TB, beta0, G, alpha0){

  LowerBound <- max(beta)-quantile(TB, alpha0)

  gamma0 <- beta0[G]

  result <- list(coverage = LowerBound <= max(gamma0), #coverage
                 LowerBound = LowerBound,
                 Length = abs(LowerBound - max(gamma0)),
                 Bias = max(beta)-mean(TB)-max(gamma0))

  return(result)
}
