#' Coverage and distance for the naive methods
#'
#' This function computes the coverage, lower bound and length of the
#' naive methods
#'
#' @param Lowerbound: the lower bound returned from naive methods
#' @param beta0: true beta values
#' @param G: indices of subgroups
#' @return a list of coverage, lower bound and bound length
#'
BSciCoverfunNaive <- function(Lowerbound, beta0, G){

  gamma0 <- beta0[G] # true gammas

  result <- list(coverage = Lowerbound <= max(gamma0),
                 LowerBound = Lowerbound,
                 Length = abs(Lowerbound - max(gamma0)))
  return(result)
}
