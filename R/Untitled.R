#'variance calculation
#' @export


IFvarestbiascorr <- function(Y.count, alpha.est, n, split.size){

  n <- length(Y.count[1,]) # only filtered beta's

  Brep <- length(alpha.est)

  n2 <- n-split.size

  cov.vec <- t(sweep(Y.count, 2, apply(Y.count, 2, mean))) *
    (alpha.est - mean(alpha.est))/Brep

  var.bias <- sum((alpha.est - mean(alpha.est))^2 )/Brep^2 * n  * n2/(n-n2)

  var.est <- (sum(cov.vec^2))

  r.corr <- (n/split.size)^2 * ((n-1)/n)^2

  var.est0 <- var.est * r.corr

  return(sqrt(abs(var.est0- var.bias)))
}
