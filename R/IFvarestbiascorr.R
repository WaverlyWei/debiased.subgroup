#' Cross-validation metric
#'
#'@param Ycount Y
#'@param alphaEst estimated values
#'@param n sample size
#'@param splitSize size of each split
#'@return
#'mean squared error
#'@export

IFvarestbiascorr <- function(Ycount,
                             alphaEst,
                             n = NULL,
                             splitSize = NULL){

    if(is.null(splitSize)){

        stop("Specify the split size")
      }


  n <- length(Ycount[1,])

  Brep <- length(alphaEst)

  n2 <- n - splitSize

  cov.vec <- t(sweep(Ycount, 2, apply(Ycount, 2, mean))) *
    (alphaEst - mean(alphaEst))/Brep

  var.bias <- sum((alphaEst - mean(alphaEst))^2 )/Brep^2 * n  * n2/(n-n2)

  var.est <- (sum(cov.vec^2))

  r.corr <- (n/splitSize)^2 * ((n-1)/n)^2

  var.est0 <- var.est * r.corr

  return(sqrt(abs(var.est0- var.bias)))
}
