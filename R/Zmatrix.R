#' Generate the nodewise Lasso matrix used in desparsified Lasso

#'@param x nodewise confounder matrix
#'@param G indices of subgroups
#'@return
#'\item{Z}{nodewise Lasso matrix}
#' @export

Zmatrix <- function(x,
                    G = NULL){

  if(is.null(G)){

    stop("Number of subgroups is missing.")

  }


  n <- length(x[,1]) #number of subjects

  k <- length(G) # number of subgroups

  Z <- matrix(0,n,k)

  for(i in 1:k){

    index <- G[i]

    fit.nodelasso <-  cv.glmnet(x = x[,-index], y = x[,index])

    lambda <- fit.nodelasso$lambda.1se # obtain lasso tuning

    lambda <- lambda * 1.1 # scale up the lasso tuning parameter

    eta <- coef(fit.nodelasso, s = lambda) # lasso coefficients

    predi <- eta[1] + x[,-index] %*% eta[-1]

    # calculate the residual for the i-th column
    Z[,i] <-  x[,index]-predi
  }
  return(Z)
}
