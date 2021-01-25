#' Generate covariance matrix
#'
#' This function generates different types of covariance matrices
#'
#' @param p dimension of confounders
#' @param type type of matrix
#' @return \item{Sigma}{A covariance matrix}
#' @export

sigmaMatNew <- function(p, type = NULL){

  if(is.null(type)){

    stop("Specify matrix type.")

  }

  if(type == "ind"){

    Sigma <- matrix(data = 0, nrow = p, ncol = p)

    diag(Sigma) <- 1

  }

  if(type == "toeplitz3"){

    Sigma <- matrix(data = 0, nrow = p, ncol = p)

    p.vec <- 1:p

    Sigma <- 0.3^(toeplitz(p.vec)-1)

  }

  if(type == "toeplitz5"){

    Sigma <- matrix(data = 0, nrow = p, ncol = p)

    p.vec <- 1:p

    Sigma <- 0.5^(toeplitz(p.vec)-1)

  }

  if(type == "toeplitz9"){

    Sigma <- matrix(data = 0, nrow = p, ncol = p)

    p.vec <- 1:p

    Sigma <- 0.9^(toeplitz(p.vec)-1)

  }

  if(type == "block"){

    Sigma <- matrix(data = 0, nrow = p, ncol = p)

    Sigma[1:20, 1:20] <- 0.5
    Sigma[51:70, 51:70] <- 0.5
    Sigma[101:120, 101:120] <- 0.5
    Sigma[151:170, 151:170] <- 0.5
    Sigma[201:220, 201:220] <- 0.5
    Sigma[251:270, 251:270] <- 0.5

    for(i in 400:500){
      for(j in 400:500){

        Sigma[i,j] <- 0.5^{abs(i-j)}
      }
    }

    diag(Sigma) <- 1
  }

  return(Sigma)
}
