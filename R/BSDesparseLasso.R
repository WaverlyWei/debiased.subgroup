#' Bootstrap-calibrated Desparsified Lasso
#'
#'This method first constructs the debiased estimator of \eqn{\beta} via the
#'desparsified Lasso procedure. Then it calculates the calibration term
#' \eqn{\hat{b}_{max} =(1-n^{r-0.5})(\hat{\beta}_{max}-\hat{\beta}_{j,lasso})}. Through B bootstrap iterations,
#' it recalibrates the bootstrap statistic \eqn{T_b}. The bias-reduced estimate
#' is computed as: \eqn{\hat{b}_{max}-\frac{1}{B}\sum_{b=1}^BT_b}.
#'
#'
#' @param y response
#' @param x design matrix
#' @param r tuning parameter
#' @param G subgroup indicator
#' @param B bootstrap iterations
#' @param alpha level of CI
#' @return
#' \item{LowerBound}{lower confidence bound}
#' \item{UpperBound}{upper confidence bound}
#' \item{betaMax}{bias-reduced maximum beta estimate}
#' \item{betaEst}{debiased beta estimate for each subgroup}
#' \item{op}{optimal tuning}
#' @export
BSDesparseLasso <- function(y, x,
                            r = NULL,
                            G = NULL,
                            B = NULL,
                            alpha = 0.95,
                            fold = 3){

  if(is.null(r)){

    stop("Tuning parameter is missing.")

  }else if(is.null(G)){

    stop("Number of subgroups is missing.")

  }else if(is.null(B)){

    stop("Specify bootstrap iterations.")
  }

  p <- length(x[1,])

  n <- length(y)

  k <- length(G) # number of subgroups

  cc <- length(r) # number of candidate tuning parameters

  y <- y - mean(y)


  for(i in 1:p)
    x[,i] <- x[,i]-mean(x[,i])

  # lasso estimates
  fit.lasso <- cv.glmnet(x = x, y = y)

  lambda <- fit.lasso$lambda.1se

  lambda <- lambda *1.1

  gamma.lasso <- coef(fit.lasso, s = lambda)

  beta.lasso <- gamma.lasso[G+1]

  # calculate residual
  pred <- gamma.lasso[1] + x %*% gamma.lasso[-1]

  residual <- y - pred

  epsilion <- residual-mean(residual)

  #Desparsified  lasso estimate
  beta.Dlasso <- 0

  Z <- Zmatrix(x,G)

  for(i in 1:k){

    index <- G[i]

    beta.Dlasso[i] <- beta.lasso[i] + sum(Z[,i]*residual)/sum(Z[,i]*x[,index])
  }

  #calculate the correction term
  TB <- matrix(0, B, cc+1)

  correction <- matrix(0, cc+1, k)

  for(i in 1:cc) {
    r0 <- r[i]
    correction[i,] = (1-n^(r0-0.5))*(max(beta.lasso)-beta.lasso)
  }

  #the simulataneous one
  correction[cc+1,] <- (max(beta.lasso)-beta.lasso)

  TB_op <- matrix(0, B, cc)

  c_op <- matrix(0, cc, k)

  for(i in 1:cc) {
    r0 <- r[i]

    rp <- r0/sqrt(k/2)

    c_op[i,] <- (1-n^(rp-0.5))*(max(beta.lasso)-beta.lasso)
  }


  for(i in 1:B){

    #generate bootstrap sample
    Bepsilion <- rnorm(n)*epsilion

    By <- pred + Bepsilion

    #calculate bootstrap desparsified lasso
    Bfit.lasso <- cv.glmnet(x = x, y = By)

    Blambda <- Bfit.lasso$lambda.1se

    Blambda <- Blambda * 1.1

    Bgamma.lasso <- coef(Bfit.lasso, s = Blambda)

    Bbeta.lasso <- Bgamma.lasso[G+1]

    Bpred <- Bgamma.lasso[1] + x %*% Bgamma.lasso[-1]

    Bresidual <- By - Bpred

    Bbeta.Dlasso <- 0

    for(j in 1:k){

      Bindex <- G[j]
      Bbeta.Dlasso[j] <- Bbeta.lasso[j] + sum(Z[,j]*Bresidual)/sum(Z[,j]*x[,Bindex])
    }

    #correct the maximum quantity
    for(j in 1:(cc+1)){
      TB[i,j] <- max(Bbeta.Dlasso+correction[j,])-max(beta.lasso)
    }

    for(j in 1:(cc)){

      TB_op[i,j] <- max(Bbeta.Dlasso+c_op[j,])-max(beta.lasso)

    }
  }

  #choose the optimal tuning parameter

  op <- cvDesparse(y, x, r, G, B, fold) # index of the correction term

  # collect results
  result <- list()

  for(j in 1:(cc+1)) {

    result[j] = list(c(BSciCoverfun(beta.Dlasso, TB[,j], G, alpha),
                       betaEst = list(beta.Dlasso),
                       op = r[op]))
  }

  if(is.integer(op)){
    result[j+1] = list(c(BSciCoverfun(beta.Dlasso, TB_op[,op], G, alpha),
                         betaEst = list(beta.Dlasso),
                         op = r[op]))
  }else{
    result[j+1] = list(c(BSciCoverfun(beta.Dlasso, TB[,cc], G, alpha),
                         betaEst = list(beta.Dlasso),
                         op = r[op]))
  }
  return(result[[12]])
}
