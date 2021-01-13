#' Bootstrap-calibrated Despasified Lasso
#'
#' @param y: response
#' @param x: design matrix
#'  @param r: tuning parameter
#'  @param G: subgroup indicator
#'  @param B: bootstrap iterations
#'  @param alpha0: level of CI
#'  @return
#'  coverage: boolean value
#'  LowerBound: lower bound
#'  Length: length of the lower bound
#'  betaEst: estimated beta values
#'  modelSize: selected model size as a reference
#'  op: index of the tuning parameter
#' @export
BSDesparseLasso <- function(y, x, r, G, B, alpha0){


  # Length: length of the lower bound
  # betaEst: estimated beta values
  # modelSize: selected model size as a reference
  # op: index of the tuning parameter

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

  beta.lasso <- coef(fit.lasso, s = lambda)

  gamma.lasso <- beta.lasso[G+1]

  # calculate residual
  pred <- beta.lasso[1] + x %*% beta.lasso[-1]

  residual <- y - pred

  epsilion <- residual-mean(residual)

  #Desparsified  lasso estimate
  gamma.Dlasso <- 0

  Z <- Zmatrix(x,G)

  for(i in 1:k){

    index <- G[i]

    gamma.Dlasso[i] <- gamma.lasso[i] + sum(Z[,i]*residual)/sum(Z[,i]*x[,index])
  }

  #calculate the correction term
  TB <- matrix(0, B, cc+1)

  correction <- matrix(0, cc+1, k)

  for(i in 1:cc) {
    r0 <- r[i]
    correction[i,] = (1-n^(r0-0.5))*(max(gamma.lasso)-gamma.lasso)
  }

  #the simulataneous one
  correction[cc+1,]= (max(gamma.lasso)-gamma.lasso)

  TB_op <- matrix(0, B, cc)

  c_op <- matrix(0, cc, k)

  for(i in 1:cc) {
    r0 <- r[i]

    r_op <- r0/sqrt(k/2)

    c_op[i,] <- (1-n^(r_op-0.5))*(max(gamma.lasso)-gamma.lasso)
  }

  modelSize <- NULL

  for(i in 1:B){

    #generate bootstrap sample
    Bepsilion <- rnorm(n)*epsilion

    By <- pred + Bepsilion

    #calculate bootstrap desparsified lasso
    Bfit.lasso <- cv.glmnet(x = x, y = By)

    Blambda <- Bfit.lasso$lambda.1se

    Blambda <- Blambda * 1.1

    Bbeta.lasso <- coef(Bfit.lasso, s = Blambda)

    # collect model size
    modelSize[i] <- length(which(Bbeta.lasso!=0))

    Bgamma.lasso <- Bbeta.lasso[G+1]

    Bpred <- Bbeta.lasso[1] + x %*% Bbeta.lasso[-1]

    Bresidual <- By - Bpred

    Bgamma.Dlasso <- 0

    for(j in 1:k){

      Bindex <- G[j]
      Bgamma.Dlasso[j] <- Bgamma.lasso[j] + sum(Z[,j]*Bresidual)/sum(Z[,j]*x[,Bindex])
    }

    #correct the maximum quantity
    for(j in 1:(cc+1)){
      TB[i,j] <- max(Bgamma.Dlasso+c[j,])-max(gamma.lasso)
    }

    for(j in 1:(cc)){

      TB_op[i,j] <- max(Bgamma.Dlasso+c_op[j,])-max(gamma.lasso)

    }
  }

  #choose the optimal tuning parameter
  ll <- 3

  op <- cvDesparse(y, x, r, G, B, ll) # index of the correction term

  # collect results
  result <- list()

  for(j in 1:(cc+1)) {

    result[j] = list(c(BSciCoverfun(gamma.dest, TB[,j], G, alpha0),
                       betaEst = list(gamma.dest),
                       modelSize = list(modelSize),
                       op = op))
  }

  if(is.integer(op)){
    result[j+1] = list(c(BSciCoverfun(gamma.dest, TB_op[,op], G, alpha0),
                         betaEst = list(gamma.dest),
                         modelSize = list(modelSize),
                         op = op))
  }else{
    result[j+1] = list(c(BSciCoverfun(gamma.dest, TB[,cc], G, alpha0),
                         betaEst = list(gamma.dest),
                         modelSize = list(modelSize),
                         op = op))
  }
  return(result)
}
