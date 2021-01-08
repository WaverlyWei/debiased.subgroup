#' Find the optimal tuning parameter for bootstrap-calibrated
#' desparsified Lasso
#'
#' @param y: response
#' @param x: design matrix
#' @param r: candidate tuning parameters
#' @param G: indices of subgroups
#' @param B: bootstrap repetitions
#' @param ll: number of folds in cross-validation
#' @return
#' op: optimal tuning parameter's index
#' @export
#'
cvDesparse <- function(y, x, r, G, B, ll)
{

  p <- length(x[1,]) # dimension of the design matrix

  n <- length(y) # number of subjects

  k <- length(G) # number of subgroups

  cc <- length(r) # number of tuning parameters

  n1 <- n/ll # divide indices

  mseMat <- matrix(0,cc,k) # record MSE of each tuning parameter

  # iterate through subsamples
  for(j in 1:ll) {

    # training set
    yb <- y[-((n1*(j-1)):(n1*j))]

    xb <- x[-((n1*(j-1)):(n1*j)),]

    # test set
    yt <- y[(n1*(j-1)):(n1*j)]

    xt <- x[(n1*(j-1)):(n1*j),]

    nn <- length(yb)

    fit.lasso <- cv.glmnet(x = xb, y = yb)

    # extract lambda
    lambda <- fit.lasso$lambda.1se

    # scale lambda
    lambda <- lambda *1.1

    # lasso estimates
    beta.lasso <-  coef(fit.lasso, s = lambda)

    gamma.lasso <-  beta.lasso[G+1]

    predb <-  beta.lasso[1] + xb %*% beta.lasso[-1]

    # residual for the training set
    residualb <-  yb - predb

    epsilionb <- residualb - mean(residualb)

    # initialize desparsified lasso estimates
    gamma.Dlasso <-  0

    Zb <- Zmatrix(xb,G)

    for(i in 1:k){

      index <- G[i]

      gamma.Dlasso[i] <- gamma.lasso[i] + sum(Zb[,i]*residualb)/sum(Zb[,i]*xb[,index])
    }

    TB <- matrix(0, B, cc)

    c <- matrix(0,cc,k)

    for(i in 1:cc){

      r0 <- r[i]

      c[i,] <- (1-n^(r0-0.5))*(max(gamma.lasso)-gamma.lasso)
    }

    for(i in 1:B){

      Bepsilion <- rnorm(nn)*epsilionb

      By <- predb + Bepsilion

      Bfit.lasso <- cv.glmnet(x = xb, y = By)

      Blambda <- Bfit.lasso$lambda.1se

      Blambda <- Blambda * 1.1

      Bbeta.lasso <- coef(Bfit.lasso, s = Blambda)

      Bgamma.lasso <- Bbeta.lasso[G+1]

      Bpred <- Bbeta.lasso[1] + xb%*%Bbeta.lasso[-1]

      Bresidual <- By - Bpred

      # initialize desparsified lasso
      Bgamma.Dlasso <- 0

      for(l in 1:k){

        Bindex <- G[l]

        Bgamma.Dlasso[l] <- Bgamma.lasso[l] + sum(Zb[,l]*Bresidual)/sum(Zb[,l]*xb[,Bindex])
      }

      for(l in 1:cc){

        TB[i,l] <- max(Bgamma.Dlasso+c[l,])-max(gamma.lasso)
      }
    }

    brEstimate <- 0

    for(l in 1:cc){ #bias-reduced estimate
      brEstimate[l]=max(gamma.lasso)-mean(TB[,l])}

    fit.lassot = cv.glmnet(x = xt, y = yt)

    lambdat <- fit.lassot$lambda.1se

    lambdat <- lambdat * 1.1

    beta.lassot <- coef(fit.lassot, s = lambdat)

    gamma.lassot <- beta.lassot[G+1]

    predt <- beta.lassot[1] + xt %*% beta.lassot[-1]

    residualt <- yt - predt

    sd <- 0

    gamma.Dlassot <- 0

    Zt <- Zmatrix(xt,G)

    for(i in 1:k){

      index <- G[i]

      gamma.Dlassot[i] <- gamma.lassot[i] + sum(Zt[,i]*residualt)/sum(Zt[,i]*xt[,index])
    }

    for(i in 1:k){
      index <- G[i]

      sd[i] <- sqrt(var(residualt*Zt[,i])/n1)/abs(sum(Zt[,i]*xt[,index])/n1)
    }
    for(i in 1:cc){
      mseMat[i,] <- mseMat[i,]+(brEstimate[i]-gamma.Dlassot)**2-sd**2}
  }

  mse <- 0

  for(i in 1:cc){
    mse[i] <- min(mseMat[i,])
  }

  op <- which.min(mse)

  return(op)
}

