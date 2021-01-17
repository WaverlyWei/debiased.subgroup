#' Find the optimal tuning parameter for bootstrap-calibrated
#' desparsified Lasso
#'
#' @param y: response
#' @param x: design matrix
#' @param r: candidate tuning parameters
#' @param G: indices of subgroups
#' @param B: bootstrap repetitions
#' @param fold: number of folds in cross-validation
#' @return
#' op: optimal tuning parameter's index
#' @export
#'
cvDesparse <- function(y, x,
                       r = NULL,
                       G = NULL,
                       B = NULL,
                       fold =NULL){

  p <- length(x[1,]) # dimension of the design matrix

  n <- length(y) # number of subjects

  k <- length(G) # number of subgroups

  cc <- length(r) # number of tuning parameters

  n1 <- n/fold # divide indices

  mseMat <- matrix(0,cc,k) # record MSE of each tuning parameter

  # iterate through subsamples
  for(j in 1:fold) {

    # training set
    ytrain <- y[-((n1*(j-1)):(n1*j))]

    xtrain <- x[-((n1*(j-1)):(n1*j)),]

    # test set
    ytest <- y[(n1*(j-1)):(n1*j)]

    xtest <- x[(n1*(j-1)):(n1*j),]

    ntrain <- length(ytrain)

    fit.lasso <- cv.glmnet(x = xtrain, y = ytrain)

    # extract lambda
    lambda <- fit.lasso$lambda.1se

    # scale lambda
    lambda <- lambda *1.1

    # lasso estimates
    gamma.lasso <-  coef(fit.lasso, s = lambda)

    beta.lasso <-  gamma.lasso[G+1]

    predb <-  gamma.lasso[1] + xtrain %*% gamma.lasso[-1]

    # residual for the training set
    residualb <-  ytrain - predb

    epsilionb <- residualb - mean(residualb)

    # initialize desparsified lasso estimates
    beta.Dlasso <-  0

    Zb <- Zmatrix(xtrain,G)

    for(i in 1:k){

      index <- G[i]

      beta.Dlasso[i] <- beta.lasso[i] + sum(Zb[,i]*residualb)/sum(Zb[,i]*xtrain[,index])
    }

    TB <- matrix(0, B, cc)

    correction <- matrix(0,cc,k)

    for(i in 1:cc){

      r0 <- r[i]

      correction[i,] <- (1-n^(r0-0.5))*(max(beta.lasso)-beta.lasso)
    }

    for(i in 1:B){

      Bepsilion <- rnorm(ntrain)*epsilionb

      By <- predb + Bepsilion

      Bfit.lasso <- cv.glmnet(x = xtrain, y = By)

      Blambda <- Bfit.lasso$lambda.1se

      Blambda <- Blambda * 1.1

      Bgamma.lasso <- coef(Bfit.lasso, s = Blambda)

      Bbeta.lasso <- Bgamma.lasso[G+1]

      Bpred <- Bgamma.lasso[1] + xtrain%*%Bgamma.lasso[-1]

      Bresidual <- By - Bpred

      # initialize desparsified lasso
      Bbeta.Dlasso <- 0

      for(l in 1:k){

        Bindex <- G[l]

        Bbeta.Dlasso[l] <- Bbeta.lasso[l] + sum(Zb[,l]*Bresidual)/sum(Zb[,l]*xtrain[,Bindex])
      }

      for(l in 1:cc){

        TB[i,l] <- max(Bbeta.Dlasso+correction[l,])-max(beta.lasso)
      }
    }

    brEstimate <- 0

    for(l in 1:cc){ #bias-reduced estimate
      brEstimate[l] <- max(beta.lasso)-mean(TB[,l])
      }

    fit.lassot <- cv.glmnet(x = xtest, y = ytest)

    lambdat <- fit.lassot$lambda.1se

    lambdat <- lambdat * 1.1

    gamma.lassot <- coef(fit.lassot, s = lambdat)

    beta.lassot <- gamma.lassot[G+1]

    predt <- gamma.lassot[1] + xtest %*% gamma.lassot[-1]

    residualt <- ytest - predt

    sd <- 0

    beta.Dlassot <- 0

    Zt <- Zmatrix(xtest,G)

    for(i in 1:k){

      index <- G[i]

      beta.Dlassot[i] <- beta.lassot[i] + sum(Zt[,i]*residualt)/sum(Zt[,i]*xtest[,index])
    }

    for(i in 1:k){
      index <- G[i]

      sd[i] <- sqrt(var(residualt*Zt[,i])/n1)/abs(sum(Zt[,i]*xtest[,index])/n1)
    }
    for(i in 1:cc){
      mseMat[i,] <- mseMat[i,]+(brEstimate[i]-beta.Dlassot)**2-sd**2}
  }

  mse <- 0

  for(i in 1:cc){
    mse[i] <- min(mseMat[i,])
  }

  op <- which.min(mse)

  return(op)
}

