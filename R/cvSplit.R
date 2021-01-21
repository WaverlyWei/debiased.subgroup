#' Choose the optimal tuning for bootstrap-calibrated R-Split
#'
#' @param y: response
#' @param x: design matrix
#' @param r: candidate tuning parameter
#' @param G: subgroup indicator
#' @param B: bootstrap iterations
#' @param BB: split number
#' @param ratio: ratio of data splitting
#' @param fold: number of folds in cross-validation
#' @return
#' op: index of the optimal tuning
#' @export

cvSplit <- function(y, x,
                    r = NULL,
                    G = NULL,
                    B = NULL,
                    BB = NULL,
                    ratio = NULL,
                    fold = 2){

  if(is.null(r)){

    stop("Tuning parameter is missing.")

  }else if(is.null(G)){

    stop("Number of subgroups is missing.")

  }else if(is.null(B)){

    stop("Specify bootstrap iterations.")

  }else if(is.null(BB)){

    stop("Specify refitting iterations.")

  }

  p <- length(x[1,])

  n <- length(y)

  k <- length(G)

  cc <- length(r)

  penalty <- rep(0,p)

  penalty[-G] <- 1

  op <- cc

  n1 <- n/fold   # sample split

  mseMat <- matrix(0,cc,k)

  for(j in 1:fold){    #iterate through subsample

    ytrain <- y[-((n1*(j-1)+1):(n1*j))]

    xtrain <- x[-((n1*(j-1)+1):(n1*j)),]

    ytest <- y[(n1*(j-1)+1):(n1*j)]

    xtest <- x[(n1*(j-1)+1):(n1*j),]

    ntrain <- length(ytrain) # nn -> ntrain

    nsub <- ntrain*ratio

    beta.lassob <- 0
    Deltab <- 0

    est <- matrix(NA,nrow = BB, ncol = k )

    for(b in 1:BB){

      index <- sample(1:ntrain,nsub)

      # selection set
      ytrain.s <- ytrain[index]

      xtrain.s <- xtrain[index,]

      # refitting set
      ytrain.r <- ytrain[-index]

      xtrain.r <- xtrain[-index,]

      # adaptive lasso
      fit.ridge <- cv.glmnet(x = xtrain.s, y=ytrain.s,
                             penalty.factor = penalty,
                             alpha=0)

      beta.ridge <- coef(fit.ridge, s= "lambda.min")[-1]

      fit.lasso <- cv.glmnet(x = xtrain.s,
                             y = ytrain.s,
                             penalty.factor = penalty/abs(beta.ridge))

      # set bounds for model size
      modIdx <- (fit.lasso$nzero>(1+k))&(fit.lasso$nzero<(k+5+k))

      ss1 <- fit.lasso$lambda[modIdx]

      mcvError <- fit.lasso$cvm[modIdx] # mean cv error

      #the optimal tuning
      ss <- ss1[which.min(mcvError)]

      # selection set
      gamma.lasso <- coef(fit.lasso, s = ss)[2:(p+1)]

      set1 <- setdiff(which(gamma.lasso!=0), G)

      set1 <- sort(set1)

      #refit estimate
      refit_est <- lm(ytrain.r~xtrain.r[,G]+xtrain.r[,set1])$coef[2:(k+1)]

      est[b,] <- refit_est

      ZZtrain <- cbind(xtrain.r[,G],xtrain.r[,set1])

      Delta0b <-  matrix(0,k,p)

      Delta0b[,union(G,set1)] <- solve(t(ZZtrain) %*% ZZtrain/(ntrain-nsub))[G,]

      Deltab <- Deltab + Delta0b
    }

    # filter results
    filtered_index <- apply( est, 2, function(x) order(x)[ceiling(( (1-0.95)*nrow(est))):(floor((1-0.05)*nrow(est)))])

    filtered_beta <- NULL

    for(i in 1:k){

      filtered_beta[i]<- mean(est[filtered_index[,i],i])
    }

    beta.lassob <- filtered_beta

    Deltab <- Deltab/BB

    fit.lassob <- cv.glmnet(x = xtrain, y = ytrain)

    gamma.lassob <- coef(fit.lassob, s = "lambda.min")

    predb <- gamma.lassob[1] + xtrain%*%gamma.lassob[-1]

    residualb <- ytrain - predb

    epsilionb <- residualb-mean(residualb)

    TBB <- matrix(0, B, cc)

    correction <- matrix(0,cc,k) # correction matrix

    for(i in 1:cc){

      r0 <- r[i]

      correction[i,] <- (1-n^(r0-0.5))*(max(beta.lassob)-beta.lassob)
    }

    for(i in 1:B){

      Bepsilionb <- as.matrix(rnorm(ntrain)*epsilionb,ntrain,1)

      Bbeta.lassob <- beta.lassob + Deltab %*% t(cbind(xtrain[,G],xtrain[,-G])) %*% Bepsilionb/ntrain

      for(l in 1:cc)
        TBB[i,l] <- max(Bbeta.lassob+correction[l,])-max(beta.lassob)
    }

    brestimate <- 0

    for(l in 1:cc) #bias-reduced estimate
      brestimate[l] <- max(beta.lassob)-mean(TBB[,l])

    beta.lassott <- matrix(0,BB,k)

    nsub2 <- (n-ntrain)*ratio # rest of the subsamples

    Y.count <- matrix(data = NA, nrow = BB, ncol = n-ntrain)

    Delta0b <- 0

    est <- matrix(NA,nrow = BB, ncol = k )

    for(b in 1:BB){

      index <- sample(1:(n-ntrain),nsub2)

      ytest.s <- ytest[index]

      xtest.s <- xtest[index,]

      ytest.r <- ytest[-index]

      xtest.r <- xtest[-index,]

      # adaptive lasso
      fit.ridge <- cv.glmnet(x = xtest.s, y=ytest.s, penalty.factor = penalty, alpha=0)

      beta.ridge <- coef(fit.ridge, s= "lambda.min")[-1]

      fit.lasso <- cv.glmnet(x = xtest.s,
                             y = ytest.s,
                             penalty.factor = penalty/abs(beta.ridge))

      # set bounds for model size
      modIdx <- (fit.lasso$nzero>(1+k))&(fit.lasso$nzero<(k+5+k))

      ss1 <- fit.lasso$lambda[modIdx]

      mcvError <- fit.lasso$cvm[modIdx]

      # the optimal tuning
      ss <- ss1[which.min(mcvError)]

      #selection set
      gamma.lasso <- coef(fit.lasso, s = ss)[2:(p+1)]

      set1 <- setdiff(which(gamma.lasso!=0), G)

      set1 <- sort(set1)

      #refit estimate
      refit_est <- lm(ytest.r~xtest.r[,G]+xtest.r[,set1])$coef[2:(k+1)]

      est[b,] <- refit_est

      Y.count[b,] <- tabulate(index,n-ntrain)
    }

    # filter results
    filtered_index <- apply( est, 2, function(x) order(x)[ceiling(( (1-0.95)*nrow(est))):(floor((1-0.05)*nrow(est)))])

    filtered_beta <- NULL

    beta.lassott <- matrix(NA, nrow = nrow(filtered_index), ncol = k)

    for(i in 1:k){

      beta.lassott[,i] <- est[filtered_index[,i],i]

      filtered_beta[i]<- mean(est[filtered_index[,i],i])
    }

    beta.lassot <- filtered_beta

    good_iter <- order(table(c(filtered_index)))[-c(1: (BB-nrow(filtered_index)))]

    Y.count <- Y.count[good_iter,]

    sd <- 0

    for(i in 1:k){

      sd[i] <- IFvarestbiascorr(Y.count, beta.lassott[,i], n-ntrain, splitSize = nsub2)
    }

    for(i in 1:cc)
      mseMat[i,] <- mseMat[i,]+(brestimate[i]-beta.lassot)**2-sd**2
  }

  mse <- 0

  for(i in 1:cc)
    mse[i] <- min(mseMat[i,])

  op <- which.min(mse)

  return(op)
}

