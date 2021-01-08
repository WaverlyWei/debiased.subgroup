#' Choose the optimal tuning for bootstrap-calibrated R-Split
#'
#' @param y: response
#' @param x: design matrix
#' @param r: candidate tuning parameter
#' @param G: subgroup indicator
#' @param B: bootstrap iterations
#' @param BB: split number
#' @param ratio: ratio of data splitting
#' @param ll: number of folds in cross-validation
#' @return
#' op: index of the optimal tuning
#' @export

cvSplit <- function(y, x, r, G, B, BB, ratio, ll)
{

  p <- length(x[1,])

  n <- length(y)

  k <- length(G)

  cc <- length(r)

  penalty <- rep(0,p)

  penalty[-G] <- 1

  op <- cc

  n1 <- n/ll   # sample split

  mseMat <- matrix(0,cc,k)

  for(j in 1:ll){    #iterate through subsample

    yb <- y[-((n1*(j-1)+1):(n1*j))]

    xb <- x[-((n1*(j-1)+1):(n1*j)),]

    yt <- y[(n1*(j-1)+1):(n1*j)]

    xt <- x[(n1*(j-1)+1):(n1*j),]

    nn <- length(yb)

    nsub <- nn*ratio

    gamma.lassob <- 0
    Deltab <- 0

    est <- matrix(NA,nrow = BB, ncol = k )

    for(b in 1:BB){

      index <- sample(1:nn,nsub)

      ybs <- yb[index]

      xbs <- xb[index,]

      ybf <- yb[-index]

      xbf <- xb[-index,]

      # adaptive lasso
      fit.ridge <- cv.glmnet(x = xbs, y=ybs, penalty.factor = penalty, alpha=0)

      beta.ridge <- coef(fit.ridge, s= "lambda.min")[-1]

      fit.lasso <- cv.glmnet(x = xbs,
                             y = ybs,
                             penalty.factor = penalty/abs(beta.ridge))

      # set bounds for model size
      modIdx <- (fit.lasso$nzero>(1+k))&(fit.lasso$nzero<(k+5+k))

      ss1 <- fit.lasso$lambda[modIdx]

      mcvError <- fit.lasso$cvm[modIdx] # mean cv error

      #the optimal tuning
      ss <- ss1[which.min(mcvError)]

      # selection set
      beta.lasso <- coef(fit.lasso, s = ss)[2:(p+1)]

      set1 <- setdiff(which(beta.lasso!=0), G)

      set1 <- sort(set1)

      #refit estimate
      refit_est <- lm(ybf~xbf[,G]+xbf[,set1])$coef[2:(k+1)]

      est[b,] <- refit_est

      ZZb <- cbind(xbf[,G],xbf[,set1])

      Delta0b <-  matrix(0,k,p)

      Delta0b[,union(G,set1)] <- solve(t(ZZb) %*% ZZb/(nn-nsub))[G,]

      Deltab <- Deltab + Delta0b
    }

    # filter results
    filtered_index <- apply( est, 2, function(x) order(x)[ceiling(( (1-0.95)*nrow(est))):(floor((1-0.05)*nrow(est)))])

    filtered_gamma <- NULL

    for(i in 1:k){

      filtered_gamma[i]<- mean(est[filtered_index[,i],i])
    }

    gamma.lassob <- filtered_gamma

    Deltab <- Deltab/BB

    fit.lassob <- cv.glmnet(x = xb, y = yb)

    beta.lassob <- coef(fit.lassob, s = "lambda.min")

    predb <- beta.lassob[1] + xb%*%beta.lassob[-1]

    residualb <- yb - predb

    epsilionb <- residualb-mean(residualb)

    TBB <- matrix(0, B, cc)

    correction <- matrix(0,cc,k) # correction matrix

    for(i in 1:cc){

      r0 <- r[i]

      correction[i,] <- (1-n^(r0-0.5))*(max(gamma.lassob)-gamma.lassob)
    }

    for(i in 1:B){

      Bepsilionb <- as.matrix(rnorm(nn)*epsilionb,nn,1)

      Bgamma.lassob <- gamma.lassob + Deltab %*% t(cbind(xb[,G],xb[,-G])) %*% Bepsilionb/nn

      for(l in 1:cc)
        TBB[i,l] <- max(Bgamma.lassob+correction[l,])-max(gamma.lassob)
    }

    brestimate <- 0

    for(l in 1:cc) #bias-reduced estimate
      brestimate[l] <- max(gamma.lassob)-mean(TBB[,l])

    gamma.lassott <- matrix(0,BB,k)

    nsub2 <- (n-nn)*ratio # rest of the subsamples

    Y.count <- matrix(data = NA, nrow = BB, ncol = n-nn)

    Delta0b <- 0

    est <- matrix(NA,nrow = BB, ncol = k )

    for(b in 1:BB){

      index <- sample(1:(n-nn),nsub2)

      yts <- yt[index]

      xts <- xt[index,]

      ytf <- yt[-index]

      xtf <- xt[-index,]

      # adaptive lasso
      fit.ridge <- cv.glmnet(x = xts, y=yts, penalty.factor = penalty, alpha=0)

      beta.ridge <- coef(fit.ridge, s= "lambda.min")[-1]

      fit.lasso <- cv.glmnet(x = xts,
                             y = yts,
                             penalty.factor = penalty/abs(beta.ridge))

      # set bounds for model size
      modIdx <- (fit.lasso$nzero>(1+k))&(fit.lasso$nzero<(k+5+k))

      ss1 <- fit.lasso$lambda[modIdx]

      mcvError <- fit.lasso$cvm[modIdx]

      # the optimal tuning
      ss <- ss1[which.min(tt1)]

      #selection set
      beta.lasso <- coef(fit.lasso, s = ss)[2:(p+1)]

      set1 <- setdiff(which(beta.lasso!=0), G)

      set1 <- sort(set1)

      #refit estimate
      refit_est <- lm(ytf~xtf[,G]+xtf[,set1])$coef[2:(k+1)]

      est[b,] <- refit_est

      Y.count[b,] <- tabulate(index,n-nn)
    }

    # filter results
    filtered_index <- apply( est, 2, function(x) order(x)[ceiling(( (1-0.95)*nrow(est))):(floor((1-0.05)*nrow(est)))])

    filtered_gamma <- NULL

    gamma.lassott <- matrix(NA, nrow = nrow(filtered_index), ncol = k)

    for(i in 1:k){

      gamma.lassott[,i] <- est[filtered_index[,i],i]

      filtered_gamma[i]<- mean(est[filtered_index[,i],i])
    }

    gamma.lassot <- filtered_gamma

    good_iter <- order(table(c(filtered_index)))[-c(1: (BB-nrow(filtered_index)))]

    Y.count <- Y.count[good_iter,]

    sd <- 0

    for(i in 1:k){

      sd[i] <- IFvarestbiascorr(Y.count, gamma.lassott[,i], n-nn, split.size = n33)
    }

    for(i in 1:cc)
      mseMat[i,] <- mseMat[i,]+(brestimate[i]-gamma.lassot)**2-sd**2
  }

  mse <- 0

  for(i in 1:cc)
    mse[i] <- min(mseMat[i,])

  op <- which.min(hh)

  return(op)
}

