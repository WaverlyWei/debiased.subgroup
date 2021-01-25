#' Bootstrap-calibrated R-split method
#'
#' @param y response
#' @param x design matrix
#' @param r tuning parameter
#' @param G subgroup indicator
#' @param B bootstrap number
#' @param BB split number
#' @param alpha level  ## change other places
#' @param splitRatio split ratio
#' @param fold cross validation fold
#' @return
#' \item{LowerBound}{lower confidence bound}
#' \item{UpperBound}{upper confidence bound}
#' \item{op}{optimal tuning}
#' @export
BSSplitLasso <- function(y, x,
                         r = NULL,
                         G = NULL,
                         B = NULL,
                         BB = NULL,
                         alpha = 0.95,
                         splitRatio = 0.6,
                         fold = 2){

  if(is.null(r)){

    stop("Tuning parameter is missing.")

  }else if(is.null(G)){

    stop("Number of subgroups is missing.")

  }else if(is.null(B)){

    stop("Specify bootstrap iterations.")

  }else if(is.null(BB)){

    stop("Specify bootstrap iterations.")
  }

  p <- length(x[1,])

  n <- length(y)

  k <- length(G)

  cc <- length(r)

  nsub <- round(n*splitRatio)

  #only select the confounders
  penalty <- rep(0,p)

  penalty[-G] <- 1

  y <- y - mean(y)

  for(i in 1:p)
    x[,i] <- x[,i]-mean(x[,i])

  beta.lasso <- 0 # change this to beta?

  Delta <- 0

  # keep track of model size
  modelSize <- NULL

  est <- matrix(NA,nrow = BB, ncol = k )

  for(b in 1:BB) {

    #random split
    index <- sample(1:n,nsub)

    #selection sample
    ytrain <- y[index]

    xtrain <- x[index,]

    #refitting sample  # more natural choice
    ytest <- y[-index]

    xtest <- x[-index,] #120

    #adaptive lasso
    fit.ridge <- cv.glmnet(x = xtrain,
                           y=ytrain,
                           penalty.factor = penalty, alpha=0)

    beta.ridge <- coef(fit.ridge, s= "lambda.min")[-1]

    fit.lasso <- cv.glmnet(x = xtrain, y = ytrain, penalty.factor = penalty/abs(beta.ridge))

    # model size bounds
    modIdx <- (fit.lasso$nzero>(1+k))&(fit.lasso$nzero<(k+5+k))

    ss1 <- fit.lasso$lambda[modIdx]

    mcvError <- fit.lasso$cvm[modIdx]

    # selection set
    ss <- ss1[which.min(mcvError)]

    # selection set
    gamma.lasso <- coef(fit.lasso, s = ss)[2:(p+1)] # remove intercept

    set1 <- setdiff(which(gamma.lasso!=0), G)

    set1 <- sort(set1)

    modelSize[b] <- length(set1) # record model size

    #refit estimate
    refit_est <- lm(ytest~xtest[,G]+xtest[,set1])$coef[2:(k+1)]

    est[b,] <- refit_est

    # collect refit estimate,and take out large coeffecients
    # add a filter here
    ZZ <- cbind(xtest[,G],xtest[,set1])

    Delta0 <- matrix(0,k,p)

    Delta0[,union(G,set1)] <- solve(t(ZZ)%*%ZZ/(n-nsub))[G,]

    Delta <- Delta + Delta0
  }

  # filter results
  filtered_index <- apply( est, 2, function(x) order(x)[ceiling(( (1-0.95)*nrow(est))):(floor((1-0.05)*nrow(est)))])

  filtered_beta <- NULL

  for(i in 1:k){
    filtered_beta[i]<- mean(est[filtered_index[,i],i])
  }

  Delta <- Delta/BB

  beta.lasso <- filtered_beta

  #get the residual for boostrap
  fit.lasso <- cv.glmnet(x = x, y = y)

  gamma.lasso <- coef(fit.lasso, s = "lambda.min")

  pred <- gamma.lasso[1] + x%*%gamma.lasso[-1]

  residual <- y - pred

  epsilion <- residual-mean(residual)

  TB <- matrix(0, B, cc+1)

  correction <- matrix(0, cc+1, k)

  for(i in 1:cc) {

    r0 <- r[i]

    correction[i,] <- (1-n^(r0-0.5))*(max(beta.lasso)-beta.lasso)
  }

  # simultaneous one for R-split
  correction[cc+1,] <- (max(beta.lasso)-beta.lasso)

  TB_op <- matrix(0, B, cc)

  c_op <- matrix(0, cc, k)

  for(i in 1:cc) {

    r0 <- r[i]

    rp <- r0/sqrt(k/2) # NOTE: change to r.p

    c_op[i,] <- (1-n^(rp-0.5))*(max(beta.lasso)-beta.lasso)
  }

  for(i in 1:B) {

    #generate bootstrap estimate
    Bepsilion <- as.matrix(rnorm(n)*epsilion,n,1)

    Bbeta.lasso <- beta.lasso+Delta %*% t(cbind(x[,G],x[,-G])) %*% Bepsilion/n

    #correct maximum quantity
    for(j in 1:(cc+1)){

      TB[i,j] <- max(Bbeta.lasso+correction[j,])-max(beta.lasso)
    }

    for(j in 1:(cc)){
      TB_op[i,j] <- max(Bbeta.lasso+c_op[j,])-max(beta.lasso)
    }
  }


  op <- cvSplit(y, x, r, G, B, BB, splitRatio, fold)

  result <- list()

  for(j in 1:(cc+1)) {
    result[j] <- list(c(BSciCoverfun(beta.lasso, TB[,j], G, alpha),
                        #betaEst = list(beta.lasso),
                        op = r[op]))
  }

  if(is.integer(op) && length(op)==1){

    result[j+1] <- list(c(BSciCoverfun(beta.lasso, TB_op[,op],G, alpha),
                          #betaEst = list(beta.lasso),
                          op = r[op]))
  }else{
    result[j+1] = list(c(BSciCoverfun(beta.lasso, TB[,cc], G, alpha),
                         #betaEst = list(beta.lasso),
                          op = r[op]))
  }
  return(result[[12]])
}

