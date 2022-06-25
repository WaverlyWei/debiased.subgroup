# debiased.subgroup

To install this package in R, run the following commands:

```R
library(devtools) 
install_github("WaverlyWei/debiased.subgroup")
```

This package implements bootstrap-assisted desparsified Lasso and bootstrap-assisted R-split estimators on selected subgroup's treatment effect estimation. The implemented estimators remove the subgroup selection bias and the regularization bias induced by high-dimensional covariates. 

- BSDesparseLasso implementes the bootstrap-assisted desparsified Lasso method
- BSSplitLasso implements the bootstrap-assisted R-Split method

The optimal tuning parameter for each method are selected via cross-validation:

- cvDesparse for desparsified Lasso
- cvSplit for R-Split

## Example usage:

```R
library(debiased.subgroup)
p <- 200
n <- 100 
ngroups <- 2 
s0 <- 4
m <- ngroups
Sigma <- matrix(0,p,p)
for (i in 1:n){
      for(j in 1:p){
        Sigma[i,j] <- 0.5^(abs(i-j))}
}

X <- mvrnorm( n = n, mu = rep(0,p), Sigma = Sigma )
Z <- matrix(0,n,m)
for(i in 1:n){
      for(j in 1:m){
        Z[i,j] <- rbinom(1,1,exp(X[i,2*j-1]   +X[i,2*j])/(1+exp(X[i,2*j-1]+ X[i,2*j])))}
}
noise.y <- 1
betas <- 1
w.index <- seq(1, m, 1) 
x <- cbind(Z,X)

beta <- c(rep(0,m-1),betas) 
gamma <- c(rep(1, s0), rep(0, p-s0)) 
noise <- mvrnorm( n = 1, mu = rep(0,n), Sigma = diag(n) * noise.y )

Y <- 0.5 + x %*% beta + noise
r <- 1/(3*1:10)
    
desparse_res <- BSDesparseLasso(y = Y,
                                x = x, 
                                r = r, 
                                G = w.index,
                                B = 5)
desparse_res
    
rsplit_res <- BSSplitLasso(y = Y,
                           x = x, 
                           r = r, 
                           G = w.index,
                           B = 5, BB = 10)
rsplit_res                     
```

### References
- Wang, J., He, X., and Xu, G. (2019). Debiased inference on treatment effect in a high dimensional model. Journal of the American Statistical Association.

- Zhang, C.-H. and Zhang, S. S. (2014). Confidence intervals for low dimensional parameters in high dimensional linear models. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 76(1):217-242.

- Zhang, X. and Cheng, G. (2017). Simultaneous inference for high-dimensional linear models. Journal of the American Statistical Association, 112(518):757-768.


