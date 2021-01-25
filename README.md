# subDebiased
This package implements bootstrap-assisted desparsified Lasso and bootstrap-assisted R-split estimators on selected subgroup's treatment effect estimation. The implemented estimators remove the subgroup selection bias and the regularization bias induced by high-dimensional covariates. 

- BSDesparseLasso implementes the bootstrap-assisted desparsified Lasso method
- BSSplitLasso implements the bootstrap-assisted R-Split method

The optimal tuning parameter for each method are selected via cross-validation:

- cvDesparse for desparsified Lasso
- cvSplit for R-Split

## References

- Zhang, C.-H. and Zhang, S. S. (2014). Confidence intervals for low dimensional parameters in
high dimensional linear models. Journal of the Royal Statistical Society: Series B (Statistical
Methodology), 76(1):217-242.

- Zhang, X. and Cheng, G. (2017). Simultaneous inference for high-dimensional linear models. Jour-
nal of the American Statistical Association, 112(518):757-768.

- Wang, J., He, X., and Xu, G. (2019). Debiased inference on treatment effect in a high dimensional
model. Journal of the American Statistical Association.

