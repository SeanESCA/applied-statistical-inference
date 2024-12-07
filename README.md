# Applied Statistical Inference

## Useful Functions

### Unpacking `deriv` Functionality

```{r deriv_pack}
deriv_pack = function(func_expr, 
                      namevec = c("alpha", "beta", "theta"),
                      data_arg_vec = c("x", "y"),
                      theta_init = rep(1, 3)) {
  # Returns functions that separately evaluate the attributes of deriv.
  # Inputs: 
  # - func_expr: Expr passed to deriv function.
  # - namevec: Non-empty str vector that contains the names of the unknown
  #            parameters, passed to deriv function.
  # - data_arg_vec: Non-empty string vector that contains the names of the data 
  #                 variables, passed to deriv function.
  # - theta_init: Numerical vector of equal length to namevec that represents
  #               the default argument of the functions in the output.
  # Output: List of length 2 containing two functions:
  # - fn: Function that evaluates func_expr.
  # - gr: Function that evaluates the gradient of func_expr.
  # - hess: Function that evaluates the hessian of func_expr.
  
  deriv_res <- deriv(
    func_expr,
    namevec = namevec,
    function.arg = c(namevec, data_arg_vec)
  )
  
  fn <- function(theta = theta_init, data_list = list("x" = 1, "y" = 1), 
                 apply.sum = T) {
    # Evaluates func_expr.
    # Input:
    # - theta: Numerical vector representing the parameters.
    # - data_list: Data.frame or list representing the data set.
    # - apply.sum: Logical representing whether the results from evaluating
    #              func_expr on each data point in data_list are summed.
    # Output: if apply.sum = T, then a float is returned; otherwise, a matrix
    # is returned.
    theta_list = as.list(theta)
    names(theta_list) = namevec
    aux  <- do.call(deriv_res, c(theta_list, data_list))
    if(apply.sum) {
      sum(as.numeric(aux))
    } else {
      c(aux)
    } 
  }
  
  gr <- function(theta = theta_init, data_list = list("x" = 1, "y" = 1), 
                 apply.sum = T) {
    # Evaluates the gradient of func_expr.
    # Input:
    # - theta: Numerical vector representing the parameters.
    # - data_list: Data.frame or list representing the data set.
    # - apply.sum: Logical representing whether the results from evaluating
    #              func_expr on each data point in data_list are summed.
    # Output: if apply.sum = T, then a float is returned; otherwise, a matrix
    # is returned.
    theta_list = as.list(theta)
    names(theta_list) = namevec
    aux  <- do.call(deriv_res, c(theta_list, data_list))
    if(apply.sum) {
      apply(attr(aux, "gradient"), 2, sum)
    } else {
      attr(aux, "gradient")
    }
  }
  
  hess <- function(theta = theta_init, data_list = list("x" = 1, "y" = 1)) {
    theta_list = as.list(theta)
    names(theta_list) = namevec
    aux  <- do.call(deriv_res, c(theta_list, data_list))
    apply(attr(aux,"hessian"), c(2,3), sum)
  }
  
  list(fn = fn, gr = gr, hess = hess)
}
``` 

### Shorthand `optim`

```{r opt}
opt = function(nll_pack, par = rep(1, 3), data_list = list("y" = 1)) {
  # Shorthand function for optim with method = "BFGS" and hessian = T.
  optim(par = par,
        fn = nll_pack$fn,
        gr = nll_pack$gr,
        method = "BFGS",
        hessian = T,
        data_list = data_list)
}
```

### Running Multiple Optimisations

```{r fit_optim}
fit_optim <- function(par_len, nll_pack, data, mean = 0, sd = 10, 
                      method = "BFGS", hessian = T, N_samples = 100, 
                      silent = T){
  # Runs N_samples optimisations and returns the best result.
  # Inputs: 
  # - par_len: Int that represents the length of the unknown parameter vector
  # - nll_pack: List that contains two functions:
  #    - fn: Function to be minimised.
  #    - gr: Function that evaluates the gradient of fn.
  # - data: Data.frame or list that represents the data set.
  # - mean: Float that represents the mean of the normal distribution the 
  #         initial points for the optimisations are generated from.
  # - sd: Float that represents the sd of the normal distribution the 
  #       initial points for the optimisations are generated from.
  # - method: Str that represents the optimisation method to be used.
  # - hessian: Logical that represents whether the hessian is estimated for all
  #            the optimisations
  # - N_samples: Int that represents the number of optimisations to run.
  # - silent: Logical that represents whether error messages are suppressed.
  # Outputs: the same output as optim
  
  fit <- vector("list", length = N_samples)
  for (i in 1:N_samples){
    
    stop_crit <-T
    
    # Attempt an optimisation until there is no upfront error. 
    # This guarantees that N_samples optimisations are done.
    while(stop_crit){
      
      # try is used as the initial point may occassionally be too far, which 
      # throws a straight error.
      fit[[i]] <- try( 
        optim(par = rnorm(par_len, mean = mean, sd = sd),
              fn  = nll_pack$fn,
              gr  = nll_pack$gr,
              data_list = data,
              method  = method ,
              hessian = hessian),
        silent = silent)
      
      if(inherits(fit[[i]], "try-error")){
        
        stop_crit <-T # if there is an error, try again
        
      } else {
        
        stop_crit <-F # if there is no error, continue to the next step.
      }
    }
    
    # Check if there is no numerical convergence.
    no_convergence <- fit[[i]]$convergence > 0
    
    # Check if asymptotic variances are obtainable.
    no_variance <- inherits(try(solve(fit[[i]]$hessian),
                                silent = T), 
                            "try-error")
    
    null_variance <- F
    NA_variance   <- F
    
    if (!no_variance){
      #  Check if any of the asymptotic variances are NaN.
      NA_variance <- as.logical(sum(is.nan(diag(solve(fit[[i]]$hessian)))))
      
      if(!NA_variance){
        # Check if any of the asymptotic variances are zero up to machine precision
        null_variance <- as.logical(sum(diag(solve(fit[[i]]$hessian)) < 
                                          .Machine$double.eps ^ 0.5))
      }
    }
    
    fail <- no_variance | no_convergence | NA_variance | null_variance 
    if (fail){
      fit[[i]]$value <- NA
    }
    
  } 
  
  extract_negloglik <- function(optim_object){
    optim_object$value
  }
  # Select the optimisation with minimum negative log-likelihood.
  nll_vals <- lapply(X = fit, FUN = extract_negloglik)
  
  # Return the selected optimisation.
  fit[[which.min(nll_vals)]] 
}
```

### GLRT

```{r glrt}
glrt_test <- function(params_h0, params_h1, data, df = 1, alpha = 0.05) {
  # Requires opt function.
  opt_res_h0 <- do.call(opt, params_h0)
  opt_res_h1 <- do.call(opt, params_h1)
  test_stat = 2 * (opt_res_h0$value - opt_res_h1$value)
  p = pchisq(test_stat, df = df, lower.tail = F)
  c(signif(test_stat, 3), p)
}
```

```{r nic}
nic = function(nll_pack, opt_res, data) {
  # Computes the NIC.
  # Input: 
  # - nll_pack: deriv_pack output.
  # - opt_res: optim output.
  # - data: Data.frame or list representing the data set.
  # Output: Float representing the NIC.
  
  aux = nll_pack$gr(opt_res$par, data, apply.sum = F)
  2 * (opt_res$value + sum(diag(solve(opt_res$hessian) %*% crossprod(aux))))
}
```

### Delta Method

> [!TIP]
> The `emdbook` [package](https://cran.r-project.org/web/packages/emdbook/index.html) 
> has the functions `deltamethod` and `deltavar`  for implementing the delta 
> method. It is believed that the code below agrees with those functions.

```{r mean_delta_method}
mean_delta_method = function(bands, mean_expr, namevec, par, cov_mat) {
  # Calculates the mean and its 95% CI using the delta method.
  # Inputs:
  # - bands: Matrix with the first three columns representing the lower bound,
  #          of the CI, the estimated mean, and the upper bound of the CI.
  # - mean_expr: Expr of the mean function.
  # - namevec: Non-empty str vector containing the names of the estimated
  #            parameters in the mean function.
  # - par: Numerical vector containing the estimated parameter values.
  # - cov_mat: Covariance matrix for the estimated parameters.
  # Output: bands matrix
  
  temp = bands
  mean_pack = deriv_pack(mean_expr, namevec, c("x"))
  for (i in 1:n_grid) {
    est = mean_pack$fn(par, list("x" = x[i]))
    J = mean_pack$gr(par, list("x" = x[i]))
    se = sqrt(J %*% cov_mat %*% J)
    temp[i,] = est + c(- 1.96 * se, 0, 1.96 * se)
  }
  temp
}
```

```{r set_P22_bands}
set_P22_bands = function(P22_bands, mean_expr, namevec, opt_res) {
  # Calculates the mean and its 95% CI using the delta method with the hessian  
  # of the NLL as the covariance matrix.
  # Inputs:
  # - bands: P22_bands_full, P22_bands_1, P22_bands_2, P22_bands_3 matrix. 
  # - mean_expr: Expr of the mean function.
  # - namevec: Non-empty str vector containing the names of the estimated
  #            parameters in the mean function.
  # - opt_res: optim output.
  # Output: bands matrix
  
  mean_delta_method(P22_bands, mean_expr, namevec, opt_res$par, 
                    solve(opt_res$hessian))
}
```

#### Misspecified Models

```{r P23_bands}
P23_bands = function(bands, mean_expr, namevec, nll_pack, opt_res) {
  # Calculates the mean and its 95% CI under a misspecified model using the 
  # delta method.
  # Inputs:
  # - bands: P22_bands_full, P22_bands_1, P22_bands_2, P22_bands_3 matrix. 
  # - mean_expr: Expr of the mean function.
  # - namevec: Non-empty str vector containing the names of the estimated
  #            parameters in the mean function.
  # - nll_pack: deriv_pack output.
  # - opt_res: optim output.
  # Output: bands matrix
  
  # Calculate the estimated covariance matrix.
  S = solve(opt_res$hessian)
  nll_gr = nll_pack$gr(opt_res$par, dataQ2, apply.sum = F)
  Khat = crossprod(nll_gr)
  cov_mat = S %*% Khat %*% S
  
  # Return the mean and its 95% CI.
  mean_delta_method(bands, mean_expr, namevec, opt_res$par, cov_mat)
}
```

### KL Divergence

```{r KL_pack}
KL_pack = function(nll_pack) {
  # Defines functions for computing the KL divergence and its gradient.
  # Input: 
  # - nll_pack: Output from deriv_pack
  # Output: List of two functions:
  # - fn: Function that evaluates the KL divergence.
  # - gr: Function that calculates the gradient of the KL divergence.
  
  fn_integrand = function(y, x, theta) {
    # Evaluates the expression that is integrated with respect to y to 
    # calculate the KL divergence.
    # Input:
    # - y: Numerical vector
    # - x: Float
    # - theta: Numerical vector representing the parameter vector
    # Output: Float
    
    (nll_pack$fn(theta, list("x" = x, "y" = y), apply.sum = F) +
      dgamma(y, 2, (4 - 2 * x)/3, log = T)) * 
      dgamma(y, 2, (4 - 2 * x)/3)
  }
  gr_integrand = function(y, x, theta, col_int = 1) {
    # Evaluates the expression that is integrated with respect to y 
    # to calculate the gradient of the KL divergence.
    # Input:
    # - y: Numerical vector
    # - x: Float
    # - theta: Numerical vector representing the parameter vector
    # - col_int: Int representing the parameter the gradient is 
    #            calculated with respect to.
    # Output: Float
    aux = nll_pack$gr(theta, list("x" = x, "y" = y), apply.sum = F) *
      dgamma(y, 2, (4 - 2 * x)/3)
    aux[, col_int]
  }
  KL_fn = function(theta, data_list = dataQ2) {
    # Evaluates the KL divergence.
    # Input: 
    # - theta: Numerical vector representing the parameter vector.
    # - data_list: Unused argument for compatibility with fit_optim.
    # Output: Float
    
    res = 1:10
    for(i in res) {
      res[i] = integrate(fn_integrand, 0, 10, x = i/10, 
                         theta = theta)$value
    }
    sum(res)
  }
  
  KL_gr = function(theta, data_list = dataQ2) {
    # Calculates the gradient of the KL divergence.
    # Input: 
    # - theta: Numerical vector representing the parameter vector.
    # - data_list: Unused argument for compatibility with fit_optim.
    # Output: Float
    
    n = length(theta)
    res = matrix(NA, nrow = 10, ncol = n)
    for(i in 1:10) {
      for(j in 1:n) {
        res[i, j] = integrate(gr_integrand, 0, 10, x = i/10, 
                              theta = theta, col_int = j)$value
      }
    }
    apply(res, 2, sum)
  }
  
  list(fn = KL_fn, gr = KL_gr)
}
```

### Coursework Functions

```{r Q1.2}
n_grid     <- 1000 
x          <- seq(1, 15, length=n_grid)
P2_bands_1 <- P2_bands_2 <- 
  P2_bands_3 <- P2_bands_4 <- matrix(NA,nrow=n_grid,ncol=3)
colnames(P2_bands_1) <- colnames(P2_bands_2) <- 
  colnames(P2_bands_3) <- colnames(P2_bands_4) <- c("lower","est","upper")
head(P2_bands_1)
```

```{r set P2_bands}
set_P2_bands = function(P2_bands, opt_res) {
  # Calculates the mean and its 95% CI for Q1.2.
  # Inputs:
  # - P2_bands: P2_bands_1, P2_bands_2, P2_bands_3, or P2_bands_4 matrix.
  # - opt_res: optim output
  # Output: P2_bands
  temp = P2_bands
  S <- solve(opt_res$hessian) 
  for (i in 1:n_grid){
    vec.x <- c(1, x[i])
    est <- crossprod(vec.x, opt_res$par[1:2])
    se <- sqrt(crossprod(vec.x, S[1:2, 1:2]) %*% vec.x)
    temp[i,] <- exp(c(est - 1.96 * se, est, est + 1.96 * se))
  }
  temp
}
```