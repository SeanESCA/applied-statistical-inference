# Applied Statistical Inference

## Useful Functions

```{r nll_init}
library(rlang)
library(glue)

deriv_pack = function(func_expr, 
                    namevec = c("theta1", "theta2", "theta3"),
                    data_arg_vec = c("y"),
                    theta_init = rep(1, 3)) {
  deriv_res <- deriv(
    func_expr,
    namevec = namevec,
    function.arg = c(namevec, data_arg_vec)
  )
  
  fn <- function(theta = theta_init, data_list = list("y" = 1), apply.sum = T) {
    theta_list = as.list(theta)
    names(theta_list) = namevec
    aux  <- do.call(deriv_res, c(theta_list, data_list))
    if(apply.sum) {
      sum(as.numeric(aux))
    } else {
      c(aux)
    } 
  }
  
  gr <- function(theta = theta_init, data_list = list("y" = 1), apply.sum = T) {
    theta_list = as.list(theta)
    names(theta_list) = namevec
    aux  <- do.call(deriv_res, c(theta_list, data_list))
    if(apply.sum) {
      apply(attr(aux, "gradient"), 2, sum)
    } else {
      attr(aux, "gradient")
    }
  }
  
  list(fn = fn, gr = gr)
}
```

```{r opt}
opt = function(nll, par = rep(1, 4), data_list = list("y" = 1)) {
  optim(par = par,
        fn = nll$fn,
        gr = nll$gr,
        method = "BFGS",
        hessian = TRUE,
        data_list = data_list)
}
```

```{r glrt}
# Requires opt function.
glrt_test <- function(params_h0, params_h1, data, df = 1, alpha = 0.05) {
  opt_res_h0 <- do.call(opt, params_h0)
  opt_res_h1 <- do.call(opt, params_h1)
  test_stat = 2 * (opt_res_h0$value - opt_res_h1$value)
  p = pchisq(test_stat, df = df, lower.tail = F)
  c(signif(test_stat, 3), p)
}
```

```{r fit_optim}
fit_optim<- function(par       = c(0,0,0,0),
                     fn ,
                     gr ,
                     method = "BFGS",
                     hessian   = T,
                     y,
                     sd        =  c(1,1,1,1),
                     N_samples = 100,
                     seed      = 3141592){
  
  set.seed(seed)
  # use set.seed for reproducibility purposes only. But always better to try without it and then settle at the end for reproducibility
  
  fit <- vector("list",
                length = N_samples)
  
  for (i in 1:N_samples){
    
    stop_crit <-T
    
    #tries until it finds an optimisation with no upfront error . this guarantees 'N_samples' optimisation are done
    
    while(stop_crit){
  
       # sometimes initial point is too far, which throws an straight error, hence the use of 'try'
       fit[[i]]<-try( 
            optim(par = c(rnorm(1,mean = par[1],sd = sd[1]),
                          rnorm(1,mean = par[2],sd = sd[2]),
                          rnorm(1,mean = par[3],sd = sd[3]),
                          rnorm(1,mean = par[4],sd = sd[4])),
                  fn  = fn,
                  gr  = gr,
                  y =y,
                  method  = method ,
                  hessian = hessian),
                      silent=T) # this supresses the red error messages
           
      if(inherits(fit[[i]], "try-error")){
      
        stop_crit <-T # if error, tries again
      
      }else{
        
         stop_crit <-F # if no error, continues to next step
      }
    }
    
  
      
    # check for numerical convergence first 
    no_convergence <- fit[[i]]$convergence > 0
    
    # checks if asymptotic variances are possible to obtain
    no_variance <- inherits(try(solve(fit[[i]]$hessian),
                                silent = T), 
                            "try-error")
  
    null_variance <- F
    NA_variance   <- F
  
    if (!no_variance){
      # checks if asymptotic variance are NaN
      NA_variance <- as.logical(sum(is.nan(diag(solve(fit[[i]]$hessian)))))
    
    if(!NA_variance){
      # checks if asymptotic variance are zero up to machine precision
      null_variance <- as.logical(sum(diag(solve(fit[[i]]$hessian))< .Machine$double.eps ^ 0.5))
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
  # selects the optimisationm with minimum negative loglikelihood
  nll_vals <- lapply(X = fit, FUN = extract_negloglik)
  fit[[which.min(nll_vals)]] # return the final selected optimisation
}
```

