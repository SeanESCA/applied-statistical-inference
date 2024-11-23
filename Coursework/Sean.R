library(rlang)
load(url("https://people.bath.ac.uk/kai21/ASI/data/CW24/dataQ1.RData"))


deriv_pack = function(func_expr, 
                      namevec = c("alpha", "beta", "theta"),
                      data_arg_vec = c("x", "y"),
                      theta_init = rep(1, 3)) {
  deriv_res <- deriv(
    func_expr,
    namevec = namevec,
    function.arg = c(namevec, data_arg_vec)
  )
  
  fn <- function(theta = theta_init, data_list = list("x" = 1, "y" = 1), apply.sum = T) {
    theta_list = as.list(theta)
    names(theta_list) = namevec
    aux  <- do.call(deriv_res, c(theta_list, data_list))
    if(apply.sum) {
      sum(as.numeric(aux))
    } else {
      c(aux)
    } 
  }
  
  gr <- function(theta = theta_init, data_list = list("x" = 1, "y" = 1), apply.sum = T) {
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

fit_optim <- function(par_len, deriv_pack, data, mean = 0, sd = 1, method = "BFGS", 
                      hessian = T, N_samples = 100, seed = 3141592){
  
  # set.seed(seed)
  # use set.seed for reproducibility purposes only. But always better to try without it and then settle at the end for reproducibility
  
  fit <- vector("list", length = N_samples)
  for (i in 1:N_samples){
    
    stop_crit <-T
    
    #tries until it finds an optimisation with no upfront error . this guarantees 'N_samples' optimisation are done
    
    while(stop_crit){
      
      # sometimes initial point is too far, which throws an straight error, hence the use of 'try'
      fit[[i]] <- try( 
        optim(par = rnorm(par_len, mean = mean, sd = sd),
              fn  = deriv_pack$fn,
              gr  = deriv_pack$gr,
              data_list = data,
              method  = method ,
              hessian = hessian),
        silent = F) # this supresses the red error messages
      
      if(inherits(fit[[i]], "try-error")){
        
        stop_crit <-T # if error, tries again
        
      } else {
        
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
  # selects the optimisation with minimum negative loglikelihood
  nll_vals <- lapply(X = fit, FUN = extract_negloglik)
  fit[[which.min(nll_vals)]] # return the final selected optimisation
}

# Question 1.1

## Model 1

### Reparameterise sigma to exp(theta)
nll_expr_q1_1 = expr(0.5 * log(2 * pi) + theta + alpha + beta * x + 0.5 * (y - exp(alpha + beta * x)) ^ 2 / exp(2 * (theta + alpha + beta * x)))
nll_pack_q1_1 = deriv_pack(nll_expr_1)

## Model 2
nll_expr_q1_2 = expr(0.5 * log(2 * pi) + theta + 0.5 * (log(y) - alpha - beta * x) ^ 2 / exp(2 * theta))
nll_pack_q1_2 = deriv_pack(nll_expr_2)

## Model 3
z_expr_q1_3 = expr((log(y) - alpha - beta * x) / exp(theta))
nll_expr_q1_3 = expr(!!z_expr_q1_3 + 2 * log((1 + exp(-!!z_expr_q1_3))) + theta)
nll_pack_q1_3 = deriv_pack(nll_expr_3)

## Model 4
nll_expr_q1_4 = expr(-log(log(2)) - theta + exp(theta) * (alpha + beta * x) + (1 - exp(theta)) * log(y) + log(2) * (y / exp(alpha + beta * x)) ^ exp(theta))
nll_pack_q1_4 = deriv_pack(nll_expr_4)

## Test
alpha = rnorm(1)
beta = rnorm(1)
theta = rnorm(1)
x = rnorm(1)
y = rnorm(1)

dnorm(y, exp(alpha + beta * x), exp(theta + alpha + beta * x), log = T)
nll_pack_q1_1$fn(c(alpha, beta, theta), list("x" = x, "y" = y))

dnorm(log(y), alpha + beta * x, exp(theta), log = T)
nll_pack_q1_2$fn(c(alpha, beta, theta), list("x" = x, "y" = y))

dlogis(log(y), alpha + beta * x, exp(theta), log = T)
nll_pack_q1_3$fn(c(alpha, beta, theta), list("x" = x, "y" = y))

## Question 1
opt_res_q1_1 = fit_optim(3, nll_pack_q1_1, dataQ1)
opt_res_q1_2 = fit_optim(3, nll_pack_q1_2, dataQ1)
opt_res_q1_3 = fit_optim(3, nll_pack_q1_3, dataQ1)
opt_res_q1_4 = fit_optim(3, nll_pack_q1_4, dataQ1)

nll_pack_q1_1$gr(opt_res_q1_1$par, dataQ1)
nll_pack_q1_2$gr(opt_res_q1_2$par, dataQ1)
nll_pack_q1_3$gr(opt_res_q1_3$par, dataQ1)
nll_pack_q1_4$gr(opt_res_q1_4$par, dataQ1)

aic_q1_1 = 2 * (nll_pack_q1_1$fn(opt_res_q1_1$par, dataQ1) + 3)
aic_q1_2 = 2 * (nll_pack_q1_2$fn(opt_res_q1_2$par, dataQ1) + 3)
aic_q1_3 = 2 * (nll_pack_q1_3$fn(opt_res_q1_3$par, dataQ1) + 3)
aic_q1_4 = 2 * (nll_pack_q1_4$fn(opt_res_q1_4$par, dataQ1) + 3)

# Model 1 is the best model as it has the lowest AIC.

# Question 1.2

plot_confint_q1_2 = function(opt_res, f = identity, main,
                             data = dataQ1, xlim = c(1, 15), ngrid = 100) {
  x <- seq(xlim[1], xlim[2],length = ngrid)
  S <- solve(opt_res$hessian) 
  ci <- matrix(NA,
               nrow = ngrid,
               ncol = 2)
  
  for (i in 1:ngrid){
    vec.x <- c(1, x[i])
    est <- crossprod(vec.x, opt_res$par[1:2])
    se <- sqrt(crossprod(vec.x, S[1:2, 1:2]) %*% vec.x)
    ci[i,] <- c(est - 1.96 * se, est + 1.96 * se) 
  }
  
  plot(data, pch = 19, main = main)
  lines(x, f(opt_res$par[1] + opt_res$par[2] * x), col = "red")
  lines(x, f(ci[,1]), col = "red", lty = 2)
  lines(x, f(ci[,2]), col = "red", lty = 2)
}

plot_confint_q1_2(opt_res_q1_1, exp, 
                  "Plot of the median for model 1 against x")
plot_confint_q1_2(opt_res_q1_2, 
                  main = "Plot of the median for model 2 against x")
plot_confint_q1_2(opt_res_q1_3,
                  main = "Plot of the median for model 3 against x")
plot_confint_q1_2(opt_res_q1_4, exp,
                  "Plot of the median for model 4 against x")

# Question 2
load(url("https://people.bath.ac.uk/kai21/ASI/data/CW24/dataQ2.RData"))

# Question 2.1
# nll_expr_q2_full = expr(lgamma(k) - k * log(lambda) + (1 - k) * log(y) + lambda * y)
# Reparameterise phi as exp(theta)

nll_expr_q2_full = expr(lgamma(exp(-theta)) - exp(-theta) * (log(alpha + beta * x) - theta) + (1 - exp(-theta)) * log(y) + exp(-theta) * (alpha + beta * x) * y)
nll_pack_q2_full = deriv_pack(nll_expr_q2_full)

nll_expr_q2_1 = expr(lgamma(exp(-theta1)) - exp(-theta1) * (log(theta1 - x * theta1 ^ 2) - theta1) + (1 - exp(-theta1)) * log(y) + exp(-theta1) * (theta1 - x * theta1 ^ 2) * y)
nll_pack_q2_1 = deriv_pack(nll_expr_q2_1,
                           c("theta1"),
                           theta_init = 2)

nll_expr_q2_2 = expr(lgamma(exp(-theta2)) - exp(-theta2) * (log(theta2 + tau2 * x) - theta2) + (1 - exp(-theta2)) * log(y) + exp(-theta2) * (theta2 + tau2 * x) * y)
nll_pack_q2_2 = deriv_pack(nll_expr_q2_2,
                           c("theta2", "tau2"),
                           theta_init = rep(1, 2))

nll_expr_q2_3 = expr(lgamma(exp(-theta4)) - exp(-theta4) * (log(theta3 - x * theta3 ^ 2) - theta4) + (1 - exp(-theta4)) * log(y) + exp(-theta4) * (theta3 - x * theta3 ^ 2) * y)
nll_pack_q2_3 = deriv_pack(nll_expr_q2_3,
                           c("theta3", "theta4"),
                           theta_init = rep(2, 2))

# Test
-dgamma(y, exp(-theta), (alpha + beta * x) * exp(-theta), log=T)
nll_pack_q2_full$fn(c(alpha, beta, theta), list("x" = x, "y" = y))

-dgamma(y, exp(-theta), (theta - x * theta ^ 2) * exp(-theta), log = T)
nll_pack_q2_1$fn(theta, list("x" = x, "y" = y))

-dgamma(y, exp(-theta), (theta + beta * x) * exp(-theta), log = T)
nll_pack_q2_2$fn(c(theta, beta), list("x" = x, "y" = y))

-dgamma(y, exp(-theta), (alpha - 1 * alpha ^ 2) * exp(-theta), log = T)
nll_pack_q2_3$fn(c(alpha, theta), list("x" = 1, "y" = y))

opt_res_q2_full = fit_optim(3, nll_pack_q2_full, dataQ2)
opt_res_q2_1 = fit_optim(1, nll_pack_q2_1, dataQ2)
opt_res_q2_2 = fit_optim(2, nll_pack_q2_2, dataQ2)
opt_res_q2_3 = fit_optim(2, nll_pack_q2_3, dataQ2)

nll_pack_q2_full$gr(opt_res_q2_full$par, dataQ2)
nll_pack_q2_1$gr(opt_res_q2_1$par, dataQ2)
nll_pack_q2_2$gr(opt_res_q2_2$par, dataQ2)
nll_pack_q2_3$gr(opt_res_q2_3$par, dataQ2)

nic = function(deriv_pack, opt_res, data = dataQ2) {
  aux = deriv_pack$gr(opt_res$par, data, apply.sum = F)
  K_hat = t(aux) %*% aux
  2 * (opt_res$value + sum(diag(solve(opt_res$hessian) %*% K_hat)))
}

nic_full = nic(nll_pack_q2_full, opt_res_q2_full)
nic_1 = nic(nll_pack_q2_1, opt_res_q2_1)
nic_2 = nic(nll_pack_q2_2, opt_res_q2_2)
nic_3 = nic(nll_pack_q2_3, opt_res_q2_3)

# The full model is the best as it has the lowest NIC.

# Question 2.2
ngrid = 100
x <- seq(0, 1,length = ngrid)
S_full <- solve(opt_res_q2_full$hessian)  
ci_full <- matrix(NA,
               nrow = ngrid,
               ncol = 2)

S_2 <- solve(opt_res_q2_2$hessian)  
ci_2 <- matrix(NA, nrow = ngrid, ncol = 2)

for (i in 1:ngrid){
  vec.x <- c(1, x[i])
  est_full <- crossprod(vec.x, opt_res_q2_full$par[1:2])
  se_full <- sqrt(crossprod(vec.x, S_full[1:2, 1:2]) %*% vec.x)
  ci_full[i,] <- c(est_full - 1.96 * se_full, est_full + 1.96 * se_full) 
  
  est_2 <- crossprod(vec.x, opt_res_q2_2$par[1:2])
  se_2 <- sqrt(crossprod(vec.x, S_2[1:2, 1:2]) %*% vec.x)
  ci_2[i,] <- c(est_2 - 1.96 * se_2, est_2 + 1.96 * se_2) 
}

est_1 <- opt_res_q2_1$par - x * opt_res_q2_1$par ^ 2
se_1 <- abs(1 - 2 * opt_res_q2_1$par * x) / c(opt_res_q2_1$hessian ^ 0.5)
ci_1 <- matrix(NA, nrow = ngrid, ncol = 2)
ci_1[, 1] <- est_1 - 1.96 * se_1
ci_1[, 2] <- est_1 + 1.96 * se_1

est_3 <- opt_res_q2_3$par[1] - x * opt_res_q2_3$par[1] ^ 2
se_3 <- sqrt(solve(opt_res_q2_3$hessian)[1, 1] * (1 - 2 * opt_res_q2_3$par[1] * x) ^ 2)
ci_3 <- matrix(NA, nrow = ngrid, ncol = 2)
ci_3[, 1] <- est_3 - 1.96 * se_3
ci_3[, 2] <- est_3 + 1.96 * se_3

plot_confint_q1_2(opt_res_q2_full, 
                  f = function(x) 1/x, 
                  main = "",
                  data = dataQ2, 
                  xlim = c(0, 1))

plot(dataQ2, pch = 19)
lines(x, 1 / (opt_res_q2_1$par - x * opt_res_q2_1$par ^ 2), col = "red")
lines(x, 1 / ci_1[,1], col = "red", lty = 2)
lines(x, 1 / ci_1[,2], col = "red", lty = 2)

plot_confint_q1_2(opt_res_q2_2, 
                  f = function(x) 1/x, 
                  main = "",
                  data = dataQ2, 
                  xlim = c(0, 1))

plot(dataQ2, pch = 19)
lines(x, 1 / (opt_res_q2_3$par[1] - x * opt_res_q2_3$par[1] ^ 2), col = "red")
lines(x, 1 / ci_3[,1], col = "red", lty = 2)
lines(x, 1 / ci_3[,2], col = "red", lty = 2)

# Question 2.3