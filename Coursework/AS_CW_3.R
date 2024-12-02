---
  title: "MA40198 Coursework"
author: "Group 13: Sean Soon and Shaveen Colambage"
date: "`r Sys.Date()`"
output: html_document
---
  
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Libraries 

```{r}
library(rlang)
```



# Part I

```{r}
load(url("https://people.bath.ac.uk/kai21/ASI/data/CW24/dataQ1.RData"))
```

## Question 1.1


```{r Q1.1}
P2_optim_1 <- P2_optim_2 <- P2_optim_3 <- P2_optim_4 <- vector("list", 5)
names(P2_optim_1) <- names(P2_optim_2) <- names(P2_optim_3) <- names(P2_optim_4) <- c("mle","negloglik","gradient","hessian","AIC")

P2_optim_1
```

## Here we create a function "deriv_pack" which for a function expression,
## given intial function paramaters and data values computes: 
## "fn" - the value of the function
## "gr" - the value of the gradient of the function
```{r deriv_pack}
deriv_pack = function(func_expr, 
                      namevec = c("alpha", "beta", "theta"),
                      data_arg_vec = c("x", "y"),
                      theta_init = rep(1, 3)) {
  deriv_res <- deriv(
    func_expr,
    namevec = namevec,
    function.arg = c(namevec, data_arg_vec)
  )
  
  fn <- function(theta = theta_init, data_list = list("x" = 1, "y" = 1), 
                 apply.sum = T) {
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
## Here we create a function 'fit_optim' which takes the arguments:
## "par_len" - the number of parameters, "nll_pack" - the output
## of 'deriv_pack' when the function is passed through it and "data"
## for the dataset. It returns an MLE estimate of the parameters as well 
## for the gradient and hessian after conducting 100 iterations with different
## intial values to determine the best MLE estimate
```{r fit_optim}
fit_optim <- function(par_len, nll_pack, data, mean = 0, sd = 10, 
                      method = "BFGS", hessian = T, N_samples = 100, 
                      silent = T, seed = 3141592){
  
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
              fn  = nll_pack$fn,
              gr  = nll_pack$gr,
              data_list = data,
              method  = method ,
              hessian = hessian),
        silent = silent) # this supresses the red error messages
      
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
```
## Here we create a function "aic" which for inputs:
## "nll_pack" and "opt_res", where "opt_res" is the MLE estimate obtained
## from 'fit_optim'
## Returns the AIC value of the model
```{r aic}
aic = function(nll_pack, opt_res, data = dataQ1) {
  par = opt_res$par
  2 * (nll_pack$fn(par, data) + length(par))
}
```
## The function 'set_P2_optim' for given arguments:
## "P2_optim" - the optimsation of the function, "nll_pack" and "opt_res"

## Outputs the values of the MLE estimate of the parameters; the negative
## log-likelihood value; the value of the hessian; the 
## values of the gradient and the AIC

```{r}
set_P2_optim = function(P2_optim, nll_pack, opt_res) {
  P2_optim$mle = c(opt_res$par[1:2], exp(opt_res$par[3]))
  P2_optim$negloglik = opt_res$value
  P2_optim$hessian = opt_res$hessian
  P2_optim$gradient = nll_pack$gr(opt_res$par, dataQ1)
  P2_optim$AIC = aic(nll_pack, opt_res)
  P2_optim
}
```
## For all four models, $sigma$ was reparameterised as $e^{\theta}$.
## Here we have:
## 'nll_expr_q1_1', 'nll_pack_q1_1', 'opt_res_q1_1' as the 
## expression, 'deriv_pack' and 'fit_optim' for the normal distribution 
## model with median(or mean) of (alpha1 + beta1*x) and
## standard deviation of sigma1*exp(alpha1+beta1*x)
## alpha1 is written as 'alpha', beta1 is written as 'beta' and sigma1 is 
## written as 'theta'

```{r P2_optim_1, message=FALSE, warning=FALSE}
nll_expr_q1_1 = expr(0.5 * log(2 * pi) + theta + alpha + beta * x + 
                       0.5 * (y - exp(alpha + beta * x)) ^ 2 / 
                       exp(2 * (theta + alpha + beta * x)))
nll_pack_q1_1 = deriv_pack(nll_expr_q1_1)
opt_res_q1_1 = fit_optim(3, nll_pack_q1_1, dataQ1)
set_P2_optim(P2_optim_1, nll_pack_q1_1, opt_res_q1_1)
```
## Here we have:
## 'nll_expr_q1_2', 'nll_pack_q1_2', 'opt_res_q1_2' as the 
## expression, 'deriv_pack' and 'fit_optim' for the normal distribution 
## model with median(or mean) of (alpha2 + beta2*x) and
## standard deviation of sigma2
## alpha2 is written as 'alpha', beta2 is written as 'beta' and sigma2 is 
## written as 'theta'

```{r P2_optim_2, message=FALSE, warning=FALSE}
nll_expr_q1_2 = expr(0.5 * log(2 * pi) + theta + 
                       0.5 * (log(y) - alpha - beta * x) ^ 2 / exp(2 * theta))
nll_pack_q1_2 = deriv_pack(nll_expr_q1_2)
opt_res_q1_2 = fit_optim(3, nll_pack_q1_2, dataQ1)
set_P2_optim(P2_optim_2, nll_pack_q1_2, opt_res_q1_2)
```
## Here we have:
## 'nll_expr_q1_3', 'nll_pack_q1_3', 'opt_res_q1_3' as the 
## expression, 'deriv_pack' and 'fit_optim' for the logistic distribution 
## model with median(or mean) of (alpha3 + beta3*x) and
## scale paramter of sigma3
## alpha3 is written as 'alpha', beta3 is written as 'beta' and sigma3 
## is written as 'theta'

```{r P2_optim_3, message=FALSE, warning=FALSE}
z_expr_q1_3 = expr((log(y) - alpha - beta * x) / exp(theta))
nll_expr_q1_3 = expr(!!z_expr_q1_3 + 2 * log((1 + exp(-!!z_expr_q1_3))) + theta)
nll_pack_q1_3 = deriv_pack(nll_expr_q1_3)
opt_res_q1_3 = fit_optim(3, nll_pack_q1_3, dataQ1)
set_P2_optim(P2_optim_3, nll_pack_q1_3, opt_res_q1_3)
```
## Here we have:
## 'nll_expr_q1_4', 'nll_pack_q1_4', 'opt_res_q1_4' as the 
## expression, 'deriv_pack' and 'fit_optim' for the weibull distribution 
## model with median of exp(alpha4 + beta4*x) and
## standard deviation of sigma4
## alpha4 is written as 'alpha', beta4 is written as 'beta' and sigma4 is
## written as 'theta'

```{r P2_optim_4, message=FALSE, warning=FALSE}
nll_expr_q1_4 = expr(-log(log(2)) - theta + exp(theta) * (alpha + beta * x) + (1 - exp(theta)) * log(y) + log(2) * (y / exp(alpha + beta * x)) ^ exp(theta))
nll_pack_q1_4 = deriv_pack(nll_expr_q1_4)
opt_res_q1_4 = fit_optim(3, nll_pack_q1_4, dataQ1)
set_P2_optim(P2_optim_4, nll_pack_q1_4, opt_res_q1_4)
```

## Based on the outputs, we see that fitted model with the lowest AIC is 
## the logistic distribution with median of (alpha3+beta3*x) and scale
## parameter of sigma3. So according to the AIC, the best model is model 3.

## Question 1.2

```{r Q1.2}
n_grid     <- 1000 
x          <- seq(1, 15, length=n_grid)
P2_bands_1 <- P2_bands_2 <- 
  P2_bands_3 <- P2_bands_4 <- matrix(NA,nrow=n_grid,ncol=3)
colnames(P2_bands_1) <- colnames(P2_bands_2) <- 
  colnames(P2_bands_3) <- colnames(P2_bands_4) <- c("lower","est","upper")
head(P2_bands_1)
```
## Here we create a function 'set_P2_bands' which:
## for a data frame "P2_bands" of 3 columns and "opt_res"
## computes the lower quantile, mean and upper quantile of a 95% 
## confidence interval  in the first, second and third columns
## respectiveley for a given value of x

```{r set P2_bands}
set_P2_bands = function(P2_bands, opt_res) {
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
## Here we create function 'plot_bands' which plots the mean and 95% 
## confidence intervals of the function

```{r plot_P2_bands}
plot_bands = function(bands, data, main) {
  plot(data, pch = 19, main = main)
  lines(x, bands[, 2], col = "red")
  lines(x, bands[, 1], col = "red", lty = 2)
  lines(x, bands[, 3], col = "red", lty = 2)
}
```
## The 95% confidence inetrval for the normal distribution with
## median(or mean)=(alpha1 + beta1*x) and 
## standard deviation=sigma1*exp(alpha1+beta1*x)

```{r P2_bands_1}
P2_bands_1 = set_P2_bands(P2_bands_1, opt_res_q1_1)
plot_bands(P2_bands_1, dataQ1, 
           main = paste(c("Plot of the median and its 95% confidence",
                          "interval (CI) for model 1 against x")))
```
## The 95% confidence inetrval for the normal distribution with
## median(or mean)=(alpha2 + beta2*x) and 
## standard deviation=sigma2

```{r P2_bands_2}
P2_bands_2 = set_P2_bands(P2_bands_2, opt_res_q1_2)
plot_bands(P2_bands_2, dataQ1, 
           main = "Plot of the median for model 2 against x")
```
## The 95% confidence inetrval for the logistic distribution with
## median(or mean)=(alpha3 + beta3*x) and 
## scale parameter=sigma3

```{r P2_bands_3}
P2_bands_3 = set_P2_bands(P2_bands_3, opt_res_q1_3)
plot_bands(P2_bands_3, dataQ1,
           main = "Plot of the median for model 3 against x")
```
## The 95% confidence interval for the weibull distribution with
## median=exp(alpha4 + beta4*x) and 
## scale parameter=sigma4

```{r P2_bands_4}
P2_bands_4 = set_P2_bands(P2_bands_4, opt_res_q1_4)
plot_bands(P2_bands_4, dataQ1,
           main = "Plot of the median for model 4 against x")
```
# Part II

```{r}
#run to save the data frame 'dataQ2'  on to the R environment
load(url("https://people.bath.ac.uk/kai21/ASI/data/CW24/dataQ2.RData"))
```



## Question 2.1

```{r Q2.1_setup}
P21_optim_1 <- P21_optim_2 <- P21_optim_3 <- P21_optim_full <- P25_optim <- vector("list", 5)
names(P21_optim_1) <- names(P21_optim_2) <- names(P21_optim_3) <- names(P21_optim_full) <- names(P25_optim) <- c("mle","negloglik","gradient","hessian","NIC")

P21_optim_1
```
## We create function 'nic' which for given "nll_pack", "opt_res" and
## "data" computes the NIC value for the criterion

```{r nic}
nic = function(nll_pack, opt_res, data = dataQ2) {
  aux = nll_pack$gr(opt_res$par, data, apply.sum = F)
  2 * (opt_res$value + sum(diag(solve(opt_res$hessian) %*% crossprod(aux))))
}
```
## For given values of the mean(mu) and variance=theta*(mu^2) to compute
## the shape parameter(alpha) and rate parameter(lambda) of the gamma
## distribution we have that the mean=alpha/lambda and 
## variance=alpha/(lamda^2) for a gamma distribution. Substituting the 
## given values and rearraging we obtain alpha=1/theta and 
## lambda=1/(theta*mu). This will be applied for the expression of each of
## the negative log likelihood of the gamma distribution model in each case.

For $\mathcal{F}_{\text{full}}$, $\phi$ is reparameterised as $e^{\theta}$.

## Here we compute the MLE, negative log-likelhood, hessian, gradient and
## NIC values for the fitted full Gamma distribution with 
## mean=(1/alpha-beta*x) and variance=phi*(mu^2)
## 'phi' is written as 'theta' here

```{r P21_optim_full, message=FALSE, warning=FALSE}
nll_expr_q2_full = expr(lgamma(exp(-theta)) - 
                          exp(-theta) * (log(alpha + beta * x) - theta) + 
                          (1 - exp(-theta)) * log(y) + 
                          exp(-theta) * (alpha + beta * x) * y)
nll_pack_q2_full = deriv_pack(nll_expr_q2_full)
opt_res_q2_full = fit_optim(3, nll_pack_q2_full, dataQ2)
P21_optim_full$mle = c(opt_res_q2_full$par[1:2], exp(opt_res_q2_full$par[3]))
P21_optim_full$negloglik = opt_res_q2_full$value
P21_optim_full$hessian = opt_res_q2_full$hessian
P21_optim_full$gradient = nll_pack_q2_full$gr(opt_res_q2_full$par, dataQ2)
P21_optim_full$NIC = nic(nll_pack_q2_full, opt_res_q2_full)
P21_optim_full
```
## Here we compute the MLE, negative log-likelhood, hessian, gradient and
## NIC values for the fitted first Gamma distribution with 
## mean=(1/theta1-(theta1^2)*x) and variance=exp(theta1)*(mu^2)
```{r P21_optim_1, message=FALSE, warning=FALSE}
nll_expr_q2_1 = expr(lgamma(exp(-theta1)) - 
                       exp(-theta1) * (log(theta1 - x * theta1 ^ 2) - theta1) +
                       (1 - exp(-theta1)) * log(y) + 
                       exp(-theta1) * (theta1 - x * theta1 ^ 2) * y)
nll_pack_q2_1 = deriv_pack(nll_expr_q2_1,
                           c("theta1"),
                           theta_init = 2)
opt_res_q2_1 = fit_optim(1, nll_pack_q2_1, dataQ2)
P21_optim_1$mle = opt_res_q2_1$par
P21_optim_1$negloglik = opt_res_q2_1$value
P21_optim_1$hessian = opt_res_q2_1$hessian
P21_optim_1$gradient = nll_pack_q2_1$gr(opt_res_q2_1$par, dataQ2)
P21_optim_1$NIC = nic(nll_pack_q2_1, opt_res_q2_1)
P21_optim_1
```
## Here we compute the MLE, negative log-likelihood, hessian, gradient and
## NIC values for the fitted second Gamma distribution with 
## mean=(1/theta2-tau2*x) and variance=exp(theta2)*(mu^2)
```{r P21_optim_2}
nll_expr_q2_2 = expr(lgamma(exp(-theta2)) - 
                       exp(-theta2) * (log(theta2 + tau2 * x) - theta2) + 
                       (1 - exp(-theta2)) * log(y) + 
                       exp(-theta2) * (theta2 + tau2 * x) * y)
nll_pack_q2_2 = deriv_pack(nll_expr_q2_2,
                           c("theta2", "tau2"),
                           theta_init = rep(1, 2))
opt_res_q2_2 = fit_optim(2, nll_pack_q2_2, dataQ2)
P21_optim_2$mle = opt_res_q2_2$par
P21_optim_2$negloglik = opt_res_q2_2$value
P21_optim_2$hessian = opt_res_q2_2$hessian
P21_optim_2$gradient = nll_pack_q2_2$gr(opt_res_q2_2$par, dataQ2)
P21_optim_2$NIC = nic(nll_pack_q2_2, opt_res_q2_2)
P21_optim_2
```

For $\mathcal{F}_{3}$, $\delta_{3}$ is reparametrised as $\exp(\theta_{4})$.
## Here we compute the MLE, negative log-likelihood, hessian, gradient and
## NIC values for the fitted third Gamma distribution with 
## mean=(1/theta3-(tau3^2)*x) and variance=delta3*(mu^2)
## delta3 is written as "theta4"

```{r P21_optim_3}
nll_expr_q2_3 = expr(lgamma(exp(-theta4)) - 
                       exp(-theta4) * (log(theta3 - x * theta3 ^ 2) - theta4) + 
                       (1 - exp(-theta4)) * log(y) + 
                       exp(-theta4) * (theta3 - x * theta3 ^ 2) * y)
nll_pack_q2_3 = deriv_pack(nll_expr_q2_3,
                           c("theta3", "theta4"),
                           theta_init = rep(2, 2))
opt_res_q2_3 = fit_optim(2, nll_pack_q2_3, dataQ2)
P21_optim_3$mle = c(opt_res_q2_3$par[1], exp(opt_res_q2_3$par[2]))
P21_optim_3$negloglik = opt_res_q2_3$value
P21_optim_3$hessian = opt_res_q2_3$hessian
P21_optim_3$gradient = nll_pack_q2_3$gr(opt_res_q2_3$par, dataQ2)
P21_optim_3$NIC = nic(nll_pack_q2_3, opt_res_q2_3)
P21_optim_3
```
## In computing the NIC of each model, we find the 'full' model has the
## lowest value so according to the NIC, it is the best model.
$\mathcal{F}_{\text{full}}$ is the best based on the NIC.

## Question 2.2 

```{r Q2.2}
n_grid     <- 1000 
x          <- seq(0, 1, length=n_grid)
P22_bands_1 <- P22_bands_2 <- 
  P22_bands_3 <- P22_bands_full <- matrix(NA,nrow=n_grid,ncol=3)
colnames(P22_bands_1) <- colnames(P22_bands_2) <- 
  colnames(P22_bands_3) <- colnames(P22_bands_full) <- c("lower","est","upper")

head(P22_bands_1)
```
## Here we store the inverse hessian of the optimsiation for each model
```{r P22_bands}
S_full <- solve(opt_res_q2_full$hessian)
S_1 <- solve(opt_res_q2_1$hessian) 
S_2 <- solve(opt_res_q2_2$hessian)
S_3 <- solve(opt_res_q2_3$hessian)


## Here we find the variance for the function mu=(1/alpha+beta*x) 
## in each case through the jacobian and inverse hessian. Since the 
## model is correctly specified we can use the formula J*I*J to find 
## the asymptotic variance where 'J' is jacobian and 'I' is the inverse 
## hessian.

## For the full model the Jacobian is:
## (-((1/(alpha+beta*x))^2), (-x*((1/(alpha+beta*x))^2),0)
```{r mean_confint}
for (i in 1:n_grid){
  vec.x <- c(1, x[i])
  est_full <- crossprod(vec.x, opt_res_q2_full$par[1:2])
  J_full <- c(-(est_full^(-2)), -(est_full^(-2))*vec.x[2],0)
  se_full <- sqrt(J_full %*% S_full %*% J_full)
  P22_bands_full[i,1] <- (1/(est_full)) + 1.96 * se_full
  P22_bands_full[i,2] <- 1/est_full
  P22_bands_full[i,3] <- (1/(est_full)) - 1.96 * se_full
}

plot(y~x, dataQ2, ylim=c(0,10), main="Gamma model with (1/(alpha-beta1*x))")
lines(x, P22_bands_full[,1], col="red")
lines(x, P22_bands_full[,2], col="green")
lines(x, P22_bands_full[,3], col="blue")
legend("topleft", legend=c("Upper bound", "Mean", "Lower bound"), lty=c(1,1,1), col=c("red","green","blue"), cex=0.5)

## For the first model the Jacobian is:
## (-((1/(theta1-(theta1^2)*x))^2)*(1-2*theta1*x))
for (i in 1:n_grid){
  vec.x <- c(1, x[i])
  est_1 <- crossprod(vec.x, c(opt_res_q2_1$par, -(opt_res_q2_1$par^2)))
  J_1 <- c(-(est_1^(-2))*(1-2*opt_res_q2_1$par*vec.x[2]))
  se_1 <- sqrt(J_1 %*% S_1 %*% J_1)
  P22_bands_1[i,1] <- (1/(est_1)) + 1.96 * se_1
  P22_bands_1[i,2] <- 1/est_1
  P22_bands_1[i,3] <- (1/(est_1)) - 1.96 * se_1
}

plot(y~x, dataQ2, ylim=c(0,10), main="Gamma model with (1/(theta1-(theta1^2*x))")
lines(x, P22_bands_1[,1], col="red")
lines(x, P22_bands_1[,2], col="green")
lines(x, P22_bands_1[,3], col="blue")
legend("topleft", legend=c("Upper bound", "Mean", "Lower bound"), lty=c(1,1,1), col=c("red","green","blue"), cex=0.5)

## For the second model the Jacobian is:
## (-((1/(theta2+tau2*x))^2), (-x*((1/(theta2+tau2*x))^2))
for (i in 1:n_grid){
  vec.x <- c(1, x[i])
  est_2 <- crossprod(vec.x, opt_res_q2_2$par[1:2])
  J_2 <- c(-(est_2^(-2)), -(est_2^(-2))*vec.x[2])
  se_2 <- sqrt(J_2 %*% S_2 %*% J_2)
  P22_bands_2[i,1] <- (1/(est_2)) + 1.96 * se_2
  P22_bands_2[i,2] <- 1/est_2
  P22_bands_2[i,3] <- (1/(est_2)) - 1.96 * se_2
}
plot(y~x, dataQ2, ylim=c(0,10), main="Gamma model with (1/(theta2-tau2*x))")
lines(x, P22_bands_2[,1], col="red")
lines(x, P22_bands_2[,2], col="green")
lines(x, P22_bands_2[,3], col="blue")
legend("topleft", legend=c("Upper bound", "Mean", "Lower bound"), lty=c(1,1,1), col=c("red","green","blue"), cex=0.5)

## For the third model the Jacobian is:
## (-((1/(theta3-(theta3^2)*x))^2)*(1-2*theta3*x),0)
for (i in 1:n_grid){
  vec.x <- c(1, x[i])
  est_3 <- crossprod(vec.x, c(opt_res_q2_3$par[1], -(opt_res_q2_3$par[1]^2)))
  J_3 <- c(-(est_3^(-2)), -(est_3^(-2))*(1-2*opt_res_q2_3$par[1]*vec.x[2]))
  se_3 <- sqrt(J_3 %*% S_3 %*% J_3)
  P22_bands_3[i,1] <- (1/(est_3)) + 1.96 * se_3
  P22_bands_3[i,2] <- 1/est_3
  P22_bands_3[i,3] <- (1/(est_3)) - 1.96 * se_3
}

plot(y~x, dataQ2, ylim=c(0,10), main="Gamma model with (1/(theta3-(theta3^2*x))")
lines(x, P23_bands_2[,1], col="red")
lines(x, P23_bands_2[,2], col="green")
lines(x, P23_bands_2[,3], col="blue")
legend("topleft", legend=c("Upper bound", "Mean", "Lower bound"), lty=c(1,1,1), col=c("red","green","blue"), cex=0.5)


## Question 2.3


```{r Q2.3}
n_grid <- 1000 
x <- seq(0,1,length=n_grid)
P23_bands_1 <- P23_bands_2 <- 
  P23_bands_3 <- P23_bands_full <- P25_bands <- matrix(NA,nrow=n_grid,ncol=3)
colnames(P23_bands_1) <- colnames(P23_bands_2) <- 
  colnames(P23_bands_3) <- colnames(P23_bands_full) <- colnames(P25_bands) <- c("lower","est","upper")

head(P23_bands_1)
```
## To compute the variance of the incorrectly specified distribution,
## we use proposition 3.7 with the fact that the least worst theta estimate
## is unknown so k-hat is found by the sum of the crossproduct of the
## gradient of log densities and J is the fisher information matrix.
```{r P23_bands_1}
S_full = solve(opt_res_q2_full$hessian)
nll_gr_q2_full = nll_pack_q2_full$gr(opt_res_q2_full$par, dataQ2, apply.sum = F)
Khat_full = crossprod(nll_gr_q2_full)
cov_full = S_full %*% Khat_full %*% S_full

for (i in 1:n_grid){
  vec.x <- c(1, x[i])
  est_full <- crossprod(vec.x, opt_res_q2_full$par[1:2])
  J_full <- c(-(est_full^(-2)), -(est_full^(-2))*vec.x[2],0)
  se_full <- sqrt(J_full %*% cov_full %*% J_full)
  P23_bands_full[i,1] <- (1/(est_full)) + 1.96 * se_full
  P23_bands_full[i,2] <- 1/est_full
  P23_bands_full[i,3] <  (1/(est_full)) - 1.96 * se_full
}

plot(y~x, dataQ2, ylim=c(0,10), main="Mispecified Gamma model with (1/(alpha-beta*x)")
lines(x, P23_bands_full[,1], col="red")
lines(x, P23_bands_full[,2], col="green")
lines(x, P23_bands_full[,3], col="blue")
legend("topleft", legend=c("Upper bound", "Mean", "Lower bound"), lty=c(1,1,1), col=c("red","green","blue"), cex=0.5)

## First model
nll_gr_q2_1 = nll_pack_q2_1$gr(opt_res_q2_1$par, dataQ2, apply.sum = F)
Khat_1 = crossprod(nll_gr_q2_1)
cov_1 = S_1 %*% Khat_1 %*% S_1
for (i in 1:n_grid){
  vec.x <- c(1, x[i])
  est_1 <- crossprod(vec.x, c(opt_res_q2_1$par, -(opt_res_q2_1$par^2)))
  J_1 <- c(-(est_1^(-2))*(1-2*opt_res_q2_1$par*vec.x[2]))
  se_1 <- sqrt(J_1 %*% cov_1 %*% J_1)
  P23_bands_1[i,1] <- (1/(est_1)) + 1.96 * se_1
  P23_bands_1[i,2] <- 1/est_1
  P23_bands_1[i,3] <- (1/(est_1)) - 1.96 * se_1
}
plot(y~x, dataQ2, ylim=c(0,10), main="Gamma model with (1/(theta1-(theta1^2)*x))")
lines(x, P23_bands_1[,1], col="red")
lines(x, P23_bands_1[,2], col="green")
lines(x, P23_bands_1[,3], col="blue")
legend("topleft", legend=c("Upper bound", "Mean", "Lower bound"), lty=c(1,1,1), col=c("red","green","blue"), cex=0.5)

## Second model
nll_gr_q2_2 = nll_pack_q2_2$gr(opt_res_q2_2$par, dataQ2, apply.sum = F)
Khat_2 = crossprod(nll_gr_q2_2)
cov_2 = S_2 %*% Khat_2 %*% S_2
for (i in 1:n_grid){
  vec.x <- c(1, x[i])
  est_2 <- crossprod(vec.x, opt_res_q2_2$par[1:2])
  J_2 <- c(-(est_2^(-2)), -(est_2^(-2))*vec.x[2])
  se_2 <- sqrt(J_2 %*% cov_2 %*% J_2)
  P23_bands_2[i,1] <- (1/(est_2)) + 1.96 * se_2
  P23_bands_2[i,2] <- 1/est_2
  P23_bands_2[i,3] <- (1/(est_2)) - 1.96 * se_2
}

plot(y~x, dataQ2, ylim=c(0,10), main="Gamma model with (1/(theta2-tau2*x))")
lines(x, P23_bands_1[,1], col="red")
lines(x, P23_bands_1[,2], col="green")
lines(x, P23_bands_1[,3], col="blue")
legend("topleft", legend=c("Upper bound", "Mean", "Lower bound"), lty=c(1,1,1), col=c("red","green","blue"), cex=0.5)

## Third model
nll_gr_q2_3 = nll_pack_q2_3$gr(opt_res_q2_3$par, dataQ2, apply.sum = F)
Khat_3 = crossprod(nll_gr_q2_3)
cov_3 = S_3 %*% Khat_3 %*% S_3
for (i in 1:n_grid){
  vec.x <- c(1, x[i])
  est_3 <- crossprod(vec.x, c(opt_res_q2_3$par[1], -(opt_res_q2_3$par[1]^2)))
  J_3 <- c(-(est_3^(-2)), -(est_3^(-2))*(1-2*opt_res_q2_3$par[1]*vec.x[2]))
  se_3 <- sqrt(J_3 %*% cov_3 %*% J_3)
  P23_bands_3[i,1] <- (1/(est_3)) + 1.96 * se_3
  P23_bands_3[i,2] <- 1/est_3
  P23_bands_3[i,3] <- (1/(est_3)) - 1.96 * se_3
}

plot(y~x, dataQ2, ylim=c(0,10), main=" Mispecified Gamma model with (1/(theta3-(theta3^2)*x))")
lines(x, P23_bands_3[,1], col="red")
lines(x, P23_bands_3[,2], col="green")
lines(x, P23_bands_3[,3], col="blue")
legend("topleft", legend=c("Upper bound", "Mean", "Lower bound"), lty=c(1,1,1), col=c("red","green","blue"), cex=0.5)

```

## Question 2.4

With $\alpha^{*} = \frac{2}{3}$, $\beta^{*} = -\frac{1}{3}$, and $\phi^{*} = \frac{1}{2}$, the correct shape and scale parameters are $\frac{1}{\phi^{*}} = 2$ and $\frac{\alpha^{*} + \beta^{*}x}{\phi^{*}} = \frac{4 - 2x}{3}$.

```{r Q2.4}
n_grid     <- 1000 
x          <- seq(0, 1,length=n_grid)
P24_bands_1 <- P24_bands_2 <- 
  P24_bands_3 <- P24_bands_full <- matrix(NA,nrow=n_grid,ncol=3)
colnames(P24_bands_1) <- colnames(P24_bands_2) <- 
  colnames(P24_bands_3) <- colnames(P24_bands_full) <- c("lower","est","upper")

head(P24_bands_1)
```
## Here, we create a function 'KL_pack' which computes the value of the
## KL function and its derivative for given values of "nll_pack" with
## the correctly specified full model as the true distribution
```{r KL_pack}
KL_pack = function(nll_pack) {
  fn_integrand = function(y, x, theta) {
    (nll_pack$fn(theta, list("x" = x, "y" = y), apply.sum = F) +
       dgamma(y, 2, (4 - 2 * x)/3, log = T)) * 
      dgamma(y, 2, (4 - 2 * x)/3)
  }
  gr_integrand = function(y, x, theta, col_int = 1) {
    aux = nll_pack$gr(theta, list("x" = x, "y" = y), apply.sum = F) *
      dgamma(y, 2, (4 - 2 * x)/3)
    aux[, col_int]
  }
  KL_fn = function(theta, data_list = dataQ2) {
    # Although unused, the data_list argument is added for compatibility with fit_optim.
    res = 1:10
    for(i in res) {
      res[i] = integrate(fn_integrand, 0, 10, x = i/10, 
                         theta = theta)$value
    }
    sum(res)
  }
  
  KL_gr = function(theta, data_list = dataQ2) {
    # Although unused, the data_list argument is added for compatibility with fit_optim.
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

##First model
```{r KL1}
KL_pack_1 = KL_pack(nll_pack_q2_1)
opt_res_q2_kl1 = fit_optim(1, KL_pack_1, dataQ2, mean = 0.5, sd = 0.25, silent = F)
z <- rep(0, n_grid)
for(i in 1:n_grid){
  vec.x <- c(1,x[i])
  theta <- opt_res_q2_kl1$par
  z1 <- crossprod(vec.x, c(theta,-(theta^2)))
  z[i] <- 1/z1
}
plot(y~x, dataQ2, ylim=c(0,10), main=" Mispecified Gamma first model - Q2_4 )")
lines(x, z, col="red")
```
## In each case we pass the values of 'KL-pack' through 'fit_optim'
## to obtain the least worst value of theta


## Second model
```{r KL2}
KL_pack_2 = KL_pack(nll_pack_q2_2)
opt_res_q2_kl2 = fit_optim(2, KL_pack_2, dataQ2, mean = 0.5, sd = 0.25, silent = F)
z <- rep(0, n_grid)
for(i in 1:n_grid){
  vec.x <- c(1,x[i])
  theta <- opt_res_q2_kl2$par
  z1 <- crossprod(vec.x, theta)
  z[i] <- 1/z1
}
plot(y~x, dataQ2, ylim=c(0,10), main=" Mispecified Gamma second model - Q2_4 )")
lines(x, z, col="red")
       ```
```
## Third model
```{r KL3}
KL_pack_3 = KL_pack(nll_pack_q2_3)
opt_res_q2_kl3 = fit_optim(2, KL_pack_3, dataQ2, mean = 0.5, sd = 0.25, silent = F)
z <- rep(0, n_grid)
for(i in 1:n_grid){
  vec.x <- c(1,x[i])
  theta <- opt_res_q2_kl3$par[1]
  z1 <- crossprod(vec.x, c(theta,-(theta^2)))
  z[i] <- 1/z1
}
plot(y~x, dataQ2, ylim=c(0,10), main=" Mispecified Gamma third model - Q2_4 )")
lines(x, z, col="red")
       ```
save.image()
```

## Question 2.5

From [Question 2.1], $\mathcal{F}_{\text{full}}$ had the best NIC of
$-1368$. Finding functions $g_{\alpha}$, $g_{\beta}$, and $g_{\phi}$ such that
$g_{\alpha}(\theta_{4}^{*}) = \alpha^{*}$, 
$g_{\beta}(\theta_{4}^{*}) = \beta^{*}$, and 
$g_{\phi}(\theta_{4}^{*}) = \phi^{*}$ would give $\mathcal{F}_{4}$ a better NIC
since the NIC would favour the simpler model $\mathcal{F}_{4}$.

$\mathcal{F}_{2}$ had the next lowest AIC of $-1359$, followed by 
$\mathcal{F}_{3}$ with an NIC of $-1224$. As such, we considered
$g_{alpha}(\theta_{4}) = \theta_{4}$, 
$g_{\beta}(\theta_{4}) = -\theta_{4}^{a}$, 
$g_{\phi}(\theta_{4}) = \exp(\theta_{4}^{b})$, where $a$ and $b$ are real
constants. By solving for the conditions mentioned in the previous paragraph, 
we found that $a \approx 0.75$ and $b \approx 0.5$. 

```{r}
a = log(-P21_optim_full$mle[2])/log(P21_optim_full$mle[1])
b = log(log(P21_optim_full$mle[3]))/log(P21_optim_full$mle[1])
print(paste("a: ", a, "; b: ", b))
```

With $g_{\beta}(\theta_{4}) = -\theta_{4}^{3/4}$ and
$g_{\phi}(\theta_{4}) = \exp(\theta_{4}^{1/2})$, $\mathcal{F}_{4}$ achieved
the lowest NIC of $-1372$.

```{r P25_optim}
nll_expr_q2_4 = expr(lgamma(exp(-theta4 ^ 0.5)) -
                       exp(-theta4 ^ 0.5) * (log(theta4 - x * theta4 ^ 0.75) - theta4 ^ 0.5) +
                       (1 - exp(-theta4 ^ 0.5)) * log(y) +
                       exp(-theta4 ^ 0.5) * (theta4 - x * theta4 ^ 0.75) * y)

nll_pack_q2_4 = deriv_pack(nll_expr_q2_4,
                           c("theta4"),
                           theta_init = 0.5)
opt_res_q2_4 = fit_optim(1, nll_pack_q2_4, dataQ2, sd = 2, silent = F)
P25_optim$mle = opt_res_q2_4$par
P25_optim$negloglik = opt_res_q2_4$value
P25_optim$hessian = opt_res_q2_4$hessian
P25_optim$gradient = nll_pack_q2_4$gr(opt_res_q2_4$par, dataQ2)
P25_optim$NIC = nic(nll_pack_q2_4, opt_res_q2_4)
P25_optim
```

```{r P25_bands}
P25_bands = mean_confint(P25_bands, 
                         expr(1/(theta4 - x * theta4 ^ 0.75)),
                         c("theta4"), 
                         opt_res_q2_4)
plot_bands(P25_bands, dataQ2, main = "Mean and 95% CI for model 4")
```

# References